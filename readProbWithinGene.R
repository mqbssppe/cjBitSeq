library('foreach')
library('doMC')


filesA <- list.files(pattern="conditionA")
filesB <- list.files(pattern="conditionB")

dir.create("tmp")
setwd("tmp")
for (sth in filesA){
	file.symlink(from = paste("../",sth,sep=""), to = sth)
}
for (sth in filesB){
	file.symlink(from = paste("../",sth,sep=""), to = sth)
}



fls<- c("A_1","B_1")
nreps1 <- 1
nreps2 <- 1

nParallelJobs <- length(filesA)+length(filesB)
nFiles <- nParallelJobs
nParallelJobs <- as.numeric(read.table("../nCores.txt"))
registerDoMC(nParallelJobs)
cat(c("number of cores: ",nParallelJobs),"\n")



n1 <- length(filesA)
n2 <- length(filesB)

nReplicate1 <- n1
nReplicate2 <- n2
for (j in 1:n1){
	system(paste("sed -i 1,5d", filesA[j]))
}
for (j in 1:n2){
	system(paste("sed -i 1,5d", filesB[j]))
}
system(paste("cat", paste(filesA,collapse=" "),"> mergedfilesA"))
system(paste("cat", paste(filesB,collapse=" "),"> mergedfilesB"))
nreadsA <- as.numeric(strsplit(system("wc -l mergedfilesA",intern=TRUE),split=" ")[[1]][1])
nLinesPerFile <- floor(nreadsA/(nParallelJobs/2) + 1)
system(paste("split --numeric-suffixes -l", nLinesPerFile, "mergedfilesA condA_"))
nreadsB <- as.numeric(strsplit(system("wc -l mergedfilesB",intern=TRUE),split=" ")[[1]][1])
nLinesPerFile <- floor(nreadsB/(nParallelJobs/2) + 1)
system(paste("split --numeric-suffixes -l", nLinesPerFile, "mergedfilesB condB_"))
for (i in filesA){
	system(paste("rm",i))
}
for (i in filesB){
	system(paste("rm",i))
}
filesA <- list.files(pattern = "condA_")
filesB <- list.files(pattern = "condB_")
nFiles <- length(filesA) + length(filesB)

# read geneToTranscript file
geneToTranscripts <-  read.table("../genesTranscripts.txt",header=TRUE,stringsAsFactors=FALSE)
#geneToTranscripts <-  read.table("../genesTranscriptsMapping.txt",header=TRUE,stringsAsFactors=FALSE)
missingGeneNames <- which(is.na(geneToTranscripts$geneID)==TRUE)
if(length(missingGeneNames)>0){
	geneToTranscripts$geneID[missingGeneNames] <- paste('MISSINGNAME',1:length(missingGeneNames),sep="")
}
clusterNames <- as.character(unique(geneToTranscripts$geneID))
nReadsPerCluster <- numeric(length(clusterNames))
names(nReadsPerCluster) <- clusterNames


readProb <- function(fileName){
	con <- file(fileName,open = "r")
	outFileName <- paste(fileName,"TMP",sep="")
	outFileNameForNReads <- paste(fileName,"nReads",sep="")
	print(files[fileIndex])
	conOut = file(outFileName,open="w")

	weirdos <- 0
	badIndex <- 0
	pseudoreads <- 0
	multiMapReads <- 0
	nReadsPerCluster <- numeric(length(clusterNames))
	names(nReadsPerCluster) <- clusterNames
	i <- 0
#	for(i in 1:5){
#		line <- readLines(con, n = 1)
#	}
	nLines<-500 ## streams nLines lines per iteration
	totalReads <- as.numeric(strsplit(system(paste("wc -l", fileName),intern=TRUE),split=" ")[[1]][1]) - 5
	mainRound <- floor(totalReads/nLines)*nLines
	remainingLines <- totalReads - floor(totalReads/nLines)*nLines
	i<-0
	#while(length(line <- readLines(con, n = nLines, warn = FALSE)) > 0){
	while ((length(line <- readLines(con, n = nLines, warn = FALSE)) > 0)&(i < mainRound)) {
		i <- i + nLines
		outLines <- c()
		for (lineIterator in 1:nLines){
			ln1 <- as.numeric(strsplit(line[lineIterator], split=" ")[[1]][-1])
			l <- 2*ln1[1] - 2
			maps <- ln1[seq(2,l,by = 2)]
			nMaps <- length(maps)
			ln1 <- ln1[1:(l+1)]
			ln1[1] <- nMaps
			if(sum(maps)==0){pseudoreads <- pseudoreads + 1}else{
				namesUnique <- unique(maps)
				lU <- length(namesUnique)
				if(l/2 != lU){
					newLL <- numeric(lU) 
					for(k in 1:lU){
						newLL[k] <- log(sum(exp(ln1[2*which(maps == namesUnique[k])+1])))
					}
					newLL[lU] <- ln1[l+1]
					checkInf <- which(is.finite(newLL)==FALSE)
					lci <- length(checkInf)
					if (lci > 0){newLL[checkInf] <- rep(minThreshold,lci)}
					perm <- c(order(newLL[-lU],decreasing=TRUE),lU)
					ln1 <- c(lU,c(rbind(namesUnique[perm],newLL[perm])))
					l <- 2*ln1[1]
					maps <- ln1[seq(2,l,by = 2)]
					nMaps <- lU
					weirdos <- weirdos + 1
				}
				geneNames <- as.character(geneToTranscripts[maps[1:(nMaps)],]$geneID)
				uniqueGeneNames <- unique(geneNames)
				if(length(uniqueGeneNames) > 1){
					myMatch <- match(geneNames,uniqueGeneNames)
					myWeight <- as.numeric(table(myMatch))[myMatch]
					myProb <- exp(ln1[seq(3,l+1,by = 2)])
					myProb[is.na(myProb)] <- 0
					myProb[is.infinite(myProb)] <- 0
					myProb <- myProb/myWeight
					if (sum(myProb) == 0){myProb = rep(1,l/2)}
					mySample <- sample(nMaps,1,prob = myProb)
					randInd <- maps[mySample]
					subSet <- which(myMatch == myMatch[mySample] )
					ln1 <- c(length(subSet),ln1[c(rbind(2*subSet,2*subSet + 1))])
					maps <- maps[subSet]
					nMaps <- length(subSet) 
					l <- 2*nMaps
					badIndex <- badIndex + 1
					geneNames <- geneNames[mySample]
				}
				clusterID <- which(clusterNames == geneNames[1])
				nReadsPerCluster[clusterID] <- 1 + nReadsPerCluster[clusterID]
				myNewLine <- c(clusterID,ln1)
				if (geneToTranscripts[maps[1],4] > 1){ # number of transcripts for this gene. if it is equal to 1, then we do not consider it.
						multiMapReads <- multiMapReads + 1
						newNames <- geneToTranscripts[maps[1:nMaps],5] - 1 #  transcriptID's within gene
						myNewLine[1 + seq(2,l,by = 2)] <- newNames
						outLines <- rbind(outLines,paste(myNewLine,collapse = " "))

				}
			}
		}
		cat(outLines,sep="\n",file = conOut,append=TRUE)
	

		if (i%%100000 == 0){
			cat(c(fileName,": ",paste(c(i,round(100*badIndex/i,3), round(100*badIndex/multiMapReads,3)))),"\n")
		}
	}
	nLines <- length(line)

	i <- i + nLines
	outLines <- c()
	for (lineIterator in 1:nLines){
		ln1 <- as.numeric(strsplit(line[lineIterator], split=" ")[[1]][-1])
		l <- 2*ln1[1] - 2
		maps <- ln1[seq(2,l,by = 2)]
		nMaps <- length(maps)
		ln1 <- ln1[1:(l+1)]
		ln1[1] <- nMaps
		if(sum(maps)==0){pseudoreads <- pseudoreads + 1}else{
			namesUnique <- unique(maps)
			lU <- length(namesUnique)
			if(l/2 != lU){
				newLL <- numeric(lU) 
				for(k in 1:lU){
					newLL[k] <- log(sum(exp(ln1[2*which(maps == namesUnique[k])+1])))
				}
				newLL[lU] <- ln1[l+1]
				checkInf <- which(is.finite(newLL)==FALSE)
				lci <- length(checkInf)
				if (lci > 0){newLL[checkInf] <- rep(minThreshold,lci)}
				perm <- c(order(newLL[-lU],decreasing=TRUE),lU)
				ln1 <- c(lU,c(rbind(namesUnique[perm],newLL[perm])))
				l <- 2*ln1[1]
				maps <- ln1[seq(2,l,by = 2)]
				nMaps <- lU
				weirdos <- weirdos + 1
			}
			geneNames <- as.character(geneToTranscripts[maps[1:(nMaps)],]$geneID)
			uniqueGeneNames <- unique(geneNames)
			if(length(uniqueGeneNames) > 1){
				myMatch <- match(geneNames,uniqueGeneNames)
				myWeight <- as.numeric(table(myMatch))[myMatch]
				myProb <- exp(ln1[seq(3,l+1,by = 2)])
				myProb[is.na(myProb)] <- 0
				myProb[is.infinite(myProb)] <- 0
				myProb <- myProb/myWeight
				if (sum(myProb) == 0){myProb = rep(1,l/2)}
				mySample <- sample(nMaps,1,prob = myProb)
				randInd <- maps[mySample]
				subSet <- which(myMatch == myMatch[mySample] )
				ln1 <- c(length(subSet),ln1[c(rbind(2*subSet,2*subSet + 1))])
				maps <- maps[subSet]
				nMaps <- length(subSet) 
				l <- 2*nMaps
				badIndex <- badIndex + 1
				geneNames <- geneNames[mySample]
			}
			clusterID <- which(clusterNames == geneNames[1])
			nReadsPerCluster[clusterID] <- 1 + nReadsPerCluster[clusterID]
			myNewLine <- c(clusterID,ln1)
			if (geneToTranscripts[maps[1],4] > 1){ # number of transcripts for this gene. if it is equal to 1, then we do not consider it.
					multiMapReads <- multiMapReads + 1
					newNames <- geneToTranscripts[maps[1:nMaps],5] - 1 #  transcriptID's within gene
					myNewLine[1 + seq(2,l,by = 2)] <- newNames
					outLines <- rbind(outLines,paste(myNewLine,collapse = " "))

			}
		}
	}
	cat(outLines,sep="\n",file = conOut,append=TRUE)

	close(con)
	close(conOut)
	write.table(nReadsPerCluster,file = outFileNameForNReads,quote=FALSE)


}



readProbOld <- function(fileName){
	con <- file(fileName,open = "r")
	outFileName <- paste(fileName,"TMP",sep="")
	outFileNameForNReads <- paste(fileName,"nReads",sep="")
	print(files[fileIndex])
	conOut = file(outFileName,open="w")

	weirdos <- 0
	badIndex <- 0
	pseudoreads <- 0
	multiMapReads <- 0
	nReadsPerCluster <- numeric(length(clusterNames))
	names(nReadsPerCluster) <- clusterNames
	i <- 0
	for(i in 1:5){
		line <- readLines(con, n = 1)
	}

	while(length(line <- readLines(con, n = 1, warn = FALSE)) > 0){
	#while ((length(line <- readLines(con, n = 1, warn = FALSE)) > 0)&(i < 10000)) {
		i <- i + 1

		ln1 <- as.numeric(strsplit(line, split=" ")[[1]][-1])
		l <- 2*ln1[1] - 2
		maps <- ln1[seq(2,l,by = 2)]
		nMaps <- length(maps)
		ln1 <- ln1[1:(l+1)]
		ln1[1] <- nMaps
		if(sum(maps)==0){pseudoreads <- pseudoreads + 1}else{
			namesUnique <- unique(maps)
			lU <- length(namesUnique)
			if(l/2 != lU){
				newLL <- numeric(lU) 
				for(k in 1:lU){
					newLL[k] <- log(sum(exp(ln1[2*which(maps == namesUnique[k])+1])))
				}
				newLL[lU] <- ln1[l+1]
				checkInf <- which(is.finite(newLL)==FALSE)
				lci <- length(checkInf)
				if (lci > 0){newLL[checkInf] <- rep(minThreshold,lci)}
				perm <- c(order(newLL[-lU],decreasing=TRUE),lU)
				#ln1 <- c(lU,paste(namesUnique[perm],newLL[perm]))
				ln1 <- c(lU,c(rbind(namesUnique[perm],newLL[perm])))
				l <- 2*ln1[1]
				maps <- ln1[seq(2,l,by = 2)]
				nMaps <- lU
				#cat(c(clusterID,lU,paste(namesUnique[perm],newLL[perm])),"\n",file=conOut,append=TRUE)
				weirdos <- weirdos + 1
			}
			geneNames <- as.character(geneToTranscripts[maps[1:(nMaps)],]$geneID)
			uniqueGeneNames <- unique(geneNames)
			if(length(uniqueGeneNames) > 1){
				myMatch <- match(geneNames,uniqueGeneNames)
				myWeight <- as.numeric(table(myMatch))[myMatch]
				myProb <- exp(ln1[seq(3,l+1,by = 2)])
				myProb[is.na(myProb)] <- 0
				myProb[is.infinite(myProb)] <- 0
				myProb <- myProb/myWeight
				if (sum(myProb) == 0){myProb = rep(1,l/2)}
				mySample <- sample(nMaps,1,prob = myProb)
				randInd <- maps[mySample]
				subSet <- which(myMatch == myMatch[mySample] )
				#subSet <- which(abs(maps[-l/2] - randInd) < 10)
				ln1 <- c(length(subSet),ln1[c(rbind(2*subSet,2*subSet + 1))])
				maps <- maps[subSet]
				nMaps <- length(subSet) 
				l <- 2*nMaps
				badIndex <- badIndex + 1
				geneNames <- geneNames[mySample]
			}
			clusterID <- which(clusterNames == geneNames[1])
			nReadsPerCluster[clusterID] <- 1 + nReadsPerCluster[clusterID]
			line <- c(clusterID,ln1)
			if (geneToTranscripts[maps[1],4] > 1){ # number of transcripts for this gene. if it is equal to 1, then we do not consider it.
				multiMapReads <- multiMapReads + 1
				newNames <- geneToTranscripts[maps[1:nMaps],5] - 1 #  transcriptID's within gene
				line[1 + seq(2,l,by = 2)] <- newNames
				cat(line,"\n",file = conOut,append=TRUE)
			}
		}
		if (i%%100000 == 0){
			cat(c(fileName,": ",paste(c(i,round(100*badIndex/i,3), round(100*badIndex/multiMapReads,3)))),"\n")
		}
	}
	close(con)
	close(conOut)
	write.table(nReadsPerCluster,file = outFileNameForNReads,quote=FALSE)

}





files <- c(filesA,filesB)
fileIndex <- 1
foreach(fileIndex=1:nFiles) %dopar% {
	fileName <- files[fileIndex]
	readProb(fileName)		
}




TMPfilesA <- paste(filesA,"TMP",sep="")
n1 <- length(TMPfilesA)
myCommand <- paste(c('cat',TMPfilesA,"> conditionA_1.panos"),collapse=" ")
system(myCommand)
TMPfilesB <- paste(filesB,"TMP",sep="")
n2 <- length(TMPfilesB)
myCommand <- paste(c('cat',TMPfilesB,"> conditionB_1.panos"),collapse=" ")
system(myCommand)


TMPfilesA <- paste(filesA,'nReads',sep='')
nReadsA <- read.table(TMPfilesA[1])
if(n1 > 1){for(i in 2:n1){
myfile <- TMPfilesA[i]
nReadsA[,1] <- nReadsA[,1] + read.table(myfile)[,1]
}
}

TMPfilesB <- paste(filesB,'nReads',sep='')
nReadsB <- read.table(TMPfilesB[1])
if(n2 > 1){for(i in 2:n2){
myfile <- TMPfilesB[i]
nReadsB[,1] <- nReadsB[,1] + read.table(myfile)[,1]
}
}


if(1>2){system('rm *.probTMP *probnReads')}
nReads <- cbind(nReadsA,nReadsB)
colnames(nReads) <- c('A','B')
write.table(nReads,file = 'nReads.txt',quote=FALSE)


ind <- which( (nReads$B + nReads$A > 0)&(nReads$B*nReads$A == 0) )
ind <- match(rownames(nReads)[ind],geneToTranscripts$geneID)
ind1 <- which(geneToTranscripts[ind,4] > 1)
ind <- ind[ind1]
if(length(ind)>0){
	for (i in ind){
		nTr <- geneToTranscripts[i,4]
		clusterID <-  which(clusterNames == geneToTranscripts[i,1])
		myAddedSentence <- c(clusterID,nTr,c(rbind((nTr - 1):0,c(rep(-20,nTr)))))
		myAddedSentence <- paste(as.character(myAddedSentence),collapse=' ')
		myCommand <- paste('echo "',myAddedSentence,'" >> conditionA_1.panos',sep = '')
		system(myCommand)
		myCommand <- paste('echo "',myAddedSentence,'" >> conditionB_1.panos',sep = '')
		system(myCommand)
	}
}



multiGenes <- which(geneToTranscripts[,4] > 1)
multiGeneNames <- unique(geneToTranscripts[multiGenes,1])
# this will apply MCMC to all expressed genes
#expressedMultiGeneNames <- names(which(rowSums(nReads[multiGeneNames,]) > 0)) # names of genes
# this will apply MCMC to all expressed genes > thresholdCount = (nreps1 + nreps2)*10
expressedMultiGeneNames <- names(which(rowSums(nReads[multiGeneNames,]) > (nReplicate1+nReplicate2)*105)) # names of genes
expressedMultiGeneID <- match(expressedMultiGeneNames,clusterNames) # clusterID of genes

expressedMultiGeneNReads <- rowSums(nReads[expressedMultiGeneID,])
myPriorityRank <- order(expressedMultiGeneNReads,decreasing=TRUE)
nJobs <- length(expressedMultiGeneNReads)
myJobMatrix <- cbind(1:nJobs,expressedMultiGeneID[myPriorityRank], expressedMultiGeneNames[myPriorityRank], expressedMultiGeneNReads[myPriorityRank])
colnames(myJobMatrix) <- c('jobID','geneID','geneNAME','nReads')


dir.create("jobs")
setwd("jobs/")
k=1
foreach(k=1:nJobs) %dopar% {
        job <- file(paste("job_",k,".sh",sep = ""),open = "w")
        cat("#!/bin/bash","\n",file = job)
#        cat(paste("cd ../clusters/cluster_",myJobMatrix[k,'geneID'],sep=""),"\n",file = job)
        cat(paste("cd ../clusters/cluster_",k,sep=""),"\n",file = job)
        cat("head -1 data.tr","\n",file = job)
	cat("withinGeneCollapsedSampler > sampler.log","\n",file = job)  # collapsed sampler
        close(job)
        system(paste("chmod u+x job_",k,".sh",sep=""))       
}
setwd('../')

dir.create("clusters")
setwd("clusters/")
k=1
foreach(k=1:nJobs) %dopar% {
	myDir <- paste("cluster_",k,sep="")
	dir.create(myDir)
	setwd(myDir)
	transcriptsOfGene <- which(geneToTranscripts$geneID == myJobMatrix[k,'geneNAME'])
	write.table(geneToTranscripts[transcriptsOfGene,],file = 'names.info',quote=FALSE)
	nComponents <- geneToTranscripts[transcriptsOfGene,4][1]
	cat(paste("# ",0," ",nComponents," ",myJobMatrix[k,'geneID']," ",myJobMatrix[k,'nReads'],sep = ""),file="data.tr",sep="\n")
	setwd('../')
}
setwd('../')

job <- file(paste("parallelGNUwithin",".bash",sep = ""),open = "w")
cat("#!/bin/bash","\n",file = job)
cat("echo \"(#M nTranscripts GeneID nReads)\"","\n",file = job)
cat("cd jobs","\n",file = job)
cat(paste("parallel --gnu ./job_{}.sh ::: {1..",nJobs,"}",sep=""),"\n",file = job)
cat("wait","\n",file = job)
cat("cd ..","\n",file = job)
cat("echo \"Processing output...\"","\n",file = job)
cat("R CMD BATCH ../getClustersWithin.R","\n",file = job)
close(job)
system(paste("chmod u+x parallelGNUwithin",".bash",sep=""))





