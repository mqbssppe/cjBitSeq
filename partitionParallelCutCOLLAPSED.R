library('Matrix')
library('foreach')
library('doMC')
registerDoMC(2)

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



#K <- 28764
maxK1 <- maxK2 <- 0
fls<- c("A_1","B_1")
nreps1 <- 1
nreps2 <- 1
nCon <- 110

filesA <- list.files(pattern="conditionA")
filesB <- list.files(pattern="conditionB")

maxDistance <- 500
minThreshold <- -720

foreach(fileIndex=1:2) %dopar% {

	if(fileIndex == 1){
		conOut = file(paste("conditionA_1New.prob",sep=""),open="w")
		clustOut1 <- file(paste("A_1First.temp",sep=""),open = "w")
		iter <- 1
		for(outFile in filesA){
		cat(paste("correcting same transcript alignments at file ", outFile,".",sep=""),sep="\n")
		con = file(outFile)
		open(con)
		if(iter == 1){
		for(i in 1:4){
			line <- readLines(con, n = 1)
			cat(line,"\n",file = conOut,append=TRUE)
		}}else{
		for(i in 1:4){
			line <- readLines(con, n = 1)
		}
		}
		i <- 4
		weirdos <- 0
		badIndex <- 0
		pseudoreads <- 0
		while (length(line <- readLines(con, n = 1, warn = FALSE)) > 0) {
		i <- i + 1

		ln1 <- as.numeric(strsplit(line, split=" ")[[1]][-1])
		l <- 2*ln1[1]
		maps <- ln1[seq(2,l,by = 2)]
		if(sum(maps)==0){pseudoreads <- pseudoreads + 1}else{
			myRange <- diff(range(maps[-l/2]))
			if(myRange > maxDistance){
				myProb <- exp(ln1[seq(3,l,by = 2)])
				myProb[is.na(myProb)] <- 0
				myProb[is.infinite(myProb)] <- 0
				if (sum(myProb) == 0){myProb = rep(1,l/2 - 1)}
				randInd <- maps[sample(l/2 - 1,1,prob = myProb)]
				subSet <- which(abs(maps[-l/2] - randInd) < 10)
				ln1 <- c(length(subSet) + 1,ln1[c(rbind(2*subSet,2*subSet + 1),l,l+1)])
				maps <- c(maps[-l/2][subSet],0)
				l <- 2*(length(subSet) + 1)
				badIndex <- badIndex + 1
				line <- c(NA, ln1)
			}

			cat(ln1[2],"\n",file=clustOut1)
			maxMaps <- max(maps) + 1

			if(maxMaps > maxK1){maxK1 <- maxMaps}
			namesUnique <- unique(maps)
			lU <- length(namesUnique)
			l <- l/2
			if(l == lU){
				cat(line,"\n",file = conOut,append=TRUE)
			}else{
				#replicates <- as.numeric(names(which(table(maps)>1)))
				#singles <- as.numeric(names(which(table(maps)==1)))
				newLL <- numeric(lU) 
				for(k in 1:lU-1){
					newLL[k] <- log(sum(exp(ln1[2*which(maps == namesUnique[k])+1])))
				}
				newLL[lU] <- ln1[2*l+1]
				checkInf <- which(is.finite(newLL)==FALSE)
				lci <- length(checkInf)
				if (lci > 0){newLL[checkInf] <- rep(minThreshold,lci)}
				perm <- c(order(newLL[-lU],decreasing=TRUE),lU)
				cat(c("NA",lU,paste(namesUnique[perm],newLL[perm])),"\n",file=conOut,append=TRUE)
				weirdos <- weirdos + 1
			}
		}
		if (i%%1000000 == 0){print(i)}
		}
		print(i)
		close(con)
		cat(paste("....found ",pseudoreads, " pseudo-reads.",sep=""),sep="\n")
		cat(paste("....found ",weirdos, " reads mapping at same transcripts.",sep=""),sep="\n")
		cat(paste("....found ",badIndex, " reads mapping at weird places.",sep=""),sep="\n")
		iter <- iter + 1
		}
		close(conOut)
		close(clustOut1)
		cat(paste(maxK1),"\n",file = "maxK1")


	}else{
		conOut = file(paste("conditionB_1New.prob",sep=""),open="w")
		clustOut1 <- file(paste("B_1First.temp",sep=""),open = "w")
		iter <- 1
		for(outFile in filesB){
		cat(paste("correcting same transcript alignments at file ", outFile,".",sep=""),sep="\n")

		con = file(outFile)
		open(con)
		if(iter == 1){
		for(i in 1:4){
			line <- readLines(con, n = 1)
			cat(line,"\n",file = conOut,append=TRUE)
		}}else{
		for(i in 1:4){
			line <- readLines(con, n = 1)
		}
		}
		i <- 4
		weirdos <- 0
		badIndex <- 0
		pseudoreads <- 0
		while (length(line <- readLines(con, n = 1, warn = FALSE)) > 0) {
		i <- i + 1

		ln1 <- as.numeric(strsplit(line, split=" ")[[1]][-1])
		l <- 2*ln1[1]
		maps <- ln1[seq(2,l,by = 2)]
		if(sum(maps)==0){pseudoreads <- pseudoreads + 1}else{
			myRange <- diff(range(maps[-l/2]))
			if(myRange > maxDistance){
				myProb <- exp(ln1[seq(3,l,by = 2)])
				myProb[is.na(myProb)] <- 0
				myProb[is.infinite(myProb)] <- 0
				if (sum(myProb) == 0){myProb = rep(1,l/2 - 1)}
				randInd <- maps[sample(l/2 - 1,1,prob = myProb)]
				subSet <- which(abs(maps[-l/2] - randInd) < 10)
				ln1 <- c(length(subSet) + 1,ln1[c(rbind(2*subSet,2*subSet + 1),l,l+1)])
				maps <- c(maps[-l/2][subSet],0)
				l <- 2*(length(subSet) + 1)
				badIndex <- badIndex + 1
				line <- c(NA, ln1)
			}

			cat(ln1[2],"\n",file=clustOut1)
			maxMaps <- max(maps) + 1
			if(maxMaps > maxK2){maxK2 <- maxMaps}
			namesUnique <- unique(maps)
			lU <- length(namesUnique)
			l <- l/2
			if(l == lU){
				cat(line,"\n",file = conOut,append=TRUE)
			}else{
				#replicates <- as.numeric(names(which(table(maps)>1)))
				#singles <- as.numeric(names(which(table(maps)==1)))
				newLL <- numeric(lU) 
				for(k in 1:lU-1){
					newLL[k] <- log(sum(exp(ln1[2*which(maps == namesUnique[k])+1])))
				}
				newLL[lU] <- ln1[2*l+1]
				checkInf <- which(is.finite(newLL)==FALSE)
				lci <- length(checkInf)
				if (lci > 0){newLL[checkInf] <- rep(minThreshold,lci)}
				perm <- c(order(newLL[-lU],decreasing=TRUE),lU)
				cat(c("NA",lU,paste(namesUnique[perm],newLL[perm])),"\n",file=conOut,append=TRUE)
				weirdos <- weirdos + 1
			}
		}
		if (i%%1000000 == 0){print(i)}
		}
		print(i)
		close(con)
		cat(paste("....found ",pseudoreads, " pseudo-reads.",sep=""),sep="\n")
		cat(paste("....found ",weirdos, " reads mapping at same transcripts.",sep=""),sep="\n")
		cat(paste("....found ",badIndex, " reads mapping at weird places.",sep=""),sep="\n")
		iter <- iter + 1
		}
		close(conOut)
		close(clustOut1)
		cat(paste(maxK2),"\n",file = "maxK2")
	}

}


for (file in filesA){
file.remove(file)
}
for (file in filesB){
file.remove(file)
}

file.rename(from = "conditionA_1New.prob",to = "conditionA_1.prob")
file.rename(from = "conditionB_1New.prob",to = "conditionB_1.prob")



K <- max(c(as.numeric(read.table(file = "maxK1"))[1],as.numeric(read.table(file = "maxK2"))[1])) + 2 #(to 2 to vazw gia na bgazei panta toulaxiston 2 dead)
maxK <- K 
cat(paste("Number of Components: ",maxK - 3,sep=""),sep="\n")







#simil <- array(data = 0, dim = c(K-1,K-1))
dyn.load("../sparse.so")
ptm <- proc.time()
.C("sparse", Kmax=as.integer(K))
tr <- read.table("triplet_sparse_matrix.txt")
simil <- sparseMatrix(tr[,1],tr[,2],x = rep(1,length(tr[,1])))
ssss <- object.size(simil)*9.53674316*1e-7
tt <- as.numeric((proc.time() - ptm))[3]
cat(paste("Size = ",round(ssss,3)," Mb. Time += ",tt,sep = ""),"\n");


m <- K-1
Dead <- 0
Singles <- 0
index <- 1:m
sumIndex <- sum(as.numeric(index))
finalClusters <- numeric(m)
dimSimil <- dim(simil)[1]
while(sumIndex > 0){

        #i <- index[1]
        i <- which(index>0)[1]
	if(i <= dimSimil){
        index2 <- which(simil[i,]>0)}else{index2 <- integer(0)}
        if(length(index2)>0){
                condition = TRUE
                k <- 0
                indexOld <- index2
                checkIndex<-index2
                while(length(checkIndex)>0){
                        #print(index2)
                        for (j in checkIndex){
                                index2 <- unique(c(index2,which(simil[j,]>0)))
                        }
                        checkIndex <- index2[which(is.na(match(index2,indexOld))==TRUE)]
                        indexOld <- index2
                        #print(checkIndex)              
                
                }

                index2 <- sort(index2)
                #print("***********cluster************** ")
                #print(index2)
                label <- min(index2)
                finalClusters[index2] <- rep(label,length(index2))
                index[index2] <- rep(0,length(index2)) 
        
        }else{
                index[i] = 0
                Dead <- c(Dead,i)
                #print(paste("dead transcript: ", i))
        }
        sumIndex <- sum(as.numeric(index))
        
}


Dead <- Dead[-1]

length(which(table(finalClusters)==1))

length(table(finalClusters))

write.table(finalClusters,file=paste("finalClusters",".txt",sep=""),col.names = FALSE,row.names = FALSE)

#simil <- 0

########################################################################################################
########################################################################################################
########################################################################################################
########################################################################################################
########################################################################################################


#compute clust Alignments New


ptm <- proc.time()
nClusters <- length(table(finalClusters))
freqClusters <- numeric(nClusters)
names(freqClusters) <- names(table(finalClusters))
empty <- names(freqClusters)
readsPerReplicate <- array(data = NA, dim = c(nClusters,length(fls)))
j<-0
for(outFile in fls){
	freqClusters <- numeric(nClusters)
	firstClust <- as.numeric(read.table(paste(outFile,"First.temp",sep=""))[,1])
	j <- j + 1
	previous <- freqClusters
	names(freqClusters) <- names(table(finalClusters))
	cat(paste("computing cluster alignments at file ", outFile,".",sep=""),sep="\n")
	myTable <- table((finalClusters[firstClust]))
	freqClusters[names(myTable)] <- myTable
	readsPerReplicate[,j] <- freqClusters # - previous
}

rownames(readsPerReplicate) <- names(freqClusters)
colnames(readsPerReplicate) <- fls
firstClust <- 0

newReadsPerReplicate <- readsPerReplicate
Clusters <- finalClusters

proc.time() - ptm

#edw emeina

#######################################################################################################
#######################################################################################################
#######################################################################################################
#######################################################################################################
#######################################################################################################


lowThreshold <- 2
smallClusters <- names(which(apply(newReadsPerReplicate,1,min)<lowThreshold))
smallClusters <- smallClusters[-1] # exclude the empty transcripts
sum(readsPerReplicate[smallClusters,]) # this is the number of observations assigned to smallClusters

if(sum(readsPerReplicate[smallClusters,]) > 0){
Grouped <- as.numeric(names(which(apply(newReadsPerReplicate,1,min)>=lowThreshold)))
nGrouped <- length(Grouped)
randGroupedLabel <- Grouped[floor(nGrouped*runif(length(smallClusters)) + 1)]
Clusters <- finalClusters
iter <- 0
newLength <- length(table(finalClusters)) - length(smallClusters)
newReadsPerReplicate <- array(data = 0, dim = c(newLength,length(fls)))
rownames(newReadsPerReplicate) <- c(0,Grouped)
colnames(newReadsPerReplicate) <- colnames(readsPerReplicate)
newReadsPerReplicate[as.character(Grouped),] <- readsPerReplicate[as.character(Grouped),]
for(k in as.numeric(smallClusters)){
	iter <- iter + 1
	smallClusterMembers <- which(finalClusters == k)
	Clusters[smallClusterMembers] <- rep(as.numeric(randGroupedLabel[iter]),length(smallClusterMembers))
	newReadsPerReplicate[as.character(randGroupedLabel[iter]),] <- newReadsPerReplicate[as.character(randGroupedLabel[iter]),] + readsPerReplicate[as.character(k),]
}
}

write.table(Clusters,file=paste("Clusters",".txt",sep=""),col.names = FALSE,row.names = FALSE)
length(table(Clusters))

#this is to relabel the transcripts inside each cluster.
# The relabelled ones will be used to the new files.
if(length(Dead)>0){
	namesClusters <- as.numeric(names(table(Clusters)))
	relabelledClusters <- Clusters
	for (k in namesClusters[-1]){
		index <- which(Clusters == k)
		relabelledClusters[index] <- 1:length(index)	
	}
}else{
	namesClusters <- as.numeric(names(table(Clusters)))
	relabelledClusters <- Clusters
	for (k in namesClusters){
		index <- which(Clusters == k)
		relabelledClusters[index] <- 1:length(index)	
	}
}



#orderedNames corresponds to the names of the files 0,1,2,3,...,length(table(Clusters))

if(length(Dead)>0){
	orderedNames <- order(namesClusters) - 1
}else{
	orderedNames <- order(namesClusters) 
}

names(orderedNames) <- namesClusters




################################################################################


rs <- rowSums(newReadsPerReplicate)
write.table(newReadsPerReplicate,file=paste("nReadsPerCluster",".txt",sep=""),col.names = FALSE,row.names = FALSE)


################################################################################


	nMaps <- colSums(newReadsPerReplicate)
	it <- 0
	# create directories and files
	nClusters <- length(table(Clusters))
	foreach(fileIter=1:2) %dopar% {
		ext <- fls[fileIter]
		it <- it + 1
		cat(paste("Creating files for sample ",ext,sep=""),sep="\n")
		#open prob file
		con = file(paste("condition",ext,".prob",sep=""))
		open(con)
		newProb <- file(paste("condition",ext,".panos",sep=""),open = "w")
		cat(readLines(con, n = 1))
		cat("\n")
		thisLine <- readLines(con, n = 1)
		#nMaps[it] <- as.numeric(strsplit(thisLine,split = " ")[[1]][3])
		cat(readLines(con, n = 1))
		cat("\n")
		cat(readLines(con, n = 1))
		cat("\n")
		while (length(line <- readLines(con, n = 1, warn = FALSE)) > 0) {
			ln1 <-  as.numeric(strsplit(line, split=" ")[[1]][-1])
			#change transcript names according to the relabelling
			old_lab <- ln1[2]
			maps <- seq(2,2*(ln1[1] - 1),by = 2)
			ln1[maps] <- relabelledClusters[ln1[maps]]
			fileID <- orderedNames[[as.character(Clusters[old_lab])]]
			ln1<-c(fileID,ln1)
			cat(ln1,"\n",file=newProb)
		}
		close(con)
		close(newProb)
	}

	# create sub-directories for each cluster with links and reads per Cluster information
	dir.create("clusters")
	setwd("clusters/")
	options(digits = 16)

	for(k in 1:(nClusters-1)){
		myDir <- paste("cluster_",k,sep="")
		dir.create(myDir)
		setwd(myDir)

		file.create("data.tr")
#		cat(paste("# M ",as.numeric(table(Clusters)[k+1])," ",k," ",rs[k+1],sep = ""),file="data.tr",sep="\n")
		pix <- as.numeric(table(Clusters)[k+1])
		cat(paste("# ",maxK-pix," ",pix," ",k," ",rs[k+1],sep = ""),file="data.tr",sep="\n")
		myFile <- "conditionA.weights"
		file.create(myFile)
		con <- file(myFile)
		open(con)
		for (rep in 1:nreps1){
			cat(as.numeric(newReadsPerReplicate[k+1,rep]/nMaps[rep]),file=myFile,append = TRUE)
			cat("\n",file=myFile,append = TRUE)
		}
		close(con)
		myFile <- "conditionB.weights"
		file.create(myFile)
		con <- file(myFile)
		open(con)
		for (rep in 1:nreps2){
			cat(as.numeric(newReadsPerReplicate[k+1,nreps1 + rep]/nMaps[nreps1 + rep]),file=myFile,append = TRUE)
			cat("\n",file=myFile,append = TRUE)
		}
		close(con)
		setwd("../")
	}


##############################################################################
##############################################################################
##############################################################################
##############################################################################
##############################################################################


setwd("../../")
rs <- rs[-1]
if(1>2){
	nJobs <- 60
	if(K < 1000){
		nJobs <- 10
	}
	# sharing jobs based on number of reads
	threshold <- sum(rs)/nJobs
	i <- 0
	groupIndex <- 0
	vec2 <- array(0,c(nJobs,2))
	j2 <- 0
	while((i < nClusters-2)){
		j1 <- i + 1
		nObs <- 0; 
		while((nObs < threshold)&(i < nClusters - 1)){
			i<-i+1;
			nObs <- nObs + rs[i]
			#cat(paste("i = ",i,", nObs = ",nObs,sep = ""),"\n")
		}
		j2 <- i
		groupIndex <- groupIndex + 1
		vec2[groupIndex,] <- c(j1,j2)
		cat(paste("group ",groupIndex,": {",j1,",...,",j2,"}. Nobs = ",nObs,sep = ""),"\n")
		cat("\n")
	}
	if (j2 < nClusters - 1){
		j1 <- j2 + 1
		j2 <- nClusters - 1
		nObs <- sum(rs[j1:j2])
		groupIndex <- groupIndex + 1
		vec2[groupIndex,] <- c(j1,j2)
		cat(paste("group ",groupIndex,": {",j1,",...,",j2,"}. Nobs = ",nObs,sep = ""),"\n")
		cat("\n")
	}
	nJobs <- groupIndex
}else{
	nJobs <- nClusters - 1
	vec2 <- array(0,c(nJobs,2))
	for(k in 1:nJobs){
		vec2[k,] <- k
	}
}




dir.create("jobs")
setwd("jobs/")
priorPerm <- order(rs,decreasing = TRUE)
for (k in 1:nJobs){
        vec <- vec2[k,]
        job <- file(paste("job_",k,".sh",sep = ""),open = "w")
        cat("#!/bin/bash","\n",file = job)
        cat(paste("cd ../tmp/clusters/cluster_",priorPerm[k],sep=""),"\n",file = job)
        cat("head -1 data.tr","\n",file = job)
	cat("collapsed > sampler.log","\n",file = job)  # collapsed sampler
        close(job)
        system(paste("chmod u+x job_",k,".sh",sep=""))
        
}



#write the JobArray for csf
job <- file(paste("jobArray",".sh",sep = ""),open = "w")
cat("#!/bin/bash","\n",file = job)
cat("#$ -S /bin/bash","\n",file = job)
cat("#$ -V","\n",file = job)
cat("#$ -cwd","\n",file = job)
#cat("#$ -pe smp.pe 2","\n",file = job)
#cat("export OMP_NUM_THREADS=$NSLOTS","\n",file = job)
cat(paste("#$ -t 1-",nJobs,sep=""),"\n",file = job)
cat("./job_${SGE_TASK_ID}.sh","\n",file = job)
close(job)
setwd("../")

deleteThreshold <- floor(nJobs/10)

if(1 < 10){
#DON'T delete files
# write the script for GNU parallel
job <- file(paste("parallelGNU",".bash",sep = ""),open = "w")
cat("#!/bin/bash","\n",file = job)
#cat("echo \"Running Reversible Jump MCMC:\"","\n",file = job)
cat("echo \"(#M nTranscripts ClusterID nReads)\"","\n",file = job)
cat("cd jobs","\n",file = job)
#cat(paste("parallel --no-notice ./job_{}.sh ::: {1..",nJobs,"}",sep=""),"\n",file = job)
cat(paste("parallel --gnu ./job_{}.sh ::: {1..",nJobs,"}",sep=""),"\n",file = job)
cat("wait","\n",file = job)
cat("cd ..","\n",file = job)
cat("echo \"Processing output...\"","\n",file = job)
cat("R CMD BATCH  getClusters.R","\n",file = job)
close(job)
system(paste("chmod u+x parallelGNU",".bash",sep=""))
}else{
# The awk procedure it is too slow => abandon for now.
## delete processed clusters after deleteThreshold finished processing
job <- file(paste("parallelGNU",".bash",sep = ""),open = "w")
cat("#!/bin/bash","\n",file = job)
#cat("echo \"Running Reversible Jump MCMC:\"","\n",file = job)
cat("echo \"(#M nTranscripts ClusterID nReads)\"","\n",file = job)
cat("cd jobs","\n",file = job)
cat(paste("parallel --gnu ./job_{}.sh ::: {1..",deleteThreshold,"}",sep=""),"\n",file = job)
cat("wait","\n",file = job)
cat("echo \"     Pause for a while to delete some unnecessary things...\"","\n",file = job)
cat("cd ../","\n",file = job)
myExp <- "ls -lah tmp/conditionA_1.panos tmp/conditionB_1.panos"
cat(myExp,"\n",file = job)
for(k in 1:deleteThreshold){
myExp <- paste("awk '$1!=",priorPerm[k],"' tmp/conditionA_1.panos > tmp/temp.panos",sep="")
cat(myExp,"\n",file = job)
cat("mv tmp/temp.panos tmp/conditionA_1.panos","\n",file = job)
myExp <- paste("awk '$1!=",priorPerm[k],"' tmp/conditionB_1.panos > tmp/temp.panos",sep="")
cat(myExp,"\n",file = job)
cat("mv tmp/temp.panos tmp/conditionB_1.panos","\n",file = job)
myExp <- paste("echo \"     Cluster ",priorPerm[k]," deleted.\"",sep="")
cat(myExp,"\n",file = job)
myExp <- "ls -lah tmp/conditionA_1.panos tmp/conditionB_1.panos"
cat(myExp,"\n",file = job)
}

cat("cd jobs","\n",file = job)

cat("echo \"     OK. Now I will continue running MCMC...\"","\n",file = job)
cat(paste("parallel --gnu ./job_{}.sh ::: {",deleteThreshold + 1,"..",nJobs,"}",sep=""),"\n",file = job)
cat("wait","\n",file = job)

cat("cd ..","\n",file = job)
cat("echo \"Processing output...\"","\n",file = job)
cat("R CMD BATCH  getClusters.R","\n",file = job)
close(job)
system(paste("chmod u+x parallelGNU",".bash",sep=""))


}





