setwd("clusters/")
numClusters <- length(list.files(pattern = "cluster_"))
ramones <- data.frame(
geneName=character(),
trName=character(), 
trID=numeric(),
nTrPerGene=numeric(),
trIDperGene=numeric(),
ProbDE=numeric(),
theta=numeric(),
w=numeric(), 
GeneProbDE=numeric(),
stringsAsFactors=FALSE) 
for(k in 1:numClusters){
	setwd(paste("cluster_",k,"/",sep=""))
	path <- "state_vector.txt";
	checkFile <- file.exists(path)
	if(checkFile == TRUE){
		myNames <- read.table("names.info")
		models <- as.numeric(read.table(path)[,1]);
		th <- as.numeric(read.table("theta.txt")[,1]);
		wou <- as.numeric(read.table("w.txt")[,1]);
		gene <- tail(read.table("acceptance.txt")[,1],1)
		ramones <- rbind(ramones,cbind(myNames,models,th,wou,rep(gene,length(models))))
	}
	setwd("../")
	#system(paste("rm -r cluster_",k,sep=""))
	if(k%%1000 == 0){print(k)}
}
ramones[,1] <- as.character(ramones[,1])
ramones[,2] <- as.character(ramones[,2])
colnames(ramones) <- c("geneName","trName","trID","nTrPerGene","trIDperGene","ProbDE","theta","w","GeneProbDE")
naIndex <- which(is.na(ramones$geneName)==TRUE)
lNA <- length(naIndex)
if( lNA > 0){ramones[naIndex,1] <- paste('missingName',1:lNA)}
setwd("../")
#########################################################################################
#	Computing filtered Probability of DTU:
#########################################################################################
geneNames <- unique(ramones$geneName)
filteredGeneProbs <- numeric(dim(ramones)[1])
newRamones <- cbind(ramones,numeric(dim(ramones)[1]),numeric(dim(ramones)[1]))
colnames(newRamones) <- c("geneName","trName","trID","nTrPerGene","trIDperGene","ProbDE","theta","w","GeneProbDE","GeneProbDE.FDRadjusted","switching")
for (gene in geneNames){
	trans <- which(ramones$geneName == gene)
	dominantTranscripts <- c(order(ramones$theta[trans],decreasing=T)[1],order(ramones$w[trans],decreasing=T)[1])
        deProb <- ramones$GeneProbDE[trans][1]
	switchingEvent <- 0
        if (diff(ramones$theta[trans][dominantTranscripts])*diff(ramones$w[trans][dominantTranscripts]) < 0){
                deProb2 <- deProb
		switchingEvent <- 1
        }else{
		if((abs(ramones[trans,]$theta[dominantTranscripts[1]]-ramones[trans,]$w[dominantTranscripts[1]]) + abs(ramones[trans,]$theta[dominantTranscripts[2]]-ramones[trans,]$w[dominantTranscripts[2]]))/2 > 0.2){deProb2 <- deProb; switchingEvent <- 1}else{deProb2 <- 0}
        }
	newRamones$GeneProbDE.FDRadjusted[trans] <- rep(deProb2,length(trans))
	newRamones$switching[trans] <- rep(switchingEvent,length(trans))	
}
ramones <- newRamones
write.table(ramones,row.names=FALSE,file = "../transcriptLevelEstimates.txt",quote = FALSE)
# filtering:
u <- abs(ramones$theta - ramones$w)
ind <- which((ramones$GeneProbDE > 0.6)&(u < 0.075))
if (length(ind) > 0){ramones$GeneProbDE[ind] <- rep(0,length(ind))}
##############################################################################################
#	Computing Probability of DTU weighted according to absolute differences (Gene ProbDE):
##############################################################################################
perm <- order(u,decreasing=TRUE)
cjPerm <- ramones[perm,]
u <- cjPerm[,1]
v <- unique(cjPerm[,1])
perm <- match(v,u)
diFF <- abs(cjPerm[perm,]$theta-cjPerm[perm,]$w)
cjPerm <- cbind(cjPerm[perm,c(1,4,9,10)],diFF,cjPerm[perm,11])

nReads <- read.table("nReads.txt")
cjPerm <- cbind(cjPerm,nReads[as.character(cjPerm$geneName),])
K <- dim(cjPerm)[1]
sigTranscripts <- numeric(K)
p <- cjPerm
orderedP <- p$GeneProbDE
K <- dim(p)[1]
myList <-  1 - orderedP[1]
k <- 1 
sigTranscripts[k] <- myList
for (k in 2:K){
	myList <- myList + 1.0*(1 - orderedP[k])
	sigTranscripts[k] <- myList/k
}
v <- 1 - sigTranscripts
cjPerm <- cbind(cjPerm,v)
filterOut <- which(cjPerm[,6] == 0)
cjPerm[,6] <- v
cjPerm[filterOut,6] <- rep(0,length(filterOut))
colnames(cjPerm)[c(4,5,6,7,8,9)] <- c('Pswitching','maxDiff','FDR','nReadsA','nReadsB','FDRraw')



#########################################################################################
write.table(cjPerm,row.names=FALSE,file = "../withinGeneEstimates.txt",quote = FALSE)
library(fields)
perm <- order(ramones$GeneProbDE.FDRadjusted, decreasing = TRUE)
pdf(file = "../de.pdf",width=6,height=6)
plot(c(0,1),c(0,1),type="n",xlab = "theta",ylab="w");
rect(par("usr")[1],par("usr")[3],par("usr")[2],par("usr")[4],col = "gray")
points(ramones[perm,c(7,8)],col = color.scale(ramones[perm,6]),pch = 16,cex=0.5);
dev.off()


