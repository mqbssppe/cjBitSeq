workdir <- getwd()
myPath <- paste(workdir,"/cluster_",sep="")
setwd("tmp/")

clusters<-read.table("Clusters.txt")
clusters<-clusters[,1]
#clusters <- clusters[1:(length(clusters)-2)]
Clusters <- clusters
namesClusters <- as.numeric(names(table(Clusters)))
relabelledClusters <- Clusters
nReads <- read.table("nReadsPerCluster.txt")
nReads <- colSums(nReads)
for (k in namesClusters[-1]){
	index <- which(Clusters == k)
	relabelledClusters[index] <- 1:length(index)    
}
orderedNames <- order(namesClusters) - 1
names(orderedNames) <- namesClusters
#### this is to compute thetaNew, wNew and RPKM

K<-length(relabelledClusters) - 2
deFinal <- numeric(K)
thetaFinal <- rep(1/nReads[1],K)
wFinal <- rep(1/nReads[2],K)
Dead <- which(clusters[1:(length(clusters)-2)] == 0)
deFinal[Dead] <- rep(0,length(Dead))
#thetaFinal[Dead] <- rep(1/nReads[1],length(Dead))
#wFinal[Dead] <- rep(1/nReads[2],length(Dead))
setwd("clusters/")
numClusters <- length(list.files(pattern = "cluster_"))
acceptanceRates <- numeric(numClusters)
for(k in 1:numClusters){
	path <- paste("cluster_",k,"/state_vector.txt",sep="");
	checkFile <- file.exists(path)
	if(checkFile == TRUE){
		models <- as.numeric(read.table(path)[,1]);
		models <- models[-c(1,length(models))];
		path <- paste("cluster_",k,"/theta.txt",sep="")
		th <- as.numeric(read.table(path)[,1]);
		th <- th[-c(1,length(th))]
		path <- paste("cluster_",k,"/w.txt",sep="")
		wou <- as.numeric(read.table(path)[,1]);
		wou <- wou[-c(1,length(wou))]
		clusterName <- as.numeric(names(orderedNames[k+1]))
		transcriptIndex <- which(Clusters == clusterName)
		deFinal[transcriptIndex] <- models
		thetaFinal[transcriptIndex] <- as.numeric(th)
		wFinal[transcriptIndex] <- as.numeric(wou)
		path <- paste("cluster_",k,"/acceptance.txt",sep="")
		acceptanceRates[k] <- mean(as.numeric(read.table(path)[,1]))
	}
	system(paste("rm -r cluster_",k,sep=""))
	if(k%%1000 == 0){print(k)}
}
setwd("../")
system("rm -r clusters")
system("rm *.panos *.prob *.temp Clusters.txt maxK1 maxK2")
setwd("../")

thetaFinal <- thetaFinal/sum(thetaFinal)
wFinal <- wFinal/sum(wFinal)
ramones <- cbind(log(thetaFinal),log(wFinal),deFinal)

sigTranscripts1 <- sigTranscripts2 <- sigTranscripts3 <- rep("non-DE",K)

# FDR: alpha = 0.01
alpha <- 0.01
p <- ramones
perm <- order(p[,3],decreasing = TRUE)
orderedP <- p[perm,3]
K <- dim(p)[1]
myList <-  1 - orderedP[1]
k <- 1 
criterion <- myList
while ((criterion < alpha)&(k < K)){
	k <- k + 1 
	myList <- myList + 1 - orderedP[k]
	criterion <- myList/k
}
if(k > 1){
	sigTranscripts1[perm[1:(k-1)]] <- rep("DE",k-1)
}
# FDR: alpha = 0.05
alpha <- 0.05
p <- ramones
perm <- order(p[,3],decreasing = TRUE)
orderedP <- p[perm,3]
K <- dim(p)[1]
myList <-  1 - orderedP[1]
k <- 1 
criterion <- myList
while ((criterion < alpha)&(k < K)){
	k <- k + 1 
	myList <- myList + 1 - orderedP[k]
	criterion <- myList/k
}
if(k > 1){
sigTranscripts2[perm[1:(k-1)]] <- rep("DE",k-1)
}

# FDR: alpha = 0.1
alpha <- 0.1
p <- ramones
perm <- order(p[,3],decreasing = TRUE)
orderedP <- p[perm,3]
K <- dim(p)[1]
myList <-  1 - orderedP[1]
k <- 1 
criterion <- myList
while ((criterion < alpha)&(k < K)){
	k <- k + 1 
	myList <- myList + 1 - orderedP[k]
	criterion <- myList/k
}
if(k > 1){
sigTranscripts3[perm[1:(k-1)]] <- rep("DE",k-1)
}
ramones <- cbind(ramones,sigTranscripts1, sigTranscripts2, sigTranscripts3 )

write.table(ramones,row.names=FALSE,col.names = c("log-Theta","log-W","Prob(D.E)","fdr_0.01","fdr_0.05","fdr_0.1"),file = "estimates.txt",quote = FALSE)
write.table(acceptanceRates,file = "acceptanceRatePerCluster.txt")

system("rm *.prob parallelGNU.bash partitionMerge2.R sparse.so")
system("rm -r jobs")




