newReadsPerReplicate <- read.table("tmp/nReadsPerCluster.txt")
rs <- rowSums(newReadsPerReplicate)
nJobs <- 1000
if(K < 1000){
        nJobs <- 10
}



rs <- rs[-1]


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
for (k in 1:nJobs){
        vec <- vec2[k,]
        job <- file(paste("job_",k,".sh",sep = ""),open = "w")
        cat("#!/bin/bash","\n",file = job)
        cat("cd tmp","\n",file = job)
        cat("echo \"Running Reversible Jump \"","\n",file = job)
        cat("cd clusters","\n",file = job)
        cat(paste("for i in {",vec[1],"..",vec[2],"}",sep = ""),"\n",file = job)
        cat("do","\n",file = job)
        cat("echo \"        Running RJMCMC at: cluster_$i\"","\n",file = job)
        cat("cd cluster_$i","\n",file = job)
        cat("echo \"              number of transcripts:\"","\n",file = job)
        cat("head -1 data.tr","\n",file = job)
        cat("rjFULL > rjmcmc.log","\n",file = job)
        cat("cd ..","\n",file = job)
        cat("done","\n",file = job)
        cat("cd ..","\n",file = job)
        close(job)
        system(paste("chmod u+x job_",k,".sh",sep=""))
        
}
job <- file(paste("jobArray",".sh",sep = ""),open = "w")
cat("#!/bin/bash","\n",file = job)
cat("#$ -S /bin/bash","\n",file = job)
cat("#$ -V","\n",file = job)
cat("#$ -cwd","\n",file = job)
cat("#$ -pe smp.pe 4","\n",file = job)
cat("export OMP_NUM_THREADS=$NSLOTS","\n",file = job)
cat(paste("#$ -t 1-",nJobs,sep=""),"\n",file = job)
cat("./job_${SGE_TASK_ID}.sh","\n",file = job)
close(job)

