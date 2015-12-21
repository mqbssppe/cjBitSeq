## this file should be executed from the command line as follows:
##
##	 R CMD BATCH '--args fastaFile gtfFile outputFile' prepareAnnotationForcjBitSeq.R
##
## where:
##
##	fastaFile   : the full path to the fasta file
##	gtfFile     : the full path to the gtf file.
##	outputFile  : output file containing the gene - transcript information for cjBitSeq.
##------------------------------------------------------------------------------------------------------------------
## the fasta file should contain the transcript names at the header section before the transcript sequence
## the gtf file should contain an entry named 'gene_id' with the gene names and an entry named 'transcript_id' with the transcript names.
##
## FASTA File header format example:
##
## >20 FBtr0300689 2L+ 7529-8116,8193-9484
## 
## GTF file format example
##2L	protein_coding	exon	7529	8116	.	+	.	 gene_id "FBgn0031208"; transcript_id "FBtr0300689"; exon_number "1"; gene_name "CG11023"; gene_biotype "protein_coding"; transcript_name "CG11023-RB"; exon_id "FBgn0031208:1";
##2L	protein_coding	CDS	7680	8116	.	+	0	 gene_id "FBgn0031208"; transcript_id "FBtr0300689"; exon_number "1"; gene_name "CG11023"; gene_biotype "protein_coding"; transcript_name "CG11023-RB"; protein_id "FBpp0289913";
##2L	protein_coding	start_codon	7680	7682	.	+	0	 gene_id "FBgn0031208"; transcript_id "FBtr0300689"; exon_number "1"; gene_name "CG11023"; gene_biotype "protein_coding"; transcript_name "CG11023-RB";
##2L	protein_coding	exon	8193	9484	.	+	.	 gene_id "FBgn0031208"; transcript_id "FBtr0300689"; exon_number "2"; gene_name "CG11023"; gene_biotype "protein_coding"; transcript_name "CG11023-RB"; exon_id "FBgn0031208:3";
##2L	protein_coding	CDS	8193	8607	.	+	1	 gene_id "FBgn0031208"; transcript_id "FBtr0300689"; exon_number "2"; gene_name "CG11023"; gene_biotype "protein_coding"; transcript_name "CG11023-RB"; protein_id "FBpp0289913";
##2L	protein_coding	stop_codon	8608	8610	.	+	0	 gene_id "FBgn0031208"; transcript_id "FBtr0300689"; exon_number "2"; gene_name "CG11023"; gene_biotype "protein_coding"; transcript_name "CG11023-RB";


myFiles <- commandArgs(trailingOnly = TRUE)

if(length(myFiles) != 3){stop("a fasta, a gtf and the name of the output file should be provided.")}
if( file.exists(myFiles[1]) == FALSE){stop("fasta file not found.")}
if( file.exists(myFiles[2]) == FALSE){stop("gtf file not found.")}
if( file.exists(myFiles[3]) == TRUE){stop("output file exists, please provide different name.")}


fasta <- myFiles[1] #'/mnt/mr01-data01/mqbssppe/spliceSimulationsZurich/drosophila/annotation/TopHatTranscriptomeIndex/Drosophila_melanogaster.BDGP5.70.dna.toplevel.fa'

gtf <- myFiles[2] #'/mnt/mr01-data01/mqbssppe/spliceSimulationsZurich/drosophila/annotation/Drosophila_melanogaster.BDGP5.70.dna.toplevel.gtf'

outputFile <- myFiles[3]

system(paste('grep ">" ',fasta,' > fastanames.txt',sep=''))#'grep ">" /mnt/mr01-data01/mqbssppe/cufflinksHumanAnnotation/transcriptome_data/known.fa > fastanames.txt')
myFile <- file("fastanames.txt",open="r")
K <- as.numeric(strsplit(system("wc -l fastanames.txt",intern = TRUE),split = " ")[[1]][1])
trNames <- numeric(K)
i <- 0
while (length(oneLine <- readLines(myFile, n = 1, warn = FALSE)) > 0){
        i <- i + 1
        trNames[i] <- strsplit(oneLine,split = " ")[[1]][2]
}
close(myFile)
write.table(trNames,file = "transcriptNames.tr",quote = FALSE,row.names = FALSE, col.names = FALSE)
#################################################################################################################

myFile <- file(gtf,open = "r")
d1 <- 100000
temp <- array(data = NA,dim = c(d1,2))
i <- 0
while(length(line <- readLines(myFile, n = 1, warn = FALSE)) > 0){
        i <- i + 1
        ln1 <- strsplit(line,split=c("gene_id"))[[1]][2]
        ln2 <- strsplit(ln1,split = "transcript_id")
        geneName <- strsplit(ln2[[1]][1],split= "\"")[[1]][2]
        trName <- strsplit(ln2[[1]][2],split = "\"")[[1]][2]    
        if(i > d1){
                cat(paste("line:",i),"\n")
                d1 <- d1 + 100000
                tempNew <- array(data = NA,dim = c(d1,2))
                tempNew[1:(i-1),] <- temp
                temp <- tempNew
                tempNew <- 0
        }
        temp[i,] <- c(geneName,trName)
}
close(myFile)
temp <- temp[1:i,]
myNames <- unique(temp)
if(length(trNames)!= dim(myNames)[1]){message("number of transcripts in fasta file is not equal to number of transcripts in gtf file")}else{message("number of transcripts in fasta file is equal to number of transcripts in gtf file")}
notFoundTranscripts <- which(is.na(match(trNames,myNames[,2]))==TRUE)
if (length(notFoundTranscripts) > 0){message("not compatible transcripts detected")}
colnames(myNames) <- c("geneID","transcriptID")
write.table(myNames,file = "geneANDtranscripts.txt",quote =FALSE)


####################################################################################################################
l <- K
finalNames <- trNames

knownIsoformNames <- myNames
v <- knownIsoformNames
knownIsoformNames <- cbind(as.character(v[,1]),as.character(v[,2]))
noNameIndex <- which(knownIsoformNames[,2] == "")
knownIsoformNames[noNameIndex,2] <- knownIsoformNames[noNameIndex[1],1]

perm <- match(finalNames,as.character(knownIsoformNames[,2]))

knownIsoformNames <- knownIsoformNames[perm,]
knownIsoformNames <- cbind(knownIsoformNames,1:l)

nTr <- numeric(l)
for(i in 1:l){
        j = knownIsoformNames[i,1]
        nTr[i] <- length(which(knownIsoformNames[,1]==j))
}
knownIsoformNames <- cbind(knownIsoformNames,nTr)
trIDperGene <- numeric(l)
tt <- table(knownIsoformNames[,1])
gMax <- length(tt)
for(i in names(tt)){
        myIndex <- which(knownIsoformNames[,1] == i)
        trIDperGene[myIndex] <- 1:length(myIndex)
}

knownIsoformNames <- cbind(knownIsoformNames,trIDperGene)
write.table(knownIsoformNames,outputFile,row.names=FALSE,col.names=c("geneID","trName","trID","nTrPerGene","trIDperGene"),quote=FALSE)



