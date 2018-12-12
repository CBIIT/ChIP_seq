###Script for transforming ROSE outputs into BED files

setwd("K:/projects/ChIP_seq")
project.folder = "projects/MCC/ROSE/"
files <- list.files(path=project.folder, pattern="*12500.combined.txt", full.names=T, recursive=FALSE)
#files <- files[grep("MCC13", files)]  #subset if needed

lapply(files, function(x) {
  
  #read in txt file
  ROSE <- read.table(x,sep="\t", header=T)
  ROSE.bed <- ROSE[,c(2,3,4,6,7,8,9,13)]
  ROSE.bed = ROSE.bed[!grepl("gl",ROSE.bed$CHROM),]
  colnames(ROSE.bed)[1] = "#CHR"
  ROSE.top5K.bed = ROSE.bed[1:5000,]
  ###output cluster information
  filename <- x
  filename <- gsub("txt","bed",filename)
  filename.top5K <- gsub("bed","top5K.bed",basename(filename))
  write.table(ROSE.bed, file=filename, sep="\t", row.names=FALSE, col.names=TRUE, quote=FALSE)
  write.table(ROSE.top5K.bed, file=paste(dirname(filename),filename.top5K,sep="/"), sep="\t", row.names=FALSE, col.names=TRUE, quote=FALSE)
})

files.top5K.beds <- list.files(path=project.folder, pattern="*top5K.bed", full.names=F, recursive=FALSE)
write.table(files.top5K.beds, file=paste(project.folder,"files.top5K.beds",sep=""), sep="\t", row.names=FALSE, col.names=F, quote=FALSE)

###Get SEs only
lapply(files, function(x) {
  
  #read in txt file
  ROSE <- read.table(x,sep="\t", header=T)
  ROSE.SE = subset(ROSE,ROSE$isSuper == 1)
  ROSE.SE <- ROSE.SE[,c(2,3,4,6,7,8,9,13)]
  ROSE.SE = ROSE.SE[!grepl("gl",ROSE.SE$CHROM),]
  colnames(ROSE.SE)[1] = "#CHR"
  
  ###output cluster information
  filename <- x
  filename.SE <- gsub("txt","SE.bed",basename(filename))
  write.table(ROSE.SE, file=paste(dirname(filename),filename.SE,sep="/"), sep="\t", row.names=FALSE, col.names=TRUE, quote=FALSE)
  
})

###Get all Enhancer beds, with enhancer rank as 4th column (for finding enhancer/gene callouts in ROSE data)
lapply(files, function(x) {
  
  #read in txt file
  ROSE <- read.table(x,sep="\t", header=T)
  ROSE.SE = subset(ROSE,ROSE$isSuper == 1)
  ROSE.SE <- ROSE.SE[,c(2,3,4,9)]
  ROSE.SE = ROSE.SE[!grepl("gl",ROSE.SE$CHROM),]
  colnames(ROSE.SE)[1] = "#CHR"
  
  ###output cluster information
  filename <- x
  filename.SE <- gsub("txt","rank.bed",basename(filename))
  write.table(ROSE.SE, file=paste(dirname(filename),filename.SE,sep="/"), sep="\t", row.names=FALSE, col.names=TRUE, quote=FALSE)
  
})