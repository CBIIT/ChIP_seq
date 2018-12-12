###Gather coltron outputs
#find the file:  assign p-value for chosen ROSE outputs; make this a new column in the Sample sheet.  Make sure coltron has been run for each chosen ROSE call
#set project global parameters
project.folder = "projects/MCC/coltron/"
Filename = "MCC_SE_TFs"

#load sample sheet from excel file
  require(XLConnect)
  setwd("K:/projects/ChIP_seq")
  wb = loadWorkbook("manage_samples/ChIP_seq_samples.xlsx")
  ChIPs = readWorksheet(wb, sheet = "ChIP_seq_samples", header = TRUE)
  #pass through txt file
  write.table(ChIPs, file="manage_samples/ChIP_seq_all_samples.txt",sep="\t", col.names = T)
  ChIPs = read.table(paste(project.folder,"ChIPseq_Samples.txt",sep=""),sep="\t", header=T)

#define coltron output paths
ChIPs = ChIPs[grep("p_", ChIPs$ChosenP_path), ]
ChIPs = ChIPs[grep("MCC", ChIPs$Project), ]  

coltron = ChIPs[,c("SampleFiles","ChosenP_path")]
coltron$SE_TF_path = paste("DATA/",coltron$SampleFiles,"/MACS_Out_",coltron$ChosenP_path,"/ROSE_out_12500/coltron/",coltron$SampleFiles,"_CANIDATE_TF_AND_SUPER_TABLE.txt",sep="")
coltron$path = paste("DATA/",coltron$SampleFiles,"/MACS_Out_",coltron$ChosenP_path,"/ROSE_out_12500/coltron/",sep="")

#pass through saved file
write.table(ChIPs, file = paste(project.folder,"ChIPseq_Samples.txt",sep=""), sep="\t", row.names=FALSE, quote=FALSE)
write.table(coltron, file = paste(project.folder,"coltronSamples.txt",sep=""), sep="\t", row.names=FALSE, quote=FALSE)
coltron <- read.table(paste(project.folder,"coltronSamples.txt",sep=""), sep="\t", header=T)

########################################################################################
######### EXTRACTING FREQ and DEGREE of TFs across a list of Coltron samples############
########################################################################################


#loop through sample list, build a list of all TFs (rbind all, then isolate the unique ones)
SampleList = as.character(coltron$SampleFiles) #create Sample list for lapply

lapply(SampleList, function(x) {
  filepath <- subset(coltron, coltron$SampleFiles %in% x) #isolate single sample
  SE_TFs <- read.table(as.character(filepath$SE_TF_path),sep="\t", header=T)
  SE_TFs = as.data.frame(SE_TFs[,2])
  write.table(SE_TFs, file = paste(project.folder,Filename,".txt",sep=""), sep="\t", row.names=FALSE, quote=FALSE, col.names=FALSE, append=TRUE)
})

All_SE_TFs = read.table(paste(project.folder,Filename,".txt",sep=""),sep="\t", header=F)
All_SE_TFs.counted<-table(All_SE_TFs$V1)
All_SE_TFs.counted<-as.data.frame(All_SE_TFs.counted)
colnames(All_SE_TFs.counted) <- c("Var1","Freq.All")
All_SE_TFs.degree = All_SE_TFs.counted

#loop through a second time, but now make a column for each sample, leave 0 if a TF is not a SE associated TF in that sample
lapply(SampleList, function(x) {
  filepath <- subset(coltron, coltron$SampleFiles %in% x) #isolate single sample
  SE_TFs <- read.table(as.character(filepath$SE_TF_path),sep="\t", header=T)
  SE_TFs <- table(SE_TFs$TF_name)
  SE_TFs <- as.data.frame(SE_TFs)
  library(plyr)
  SE_TFs_joined <- join(All_SE_TFs.counted, SE_TFs, by = "Var1")
    SE_TFs_joined[is.na(SE_TFs_joined)] <- 0
    SE_TFs_joined <- as.data.frame(SE_TFs_joined[,c("Freq")])
    colnames(SE_TFs_joined) <- c(x)
    All_SE_TFs.counted <<- cbind(All_SE_TFs.counted,SE_TFs_joined) #write outside the loop to build up matrix
  
  TF_Degree <- read.table(paste(filepath$path,filepath$SampleFiles,"_DEGREE_TABLE.txt",sep=""),sep="\t", header=T)
  
    colnames(TF_Degree)=c("Var1","In_Degree","Out_Degree","Total_Connections")
    TF_Degree$Total_Connections = TF_Degree$Total_Connections/max(TF_Degree$Total_Connections)
    
    SE_TFs_degree <- join(All_SE_TFs.degree, TF_Degree, by = "Var1")
    SE_TFs_degree[is.na(SE_TFs_degree)] <- 0
    SE_TFs_degree <- as.data.frame(SE_TFs_degree[,c("Total_Connections")])
    colnames(SE_TFs_degree) <- c(x)
    
    All_SE_TFs.degree <<- cbind(All_SE_TFs.degree,SE_TFs_degree) #write outside the loop to build up matrix
})

write.table(All_SE_TFs.counted, file = paste(project.folder,Filename,".freq.matrix.txt",sep=""), sep="\t", row.names=F, quote=F, col.names=T)
write.table(All_SE_TFs.degree, file = paste(project.folder,Filename,".degree.matrix.txt",sep=""), sep="\t", row.names=F, quote=F, col.names=T)

degree.matrix = as.matrix(All_SE_TFs.degree[,c(3:ncol(All_SE_TFs.degree))])
row.names(degree.matrix) = All_SE_TFs.degree$Var1

library(pheatmap)
library(RColorBrewer)
pheatmap(degree.matrix,scale='none')

All_SE_TFs.top = All_SE_TFs.degree
All_SE_TFs.top$rowmax = apply(All_SE_TFs.top[,c(3:ncol(All_SE_TFs.degree))], 1, FUN=max)
All_SE_TFs.top = subset(All_SE_TFs.top, All_SE_TFs.top$rowmax > 0.2)

All_SE_TFs.top$rowsum = apply(All_SE_TFs.top[,c(3:ncol(All_SE_TFs.degree))], 1, FUN=sum)
All_SE_TFs.top = subset(All_SE_TFs.top, All_SE_TFs.top$rowsum > 2)

#All_SE_TFs.top = subset(All_SE_TFs.top, All_SE_TFs.top$Var1 %in% c("SOX4"))

All_SE_TFs.top = All_SE_TFs.top[rowSums(All_SE_TFs.top==0)<=4,]

degree.top = as.matrix(All_SE_TFs.top[,c(3:(ncol(All_SE_TFs.top)-2))])
row.names(degree.top) = All_SE_TFs.top$Var1

#pheatmap(degree.top,scale='none',cluster_rows = T,  color = colorRampPalette(c("white", "dodgerblue", "purple"))(50),main=paste(Filename," coltron heatmap",sep=""))

#sort it out a plot again. from http://slowkow.com/notes/heatmap-tutorial/
mat_cluster_cols <- hclust(dist(t(degree.top)))
plot(mat_cluster_cols, main = "Unsorted Dendrogram", xlab = "", sub = "")

library(dendsort)
  sort_hclust <- function(...) as.hclust(dendsort(as.dendrogram(...)))
  mat_cluster_cols <- sort_hclust(mat_cluster_cols)
  plot(mat_cluster_cols, main = "Sorted Dendrogram", xlab = "", sub = "")
  mat_cluster_rows <- sort_hclust(hclust(dist(degree.top)))
  
pheatmap(degree.top,scale='none',cluster_cols=mat_cluster_cols,cluster_rows=mat_cluster_rows,  color = colorRampPalette(c("white", "dodgerblue", "purple"))(50),main=paste(Filename," coltron heatmap",sep=""))
  
