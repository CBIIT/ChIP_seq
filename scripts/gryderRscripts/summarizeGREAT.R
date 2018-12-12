#GREAT ontology comparisons and plots
#first, import GREAT output tables

files <- list.files(path="RMS_Epigenetics/ROSE/RMS_Tumors/P3F_InTumorEnhancers/GREAT/", pattern="*.tsv", full.names=T, recursive=FALSE)

great.summary = read.table("RMS_Epigenetics/ROSE/RMS_Tumors/P3F_InTumorEnhancers/GREAT/FP_FN_greatExportAll.tsv",sep="\t",header=F,quote = "")

great.GOlist = great.summary[,2:3]
colnames(great.GOlist) <- c("GO_ID","Description")

lapply(files, function(x) {
  
  great <<- read.table(x,sep="\t", header=F, quote="")
  
  great <<- great[,2:3]
  colnames(great) <<- c("GO_ID","Description")
  great.GOlist <<- rbind(great.GOlist,great)
  great.GOlist = great.GOlist[!duplicated(great.GOlist[,c('Description')]),] 
})

#open each TSV file, remove unwanted columns (keep p-value). Then, name columns by X and join.

lapply(files, function(x) {
  
  great <<- read.table(x,sep="\t", header=F, quote="")
  great <<- great[,3:5]
  colnames(great) <<- c("Description",paste("Rank",basename(x)),basename(x))
  library(plyr)
  great.GOlist <<- join(great.GOlist, great, by = "Description")
  
})
great.GOlist[is.na(great.GOlist)] <- 501
great.GOlist[,"Rank min"] <- apply(great.GOlist[, c(3,5,7)], 1, min)
topRank = subset(great.GOlist, great.GOlist$`Rank min` < 11)
topRank = topRank[,c(2,4,6,8,10)]
topRank[,2:5] = (-1)*log10(topRank[,2:5])
topRank.rownames=topRank$Description
topRank=topRank[,2:5]
rownames(topRank)=topRank.rownames
topRank[topRank<0] <- 0

write.table(topRank, file="RMS_Epigenetics/ROSE/RMS_Tumors/P3F_InTumorEnhancers/GREAT/topRank.GREAT.summary.txt", sep="\t", row.names=F, col.names=T, quote=FALSE)

library(pheatmap)
pheatmap(topRank[,1:3],
          show_colnames = T, show_rownames = T, cluster_rows = T,
          cluster_cols = T, legend = TRUE,
          clustering_distance_rows = "euclidean", border_color = FALSE)




