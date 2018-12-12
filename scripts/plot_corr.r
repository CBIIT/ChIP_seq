library(cluster)
library(gplots)
library(Biobase)
library(ctc)
library(ape)
library(rgl)

Args<-commandArgs()
in_file<-Args[6]
bam_list_file<-Args[7]

data = read.table(in_file, header=F, com='', sep="\t")
#data = read.table("scripts/corr/Mari/runSE/merged_cov.bed", header=F, com='', sep="\t") #testing only, Args bypass

bam_list = read.table(bam_list_file, header=F, com='', sep="\t")
#bam_list = read.table("scripts/corr/Mari/bam_list.txt", header=F, com='', sep="\t") #testing only, Args bypass

colnames(data)[1:3] <- c('Chr','Start','End')
colnames(data)[4:length(colnames(data))] <- as.character(bam_list[,1])
total_reads = bam_list[,2] #reads per million mapped reads
#total_reads = colSums(data[,4:length(data)]) #reads per million in these regions

data[,4:length(data)] = t(apply(t(data[,4:length(data)]+1), 2, function(x) x/total_reads*1000000))  #make dataframe into RPM
write.table(data, file='rpm.tsv', sep="\t", row.names = FALSE)
mat = as.matrix(data[,4:length(data)]) # convert to matrix
rowmed <- apply(mat,1,median)
mat_median = log2(mat)-log2(rowmed)
rowsize <- data[,3]-data[,2]
mat_sized <- mat/rowsize

cr = cor(mat,method="pearson")
write.table(cr, file='corr.tsv', sep="\t", row.names = T)
cr_median = cor(mat_median,method="pearson")
write.table(cr_median, file='corr_log2median.tsv', sep="\t", row.names = T)
cr_sized = abs(cor(mat_sized,method="pearson"))
write.table(cr_sized, file='corr_sized.tsv', sep="\t", row.names = T)

#mat = t(mat)
#pca <- prcomp(mat, scale=T)
#plot3d(pca$x[,1:3], type="s", size=1)
#text3d(pca$x[,1], pca$x[,2], pca$x[,3], text=rownames(pca$x),adj = -0.2 ,family="serif", font=5, cex=1)

#samples_median = hclust(as.dist(1-cr_median), method="average") # cluster samples
#samples = hclust(as.dist(1-cr), method="average") # cluster samples
pdf("correlation.pdf")
myheatcol1 <- colorRampPalette(c("black", "red"))(n = 1000)
myheatcol2 <- colorRampPalette(c("green","black", "red"))(n = 1000)
myheatcol3 <- colorRampPalette(c("black","red", "gold"))(n = 1000)
myheatcol4 <- colorRampPalette(c("dodgerblue","white", "red3"))(n = 1000)
try(heatmap.2(cr,        col = myheatcol4, scale='none', main = "Correlation, RPM", symm=TRUE, key=TRUE,density.info='none', trace='none', symkey=FALSE, Colv=TRUE,margins=c(10,10), cexCol=0.6, cexRow=0.6))
try(heatmap.2(cr_median, col = myheatcol2, scale='none', main = "Log2 Median Centered", symm=TRUE, key=TRUE,density.info='none', trace='none', symkey=FALSE, Colv=TRUE,margins=c(10,10), cexCol=0.6, cexRow=0.6))
try(heatmap.2(cr_sized,  col = myheatcol3, scale='none', main = "Correlation, Size Normalized", symm=TRUE, key=TRUE,density.info='none', trace='none', symkey=FALSE, Colv=TRUE,margins=c(10,10), cexCol=0.6, cexRow=0.6))


dev.off()

