##BED building

#bed of TES overlaps with K27ac peak calling
bedGenes = read.table("RMS_Epigenetics/EpiDrugs/NGSplots/HDACi_Pol2_K27ac/TES/TES_Ent6_K27acp-7.one.bed",header=F,sep="\t")
colnames(bedGenes) = c("CHR","START","END","gene_id")
bedGenes$source = as.factor("RH4_Ent6_K27ac_p-7")
bedGenes$length = bedGenes$END-bedGenes$START
bedGenes=subset(bedGenes,bedGenes$length>500)
library(ggplot2)
ggplot(bedConcat, aes(x=length))+geom_density(fill="blue", alpha=0.4) #+geom_histogram(binwidth = 10)

#all Gene TES file
bedAllGenes = read.table("RMS_Epigenetics/EpiDrugs/NGSplots/HDACi_Pol2_K27ac/TES/TES.bed",header=F,sep="\t")
colnames(bedAllGenes) = c("CHR","START","END","gene_id")
bedAllGenes$source = as.factor("5kb_past_TES")
bedAllGenes$length = bedAllGenes$END-bedAllGenes$START

#concatenate (to fill in the missing genes)
bedConcat = rbind(bedGenes, bedAllGenes)
bedConcat = bedConcat[!duplicated(bedConcat[,c('gene_id')]),]

#makes BEDs for each gene list in a folder
genesetpath = "RMS_Epigenetics/EpiDrugs/NGSplots/HDACi_Pol2_K27ac/TES/Spike/"
config = read.table(paste(genesetpath,"NGSconfig.setup.txt",sep=""), sep="\t", header=T)

files <- list.files(path=genesetpath, pattern="*genelist", full.names=T, recursive=FALSE)

lapply(files, function(x) {
  genelist <<- read.table(x,header=F,sep="\t")
  gspbed <<- subset(bedConcat, bedConcat$gene_id %in% genelist$V1)
  gspbed.filename <<- gsub(".genelist",".gsp.bed",basename(x))
  write.table(gspbed,file=paste(genesetpath,gspbed.filename,sep=""),sep="\t", row.names=F,col.names=F) #write beds
})


