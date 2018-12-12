###Gather read count data

setwd("K:/projects/ChIP_seq/")
samples <- list.dirs(path="DATA", full.names=T, recursive=F)

sample.folders <- samples[intersect(grep("RH4|CTR", samples),grep("021",samples,invert=FALSE))]  #subset if needed Waga|MCC|UISO|MKL|MS1|MS-1
sample.folders <- sample.folders[grep('RH4[Flag]', sample.folders)]
sample.IDs <- basename(sample.folders)

#initiate dataframe
reads.allsamples = data.frame(sample = c("reads mapped"))


lapply(sample.folders, function(x) {
  #read in txt file
  flagstat.path = paste(x,"/",basename(x),".flagstat.txt",sep="")
  flagstat = read.table(flagstat.path,sep="\t", header=F)
  reads.map = as.data.frame(flagstat[5,])
  
  library(stringr)
  reads.map = str_split_fixed(reads.map[,1], " +", 2)
  reads.map = as.data.frame(as.numeric(reads.map[,1]))
  colnames(reads.map) = basename(x)
  reads.allsamples <<- cbind(reads.allsamples,reads.map)
})

write.table(reads.allsamples, file="manage_samples/sample_stats/RMS_DNase_reads.txt", sep="\t", row.names=FALSE, col.names=TRUE, quote=FALSE)
#navigate in biowulf2: cd /home/gryderbe/khanlab/projects/ChIP_seq/manage_samples/
#open permissions with:  chgrp -R khanlab sample_stats/

###SPIKE IN read counts
#initiate dataframe
sample.spike.summary <- paste(sample.folders,"SpikeIn/spike_map_summary",sep="/")
reads.allsamples = read.table("DATA/Sample_HCT116_5FU_p53_019_C_HCLGVBGX2/SpikeIn/spike_map_summary", sep="\t", header=T)
reads.allsamples = reads.allsamples[-1,]

lapply(sample.spike.summary, function(x) {
  #read in txt file
  spike.summary <<- read.table(x,sep="\t", header=T)
  filename <<- gsub("/SpikeIn/spike_map_summary","",x)
  rownames(spike.summary) = basename(filename)
  reads.allsamples <<- rbind(reads.allsamples,spike.summary)
})

write.table(reads.allsamples, file="manage_samples/sample_stats/RH4.HDACi_021.SpikeIn.txt", sep="\t", row.names=T, col.names=TRUE, quote=FALSE)
