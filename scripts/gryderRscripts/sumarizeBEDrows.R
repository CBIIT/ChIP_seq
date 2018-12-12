#Triple Venn

files <- list.files(path="RMS_Epigenetics/EDEN/RH4_P3F/P3Fp-5_q30/", pattern="*.bed", full.names=T, recursive=FALSE)

bed.nrows.all = data.frame(count = 1,qualityMED = 2,qualityAVE = 3)

lapply(files, function(x) {
  
  bed <<- read.table(x,sep="\t", header=F)
  bed.nrows <<- as.data.frame(nrow(bed))
  medianPval <<- as.data.frame(median(bed$V5))
  averagePval <<- as.data.frame(mean(bed$V5))
  bed.nrows <<- cbind(bed.nrows,medianPval,averagePval)
  
  rownames(bed.nrows) <<- basename(x)
  colnames(bed.nrows) <<- c("count","qualityMED","qualityAVE")
  bed.nrows.all <<- rbind(bed.nrows.all,bed.nrows)
  
})

bed.nrows.all = bed.nrows.all[-1,]
write.table(bed.nrows.all, file="RMS_Epigenetics/EDEN/RH4_P3F/P3Fp-5_q30/bed.rows.summary.txt", sep="\t", row.names=T, col.names=F, quote=FALSE)
