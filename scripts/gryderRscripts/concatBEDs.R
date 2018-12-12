#concatentate all BEDs in a folder

  files <- list.files(path="RMS_Epigenetics/ROSE/RMS_Tumors/ChosenPbed/FPRMS_SEs/", pattern="*SE.bed", full.names=T, recursive=FALSE)

lapply(files, function(x) {
  
  bed <<- read.table(x,sep="\t", header=F)
  write.table(bed, file="RMS_Epigenetics/ROSE/RMS_Tumors/bed.concat.FP-RMS-SEs.bed", sep="\t", row.names=F, col.names=F, quote=FALSE, append=T)
  
})
