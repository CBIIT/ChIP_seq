###Gather BEDs
#find the file:  assign p-value for chosen MACS outputs; make this a new column in the Sample sheet.  
#load sample sheet from excel file
  require(XLConnect)
  setwd("K:/projects/ChIP_seq")
  wb = loadWorkbook("manage_samples/ChIP_seq_samples.xlsx")
  ChIPs = readWorksheet(wb, sheet = "ChIP_seq_samples", header = TRUE)

#define meta output paths
meta = ChIPs[grep("noMYCN", ChIPs$ChosenP_path), ]
meta = meta[grep("2D_C|2D_RA|4D_RA|8D_RA", meta$SampleFiles), ]  #for testing with smaller sample size
meta$path = paste("DATA/",meta$SampleFiles,"/MACS_Out_",meta$ChosenP_path,"/",sep="")

#set project global parameters
project.folder = "projects/NB_Thiele/RPMPR/"
Filename = "KCNR_K27ac"

#pass through saved file
write.table(meta, file = paste(project.folder,Filename,"metaSamples.txt",sep=""), sep="\t", row.names=FALSE, quote=FALSE)
meta <- read.table(paste(project.folder,Filename,"metaSamples.txt",sep=""), sep="\t", header=T)

################################################################
######### COMBINING BED files from ChIP-seq samples#############
################################################################

#establish a list of samples to get BEDs for, BEDtype and dataframe to build
SampleList = as.character(meta$SampleFiles) #create Sample list for lapply
BEDtype = "_peaks.narrowPeak.nobl.GREAT.bed"
allBEDs = data.frame(CHR=character(),START=numeric(),END=numeric())

#loop through a second time, but now make a column for each sample, leave 0 if a TF is not a SE associated TF in that sample
lapply(SampleList, function(x) {
  filepath <- subset(meta, meta$SampleFiles %in% x) #isolate single sample
  
  BED <- read.table(paste(filepath$path,filepath$SampleFiles,BEDtype,sep=""),sep="\t", header=T)
  
    colnames(BED)=colnames(allBEDs)
    
    allBEDs <<- rbind(allBEDs,BED) #write outside the loop to build up matrix
})

write.table(allBEDs, file = paste(project.folder,Filename,".combined.bed",sep=""), sep="\t", row.names=F, quote=F, col.names=F)

