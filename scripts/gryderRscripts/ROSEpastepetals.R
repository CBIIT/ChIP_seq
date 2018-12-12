###Extract MYCN SE from sample_id."_peaks_AllEnhancers.table.txt"; add that row to noMYCN version

#Get sample names via excel
require(XLConnect)
setwd("K:/projects/ChIP_seq")
wb = loadWorkbook("manage_samples/ChIP_seq_samples.xlsx")
ChIPs = readWorksheet(wb, sheet = "ChIP_seq_samples", header = TRUE)
write.table(ChIPs, file = "manage_samples/ChIP_seq_samples.txt", sep=",", row.names=FALSE, quote=FALSE)
#or load in from txt
ChIPs <- read.table("manage_samples/ChIP_seq_samples_NB_Ac.txt", sep="\t", header=T)

#define ROSE output locations
ROSES = ChIPs[,c("SampleFiles","PairedInput")]
ROSES = ROSES[grep("HLL3CBGXX", ROSES$SampleFiles), ]
ROSES$wMYCN.Enhancers = paste("DATA/",ROSES$SampleFiles,"/MACS_Out_p_1e-03_dup1/ROSE_out_12500/",ROSES$SampleFiles,"_peaks_AllEnhancers.table.txt",sep="")
ROSES$woMYCN.Enhancers= paste("DATA/",ROSES$SampleFiles,"/MACS_Out_p_1e-03_noMYCN/ROSE_out_12500/",ROSES$SampleFiles,"_peaks_AllEnhancers.table.txt",sep="")

##extract top SE from each sample run with MYCN, add it back to 
SampleList = ROSES$SampleFiles #create Sample list for lapply

lapply(SampleList, function(x) {
  filepaths = subset(ROSES, ROSES$SampleFiles %in% x) #isolate single sample
  AllEnhancers = read.table(filepaths$wMYCN.Enhancers, sep="\t", header=T, quote="")
    #write.table(AllEnhancers, file=paste("projects/NB_Thiele/coltron/ROSE.outputs/",x,"p-3_dup1_peaks_AllEnhancers.table.original.txt"), sep="\t", row.names=FALSE, col.names=TRUE, quote=FALSE)
  MYCNEnhancer = subset(AllEnhancers,CHROM %in% c("chr2"))
  MYCNEnhancer = subset(MYCNEnhancer, START > 16120000 & STOP < 16400000)
  CalledwoMYCN = read.table(filepaths$woMYCN.Enhancers, sep="\t", header=T, quote="")
    #write.table(CalledwoMYCN, file=paste("projects/NB_Thiele/coltron/ROSE.outputs/",x,"p-3_noMYCN_peaks_AllEnhancers.table.original.txt"), sep="\t", row.names=FALSE, col.names=TRUE, quote=FALSE)
  MYCNpastedSE = rbind(MYCNEnhancer,CalledwoMYCN)
  MYCNpastedSE$H3K27ac.delta = MYCNpastedSE[,7]-MYCNpastedSE[,8]
  MYCNpastedSE = MYCNpastedSE[order(-MYCNpastedSE[,11]),]
  MYCNpastedSE$enhancerRank = order(-MYCNpastedSE$H3K27ac.delta)
  MYCNpastedSE = MYCNpastedSE[,1:10]
    
  ##write frontloaded frame with comments for extra rows"
  buffer5lines = as.data.frame(c("####","#fill in","#to match","#ROSE2","####"))
  write.table(buffer5lines, file=paste("projects/NB_Thiele/coltron/ROSE.outputs/",x,"p-3_MYCNback_peaks_AllEnhancers.table.txt"), sep="\t", row.names=FALSE, col.names=FALSE, quote=FALSE, append = FALSE)
  write.table(MYCNpastedSE, file=paste("projects/NB_Thiele/coltron/ROSE.outputs/",x,"p-3_MYCNback_peaks_AllEnhancers.table.txt"), sep="\t", row.names=FALSE, col.names=TRUE, quote=FALSE, append = TRUE)
    #write.table(MYCNpastedSE, file=filepaths$woMYCN.Enhancers,sep="\t", row.names=FALSE, col.names=TRUE, quote=FALSE)
  
})