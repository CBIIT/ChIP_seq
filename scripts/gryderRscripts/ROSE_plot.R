###GET ENHANCER_TO_GENE###
##########################

#excel spreadsheet must be made, which has "callout" column to pick out certain enhaners.  Fill in blanks with periods "."
require(XLConnect)
wb = loadWorkbook("RMS_Epigenetics/ROSE/RMS_Tumors/ChosenP/Sample_RD_H3K27ac_008_C_HHC7JBGXX.1e-05.12500.combined.xlsx")
ROSE <- readWorksheet(wb, sheet = "SEtoGene", header = TRUE)
#write.table(ROSE, file="RMS_Epigenetics/ROSE/RMS_Tumors/ChosenPtxt/ROSE.SCMC.txt", sep="\t", row.names=FALSE, col.names=TRUE, quote=FALSE)
ROSE <- ROSE[,c(1:9,13,14)]
colnames(ROSE)<-c("REGION_ID","CHROM","START","STOP","NUM_LOCI","CONSTITUENT_SIZE","signal.bam","input.bam","enhancerRank","isSuper","callout") #include "input.bam", this time

ROSE$signal <- ROSE$signal.bam -ROSE$input.bam
ROSE <- subset(ROSE,ROSE$signal>0)
ROSE$norm.signal <- ROSE$signal/max(ROSE$signal)
ROSE.callout <- subset(ROSE,!(ROSE$callout %in% c(".")))
###rank percentile
ROSE.callout$rank.percent = 100-(ROSE.callout$enhancerRank*100/nrow(ROSE))

###output cluster information
write.table(ROSE.callout, file="RMS_Epigenetics/ROSE/RMS_Tumors/ROSE.callout.RMS216.txt", sep="\t", row.names=FALSE, col.names=TRUE, quote=FALSE)


##############
##ROSE PLOT###
##############
super.threshold <- sum(ROSE$isSuper)
library(ggplot2)
library(scales)
ggplot(ROSE, aes(x=enhancerRank, y=norm.signal)) + 
  geom_vline(xintercept = super.threshold,colour="grey", linetype = "longdash") +
  geom_point(data=ROSE.callout,color = "red") + geom_line(color="red") + 
  scale_x_reverse() + theme_bw() + scale_y_continuous(labels = comma) +
  labs(x = "Rank order enhancers") +labs(y = "H3K27ac signal (normalized RPM)") + #ggtitle("Human skeletal muscle myobtubes (HSMMtube)") +  
  geom_text(data=ROSE.callout,label=ROSE.callout$callout, hjust=1.2, size=5) + 
  theme(axis.text = element_text(size = 12), aspect.ratio=1,
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()
  )

###GET Callout Summary###
#########################

###initiate called out set
called.out = data.frame(ROSE.callout[,c("callout")])
colnames(called.out)<-c("callout")
SErank.mins = data.frame(1)

files <- list.files(path="RMS_Epigenetics/ROSE/RMS_Tumors/ChosenPtxt/", pattern="*.txt", full.names=T, recursive=FALSE)
lapply(files, function(x) {
  
  #excel spreadsheet must be made, which has "callout" column to pick out certain enhaners.  Fill in blanks with periods "."
  ROSE <- read.table(x,sep="\t", header=T)
  ROSE <- ROSE[,c("callout","enhancerRank","isSuper")]
  ROSE$rank.percent = 100-(ROSE$enhancerRank*100/nrow(ROSE))    ###rank percentile
  superROSE <- subset(ROSE, ROSE$isSuper == "1")
  super.percent.threshold = data.frame(min(superROSE$rank.percent))
  colnames(super.percent.threshold)<-c(x)
  SErank.mins <<- cbind(SErank.mins,super.percent.threshold)
  #ROSE.callout <- subset(ROSE,!(ROSE$callout %in% c(".")))
  #colnames(ROSE.callout)<-c("callout","enhancerRank","isSuper",x)
  #library(plyr)
  #called.out<<-join(called.out, ROSE.callout, by = "callout")
  
  
  ###output cluster information
  #write.table(called.out, file="RMS_Epigenetics/ROSE/RMS_Tumors/ROSE.callout.Summary2.txt", sep="\t", row.names=FALSE, col.names=TRUE, quote=FALSE)
  
})

