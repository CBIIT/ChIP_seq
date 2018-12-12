##load data.  each has 2 bam files coverage
bedcov1 <- read.table("RMS_Epigenetics/CTR_DvTram/CTR_DvT_DNase/CTR_48hr_DMSO_Tram_combined_30X_D.DvTcov.30X_D.bed", sep="\t", header=F)
  bedcov1$V4 <- bedcov1$V4*1000000/22990181
  bedcov1$V5 <- bedcov1$V5*1000000/27843174
bedcov2 <- read.table("RMS_Epigenetics/CTR_DvTram/CTR_DvT_DNase/CTR_48hr_DMSO_Tram_combined_30X_D.DvTcov.40X_D.bed", sep="\t", header=F)
  bedcov1$cov2.1 <- bedcov2$V4*1000000/10971996
  bedcov1$cov2.2 <- bedcov2$V5*1000000/18795477
bedcov3 <- read.table("RMS_Epigenetics/CTR_DvTram/CTR_DvT_DNase/CTR_48hr_DMSO_Tram_combined_30X_D.DvTcov.MYOD.bed", sep="\t", header=F)
  bedcov1$cov3.1 <- bedcov3$V4*1000000/12496880
  bedcov1$cov3.2 <- bedcov3$V5*1000000/5426759
bedcov4 <- read.table("RMS_Epigenetics/CTR_DvTram/CTR_DvT_DNase/CTR_48hr_DMSO_Tram_combined_30X_D.DvTcov.MYC.bed", sep="\t", header=F)
  bedcov1$cov4.1 <- bedcov4$V4*1000000/11590382
  bedcov1$cov4.2 <- bedcov4$V5*1000000/10982673
bedcov5 <- read.table("RMS_Epigenetics/CTR_DvTram/CTR_DvT_DNase/CTR_48hr_DMSO_Tram_combined_30X_D.DvTcov.H3K27ac.bed", sep="\t", header=F)
  bedcov1$cov5.1 <- bedcov5$V4*1000000/44431521
  bedcov1$cov5.2 <- bedcov5$V5*1000000/13016246
bedcov6 <- read.table("RMS_Epigenetics/CTR_DvTram/CTR_DvT_DNase/CTR_48hr_DMSO_Tram_combined_30X_D.DvTcov.E2A.bed", sep="\t", header=F)
  bedcov1$cov6.1 <- bedcov6$V4*1000000/39247973
  bedcov1$cov6.2 <- bedcov6$V5*1000000/41294680

colnames(bedcov1)<-c("CHROM","START","STOP","DNase30xD48","DNase30xT48","DNase40xD48","DNase40xT48","MYOD.D48","MYOD.T48","MYC.D48","MYC.T48","H3K27ac.D48","H3K27ac.T48","E2A.D48","E2A.T48") 

##filtering criteria and delta treatments
  bedcov1$DNase.D48.ave=(bedcov1$DNase30xD48+bedcov1$DNase40xD48)/2
  bedcov1$DNase.T48.ave=(bedcov1$DNase30xT48+bedcov1$DNase40xT48)/2
  library(dplyr)
  bedcov1=bedcov1 %>% mutate(DNase.max = pmax(DNase.D48.ave, DNase.T48.ave))

  bedcov1$d.DNase=bedcov1$DNase.T48.ave-bedcov1$DNase.D48.ave
  bedcov1$d.MYOD=bedcov1$MYOD.T48-bedcov1$MYOD.D48
  bedcov1$d.MYC=bedcov1$MYC.T48-bedcov1$MYC.D48
  bedcov1$d.H3K27ac=bedcov1$H3K27ac.T48-bedcov1$H3K27ac.D48
  bedcov1$d.E2A=bedcov1$E2A.T48-bedcov1$E2A.D48
  bedcov1$log2.d.DNase=log(bedcov1$DNase.T48.ave+0.1,2)-log(bedcov1$DNase.D48.ave+0.1,2)
  bedcov1$log2.d.MYOD=log(bedcov1$MYOD.T48+0.1,2)-log(bedcov1$MYOD.D48+0.1,2)

##filter
bedcov.filter=subset(bedcov1,bedcov1$DNase.max>1.5)

###subset data
DNase.cov <- bedcov1[,19:25]
DNase.cov.filter <- bedcov.filter[,19:25]
T.up.cov = subset(bedcov.filter,bedcov.filter$d.MYOD>2)
T.up.cov = subset(T.up.cov,T.up.cov$d.DNase>1)
T.up.noMYOD= subset(bedcov.filter,bedcov.filter$d.MYOD<2)
T.up.noMYOD= subset(T.up.noMYOD,T.up.noMYOD$d.DNase>1)
T.down.cov = subset(bedcov.filter,bedcov.filter$d.DNase<(-1))
T.down.cov = subset(T.down.cov,T.down.cov$d.MYOD<(0))

###raw plots
library(ggplot2)
library(scales) 
ggplot(DNase.cov,aes(x=d.DNase)) + geom_histogram() +geom_histogram(data=DNase.cov.filter,color = "red",alpha=0.1)
ggplot(DNase.cov.filter,aes(x=d.DNase,y=d.MYOD)) +geom_point(colour="gray",alpha=0.2) +theme_bw() +
  geom_point(data=T.up.cov,color = "aquamarine3",alpha=0.7) + #geom_point(data=T.up.noMYOD,color = "gold",alpha=0.2) +
  geom_point(data=T.down.cov,color = "dodgerblue",alpha=0.7) +theme(aspect.ratio=1)#+scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x), labels = trans_format("log10", math_format(10^.x))) +scale_x_log10(breaks = trans_breaks("log10", function(x) 10^x), labels = trans_format("log10", math_format(10^.x)))


###melted data
library(tidyr)
library(plyr)
T.up.less=T.up.cov[,8:17]
T.down.less=T.down.cov[,8:17]
T.up.gather<- gather(T.up.less, Condition, RPM, MYOD.D48:DNase.T48.ave)
T.down.gather<- gather(T.down.less, Condition, RPM, MYOD.D48:DNase.T48.ave)


###boxing ring
ggplot(T.down.gather, aes(x=Condition,y=RPM)) + scale_y_continuous(lim=c(0,20)) +theme_bw()+
  geom_boxplot(outlier.shape = NA,aes(fill= Condition),scale = "width") #+scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x), labels = trans_format("log10", math_format(10^.x))) #+scale_x_log10(breaks = trans_breaks("log10", function(x) 10^x), labels = trans_format("log10", math_format(10^.x)))


###matrix plots
library(GGally)
ggpairs(DNase.cov.filter, lower=list(continuous="smooth", params=c(colour="blue",alpha=0.2)),axisLabels='show') #diag=list(continuous="density")

###write txt and bed file
write.table(bedcov.filter, file="RMS_Epigenetics/CTR_DvTram/CTR_DvT_DNase/bedcov.all.DvT.txt", sep="\t", row.names=FALSE, col.names=TRUE, quote=FALSE)
T.DNaseup.MYODup.bed = T.up.cov[,1:3]
T.DNaseup.MYODdown.bed = T.up.noMYOD[,1:3]
T.DNasedown.bed = T.down.cov[,1:3]
write.table(T.DNaseup.MYODup.bed, file="RMS_Epigenetics/CTR_DvTram/CTR_DvT_DNase/T.DNaseup.MYODup.bed", sep="\t", row.names=FALSE, col.names=FALSE, quote=FALSE)
write.table(T.DNaseup.MYODdown.bed, file="RMS_Epigenetics/CTR_DvTram/CTR_DvT_DNase/T.DNaseup.MYODdown.bed", sep="\t", row.names=FALSE, col.names=FALSE, quote=FALSE)
write.table(T.DNasedown.bed, file="RMS_Epigenetics/CTR_DvTram/CTR_DvT_DNase/T.DNasedown.bed", sep="\t", row.names=FALSE, col.names=FALSE, quote=FALSE)
