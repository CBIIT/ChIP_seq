##load data
DNaseNILC_PF19 <- read.table("RMS_Epigenetics/NILC_PF19/P3F_activated_sites/P3F_cov_DNaseNILC_PF19.bed", sep="\t", header=F)
H3K27acNILC_PF19 <- read.table("RMS_Epigenetics/NILC_PF19/P3F_activated_sites/P3F_cov_H3K27acNILC_PF19.bed", sep="\t", header=F)

DNase_v_H3K27ac = DNaseNILC_PF19[,1:5]
DNase_v_H3K27ac$neglog10.qval <- DNase_v_H3K27ac$V5/10
DNase_v_H3K27ac$enhancers <- DNaseNILC_PF19$V11

DNaseNILC.reads = 19423852
DNasePF19.reads = 21110506
H3K27acNILC.reads = 28138534
H3K27acPF19.reads = 29191709

DNase_v_H3K27ac$DNase.change=DNaseNILC_PF19$V13 - DNaseNILC_PF19$V12
DNase_v_H3K27ac$DNase.change.RPM=DNaseNILC_PF19$V13*1000000/DNasePF19.reads - DNaseNILC_PF19$V12*1000000/DNaseNILC.reads
DNase_v_H3K27ac$DNase.log2fc=log(DNaseNILC_PF19$V13*1000000/DNasePF19.reads+0.00001,2)-log(DNaseNILC_PF19$V12*1000000/DNaseNILC.reads+0.00001,2)

DNase_v_H3K27ac$H3K27ac.change=H3K27acNILC_PF19$V13 - H3K27acNILC_PF19$V12
DNase_v_H3K27ac$H3K27ac.change.RPM=H3K27acNILC_PF19$V13*1000000/H3K27acPF19.reads - H3K27acNILC_PF19$V12*1000000/H3K27acNILC.reads
DNase_v_H3K27ac$H3K27ac.log2fc=log(H3K27acNILC_PF19$V13*1000000/H3K27acPF19.reads+0.00001,2)-log(H3K27acNILC_PF19$V12*1000000/H3K27acNILC.reads+0.00001,2)
###make some bins!
DNase_v_H3K27ac$H3K27acNILC.RPM = H3K27acNILC_PF19$V12*1000000/H3K27acNILC.reads
DNase_v_H3K27ac$DNaseNILC.RPM   = DNaseNILC_PF19$V12*1000000/DNaseNILC.reads
DNase_v_H3K27ac$H3K27ac.bins=cut(DNase_v_H3K27ac$H3K27acNILC.RPM, c(-Inf,1,Inf))


#DNase_v_H3K27ac=subset(DNase_v_H3K27ac,DNase_v_H3K27ac$enhancers<1)
DNase_v_H3K27ac.NILCclosed=subset(DNase_v_H3K27ac,(DNase_v_H3K27ac$H3K27acNILC.RPM<5 & DNase_v_H3K27ac$DNaseNILC.RPM<2))
DNase_v_H3K27ac.activated=subset(DNase_v_H3K27ac.NILCclosed,(DNase_v_H3K27ac.NILCclosed$H3K27ac.change.RPM>1))


###write txt and bed file
write.table(DNase_v_H3K27ac, file="RMS_Epigenetics/NILC_PF19/P3F_activated_sites/DNase_v_H3K27ac.txt", sep="\t", row.names=FALSE, col.names=TRUE, quote=FALSE)
DNase_v_H3K27ac.activated.bed = DNase_v_H3K27ac.activated[,1:5]
write.table(DNase_v_H3K27ac.activated.bed, file="RMS_Epigenetics/NILC_PF19/NGSplot/NILCclosed_and_P3Factivated2.800.bed", sep="\t", row.names=FALSE, col.names=FALSE, quote=FALSE)


###plots
library(ggplot2)
ggplot(DNase_v_H3K27ac,aes(x=DNaseNILC.RPM,y=H3K27acNILC.RPM)) +geom_point(colour="grey",alpha=0.5) +theme_bw() +geom_point(data=DNase_v_H3K27ac.NILCclosed,colour="blue",alpha=0.5) +geom_point(data=DNase_v_H3K27ac.activated,colour="red2",alpha=0.5) +scale_x_continuous(limits=c(0,10)) +scale_y_continuous(limits=c(0,15)) 
ggplot(DNase_v_H3K27ac,aes(x=H3K27acNILC.RPM,y=H3K27ac.change.RPM)) +geom_point(colour="grey",alpha=0.5) +theme_bw() +geom_point(data=DNase_v_H3K27ac.activated,colour="red2",alpha=0.5)

ggplot(DNase_v_H3K27ac,aes(x=H3K27ac.change.RPM,y=DNase.change.RPM))  +geom_point(data=DNase_v_H3K27ac.NILCclosed,colour="blue",alpha=0.3)+ #+geom_point(colour="grey",alpha=0.5)
  theme_bw() +geom_point(data=DNase_v_H3K27ac.activated,colour="red2",alpha=0.4) +coord_flip() +
  scale_x_continuous(limits=c(-10,50)) + scale_y_continuous(limits=c(-5,10)) 