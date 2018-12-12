##load data
CTR_DvT_EZH2 <- read.table("RMS_Epigenetics/CTR_DvTram/EZH2_ERK_H3K27me3/CTR_DvT_PcG_EZH2.bed", sep="\t", header=F)
CTR_DvT_ERK <- read.table("RMS_Epigenetics/CTR_DvTram/EZH2_ERK_H3K27me3/CTR_DvT_PcG_ERK.bed", sep="\t", header=F)
    colnames(CTR_DvT_EZH2) <- c("#Chr", "Start","End","Peak.Score","distTSS","EZH2.D48","EZH2.T48","poise")
    colnames(CTR_DvT_ERK) <- c("#Chr", "Start","End","Peak.Score","distTSS","ERK.D48","ERK.T48","poise")
    CTR_DvT_EZH2$log10.distTSS <- log10(abs(CTR_DvT_EZH2$distTSS))
    CTR_DvT_ERK$log10.distTSS <- log10(abs(CTR_DvT_ERK$distTSS))
    CTR_DvT_EZH2$Peak.bins=cut(abs(CTR_DvT_EZH2$distTSS), c(-Inf,2500,Inf),labels=c("promoters","enhancers"))
    CTR_DvT_ERK$Peak.bins=cut(abs(CTR_DvT_ERK$distTSS), c(-Inf,2500,Inf),labels=c("promoters","enhancers"))

##build comparison file which combines data from multiple loaded bed files
    EZH2vERK = CTR_DvT_EZH2[,1:5]
    EZH2vERK$neglog10.qval <- EZH2vERK$Peak.Score/10
    EZH2vERK$poise <- CTR_DvT_EZH2$poise
    EZH2.D48.reads = 18753007
    EZH2.T48.reads = 23957825
    ERK.D48.reads = 28774900
    ERK.T48.reads = 10399465
    EZH2vERK$EZH2.change=CTR_DvT_EZH2$EZH2.T48 - CTR_DvT_EZH2$EZH2.D48
    EZH2vERK$EZH2.change.RPM=CTR_DvT_EZH2$EZH2.T48*1000000/EZH2.T48.reads - CTR_DvT_EZH2$EZH2.D48*1000000/EZH2.D48.reads
    EZH2vERK$EZH2.change.RPKb=CTR_DvT_EZH2$EZH2.T48/(CTR_DvT_EZH2$End-CTR_DvT_EZH2$Start)*1000 - CTR_DvT_EZH2$EZH2.D48/(CTR_DvT_EZH2$End-CTR_DvT_EZH2$Start)*1000
    EZH2vERK$ERK.change=CTR_DvT_ERK$ERK.T48 - CTR_DvT_ERK$ERK.D48
    EZH2vERK$ERK.change.RPM=CTR_DvT_ERK$ERK.T48*1000000/ERK.T48.reads - CTR_DvT_ERK$ERK.D48*1000000/ERK.D48.reads
    EZH2vERK$ERK.change.RPKb=CTR_DvT_ERK$ERK.T48/(CTR_DvT_ERK$End-CTR_DvT_ERK$Start)*1000 - CTR_DvT_ERK$ERK.D48/(CTR_DvT_ERK$End-CTR_DvT_ERK$Start)*1000
    #make some bins!
    EZH2vERK$ERK.D48.RPM  = CTR_DvT_ERK$ERK.D48*1000000/ERK.D48.reads
    EZH2vERK$EZH2.D48.RPM = CTR_DvT_EZH2$EZH2.D48*1000000/EZH2.D48.reads
    EZH2vERK$ERK.T48.RPM  = CTR_DvT_ERK$ERK.T48*1000000/ERK.T48.reads
    EZH2vERK$EZH2.T48.RPM = CTR_DvT_EZH2$EZH2.T48*1000000/EZH2.T48.reads
    EZH2vERK$Peak.bins=cut(abs(EZH2vERK$distTSS), c(-Inf,2500,Inf),labels=c("promoters","enhancers"))
    EZH2vERK$EZH2.bins=cut(EZH2vERK$EZH2.change.RPKb, c(-1000,-20,20,1000),labels=c("loss","nochange","gain"))
    #make subset data frames
    EZH2vERK.up=subset(EZH2vERK,(EZH2vERK$EZH2.change.RPKb>20))
    EZH2vERK.down=subset(EZH2vERK,(EZH2vERK$EZH2.change.RPKb<(-20)))
    EZH2vERK.same=subset(EZH2vERK,EZH2vERK$EZH2.bins %in% c("nochange"))
    #callout!
    EZH2vERK.callout=subset(EZH2vERK,EZH2vERK$Start %in% c("203055226","100105860"))

###write txt and bed file
names(EZH2vERK)[names(EZH2vERK) == 'Chr'] <- '#Chr' #get it ready for the EDEN format .bed life
write.table(EZH2vERK, file="RMS_Epigenetics/EDEN/CTR_DvT_EZH2/CTR_DvT_PcG_EZH2.bed", sep="\t", row.names=FALSE, col.names=TRUE, quote=FALSE)
    EZH2vERK.down.bed = EZH2vERK.down[,1:5]
    EZH2vERK.down.promoters = subset(EZH2vERK.down, EZH2vERK.down$Peak.bins %in% c("promoters"))
    EZH2vERK.down.promoters.GREAT.bed = EZH2vERK.down.promoters[,1:3]
    write.table(EZH2vERK.down.bed, file="RMS_Epigenetics/CTR_DvTram/EZH2_ERK_H3K27me3/EZH2vERK.down.bed", sep="\t", row.names=FALSE, col.names=FALSE, quote=FALSE)
    write.table(EZH2vERK.down.GREAT.bed, file="RMS_Epigenetics/CTR_DvTram/EZH2_ERK_H3K27me3/EZH2vERK.down.GREAT.bed", sep="\t", row.names=FALSE, col.names=FALSE, quote=FALSE)

###plots
library(ggplot2)
library(gridExtra)
#basic
    pEZH2=ggplot(CTR_DvT_EZH2,aes(x=log10.distTSS,y=EZH2.Reads,color=variable)) +geom_point(aes(y=EZH2.D48,col="EZH2.D48"),colour="hotpink1",alpha=0.2) +geom_point(aes(y=EZH2.T48,col="EZH2.T48"),colour="hotpink3",alpha=0.2) +theme_bw()+coord_flip()
    pEZH2.hist=ggplot(CTR_DvT_EZH2,aes(x=log10.distTSS,fill=variable)) +geom_histogram(fill="hotpink") +theme_bw()+coord_flip()
    pERK=ggplot(CTR_DvT_ERK,aes(x=log10.distTSS,y=ERK.Reads,color=variable)) +geom_point(aes(y=ERK.D48,col="ERK.D48"),colour="chartreuse2",alpha=0.2) +geom_point(aes(y=ERK.T48,col="ERK.T48"),colour="chartreuse3",alpha=0.2) +theme_bw()+coord_flip()
    grid.arrange(pEZH2,pERK,pEZH2.hist,ncol=1)
    #distance bins
    ggplot(EZH2vERK,aes(x=log10(abs(distTSS)),fill=variable)) +geom_histogram(fill="hotpink") +theme_bw()+coord_flip() +facet_wrap(~Peak.bins)
    
    #compare seq1 and seq2
    pEZH2.D48vERK.D48=ggplot(EZH2vERK,aes(x=EZH2.D48.RPM,y=ERK.D48.RPM)) +geom_point(colour="mediumorchid",alpha=0.2) +theme_bw() +scale_x_continuous(limits=c(0,20)) +scale_y_continuous(limits=c(0,20)) #+geom_point(data=EZH2vERK.up,colour="blue",alpha=0.5) +geom_point(data=EZH2vERK.down,colour="red2",alpha=0.5) +scale_x_continuous(limits=c(0,10)) +scale_y_continuous(limits=c(0,15)) 
    pEZH2.T48vERK.T48=ggplot(EZH2vERK,aes(x=EZH2.T48.RPM,y=ERK.T48.RPM)) +geom_point(colour="mediumorchid4",alpha=0.2) +theme_bw() +scale_x_continuous(limits=c(0,20)) +scale_y_continuous(limits=c(0,20))#+geom_point(data=EZH2vERK.up,colour="blue",alpha=0.5) +geom_point(data=EZH2vERK.down,colour="red2",alpha=0.5) +scale_x_continuous(limits=c(0,10)) +scale_y_continuous(limits=c(0,15)) 
    grid.arrange(pEZH2.D48vERK.D48,pEZH2.T48vERK.T48,ncol=2)

#delta treatment
pDvT_PnE=ggplot(EZH2vERK.same,aes(x=ERK.change.RPKb,y=EZH2.change.RPKb)) +geom_point(colour="grey",alpha=0.2) +theme_bw() +coord_flip() +facet_wrap(~Peak.bins) +geom_point(data=EZH2vERK.up,colour="blue",alpha=0.2) +geom_point(data=EZH2vERK.down,colour="red2",alpha=0.2) +geom_point(data=EZH2vERK.callout,colour="black")
pDvT_PnE.hist=ggplot(EZH2vERK.same,aes(x=EZH2.change.RPKb)) +geom_histogram(fill="grey",binwidth = 20) +theme_bw() +facet_wrap(~Peak.bins) +scale_x_continuous(labels = NULL,breaks=NULL,limits=c(-600,400)) +theme(axis.title.x = element_blank()) +scale_y_continuous(limits=c(0,900)) +geom_histogram(data=EZH2vERK.down,fill="red2",alpha=0.5,binwidth = 20) +geom_histogram(data=EZH2vERK.up,fill="blue",alpha=0.5,binwidth = 20)
grid.arrange(pDvT_PnE.hist,pDvT_PnE,heights=c(1,2))

#poised
ggplot(EZH2vERK.down.promoters,aes(x=ERK.change.RPKb,y=EZH2.change.RPKb)) +geom_point(colour="purple",alpha=0.3) +
  theme_bw() +coord_flip() +facet_wrap(~poise) 

###EvP
pval.EvP.EZH2 <- signif(t.test(EZH2.change.RPKb ~ Peak.bins,EZH2vERK)$p.value,3)
ggplot(EZH2vERK, aes(x=Peak.bins,y=EZH2.change.RPKb)) + theme_bw() + stat_boxplot(geom ='errorbar',stat_params = list(width = 0.2)) +geom_boxplot(outlier.shape=NA,width=0.7) +coord_cartesian(ylim=c(-600,400)) + annotate("text",x=1.5,y=120,label=pval.EvP.EZH2)

##melt DMSO v Tram columns into single column
#CTR_DvT.melted =
#pval.EvP.EZH2 <- signif(t.test(EZH2.change.RPKb ~ Peak.bins,CTR_DvT_EZH2)$p.value,3)
#ggplot(CTR_DvT_EZH2, aes(x=Peak.bins,y=EZH2.change.RPKb)) + theme_bw() + stat_boxplot(geom ='errorbar',stat_params = list(width = 0.2)) +geom_boxplot(outlier.shape=NA,width=0.7) +coord_cartesian(ylim=c(-300,150)) + annotate("text",x=1.5,y=120,label=pval.EvP.EZH2)

ggplot(EZH2vERK,aes(x=EZH2.change.RPKb)) +geom_histogram(colour="grey",alpha=0.5) +theme_bw() +coord_flip() +facet_wrap(~EZH2.bins) #+geom_point(data=EZH2vERK.up,colour="blue",alpha=0.5) +geom_point(data=EZH2vERK.down,colour="red2",alpha=0.5)

