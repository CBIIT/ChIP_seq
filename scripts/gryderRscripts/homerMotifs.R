##load data.  
homer.string <- readLines("RMS_Epigenetics/BCHN/DNase_DvT/output1/HOMER_T.DNasedown_only/knownResults.txt")
homer <- read.table(text=homer.string,skip=1,sep="\t")
colnames(homer)<-c("Motif","Consensus","Pval","Log.Pval","Qval","Seq.w.Motif","pcentSeq.w.Motif","Backgr.w.Motif","pcentBackgr.w.Motif") 
homer$rownumber=1:nrow(homer)
homer$neglogP=(-1)*homer$Log.Pval

#subset for lables
homer.subset1 <- homer[grep("bZIP", homer$Motif), ]
homer.subset2 <- homer[grep("ETS", homer$Motif), ]
homer.subset3 <- homer[grep("bHLH", homer$Motif), ]

homer.top5<-subset(homer,homer$rownumber <6)


require(ggplot2)

ggplot(homer, aes(y=neglogP,x=-rownumber)) +geom_bar(stat="identity") +
  geom_bar(data=homer.subset3,fill="seagreen",stat="identity") +
  geom_bar(data=homer.subset2,fill="royalblue",stat="identity") +
  geom_bar(data=homer.subset1,fill="blue",stat="identity") +
  coord_flip() +theme_bw() +theme(axis.text.y=element_blank(),axis.ticks=element_blank()) +theme(aspect.ratio=1)
