library(ggplot2)
library(grid)
library(gridExtra)

##if using txt file
DF <- read.table("PAX3-FOXO1_ChipSeq/Sample_RH4_PAX3FOXO1_SRR.32.35_C_Combined_peaks.narrowPeak.nobl.bed.annotation.txt", sep="\t", header=T)

##using excel format (useful for manually cleaning intronic entries in the Annotation column)
require(XLConnect)
wb = loadWorkbook("PAX3-FOXO1_ChipSeq/Sample_RH4_PAX3FOXO1_SRR.32.35_C_Combined_peaks.narrowPeak.nobl.bed.annotation.xlsx")
DF = readWorksheet(wb, sheet = "Sample_RH4_PAX3FOXO1_SRR.32.35_", header = TRUE)


DF <- subset(DF,DF$Annotation %in% c("exon","Intergenic","intron","promoter-TSS"))
DF$neglog10.qval <- DF$Peak.Score/10
#DF$group <- cut(DF$neglog10.p.value,breaks = seq(-40, 120, by = 80))
DF$log10.distTSS <- log10(abs(DF$Distance.to.TSS))



p1 <- ggplot(DF,aes(x=neglog10.qval,y=log10.distTSS,colour=factor(Annotation))) + geom_point() +
  scale_x_continuous(expand=c(0.02,0),breaks=c(0,40,80,120,160)) +
  scale_y_continuous(expand=c(0.02,0)) +
  theme_bw() +
  theme(legend.position="none",plot.margin=unit(c(0,0,0,0),"points"))

theme0 <- function(...) theme( panel.background = element_blank(), #legend.position = "none",
                               panel.grid.major = element_blank(),
                               panel.grid.minor = element_blank(),
                               panel.margin = unit(0,"null"),
                               axis.ticks = element_blank(),
                               axis.text.x = element_blank(),
                               axis.text.y = element_blank(),
                               axis.title.x = element_blank(),
                               axis.title.y = element_blank(),
                               axis.ticks.length = unit(0,"null"),
                               axis.ticks.margin = unit(0,"null"),
                               panel.border=element_rect(color=NA),...)

p2 <- ggplot(DF,aes(x=neglog10.qval,colour=factor(Annotation),fill=factor(Annotation))) + 
  geom_histogram(alpha=0.5) + 
  scale_x_continuous(breaks=NULL,expand=c(0.02,0)) +
  scale_y_continuous(breaks=NULL,expand=c(0.02,0)) +
  theme_bw() +
  theme0(plot.margin = unit(c(1,0,0,2.2),"lines")) 

p3 <- ggplot(DF,aes(x=log10.distTSS,colour=factor(Annotation),fill=factor(Annotation))) + 
  geom_histogram(alpha=0.5) + 
  coord_flip()  + 
  scale_x_continuous(labels = NULL,breaks=NULL,expand=c(0.02,0)) +
  scale_y_continuous(labels = NULL,breaks=NULL,expand=c(0.02,0)) +
  theme_bw() +
  theme0(plot.margin = unit(c(0,1,1.2,0),"lines"))

grid.arrange(arrangeGrob(p2,ncol=2,widths=c(3,3)),
             arrangeGrob(p1,p3,ncol=2,widths=c(3,3)),
             heights=c(1,3))

median(DF$Distance.to.TSS)

###extract summary counts of each class
peak.stats<-as.data.frame(table(DF$Annotation))

###reviewer wants a scale!
ggplot(DF,aes(x=log10.distTSS,colour=factor(Annotation),fill=factor(Annotation))) + 
  geom_histogram(alpha=0.5) + 
  coord_flip()  + 
  theme_bw() +
  theme0(plot.margin = unit(c(0,1,1.2,0),"lines"))

