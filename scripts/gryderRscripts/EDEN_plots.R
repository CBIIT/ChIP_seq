##########
##Load EDEN summary xlsx worked up for gene expression
require(XLConnect)
wb = loadWorkbook("RMS_Epigenetics/SuperEnhancers/SuperEvTE_EDEN.xlsx")
SEvTE_genes <- readWorksheet(wb, sheet = "SEvTEviolin", header = TRUE)
SEvTE_genes$enhancer <- as.factor(SEvTE_genes$enhancer)
library(ggplot2)
ggplot(SEvTE_genes, aes(x=enhancer,y=nearest.log2)) + geom_boxplot(width=0.5)
#######################
#EDEN output txt files#
#######################
#this EDEN output format has input with columns for MTF co-ocupancy and SE/TE status, etc.
#

###EDEN output with 1st parameter set###

EDEN1 <- read.table("RMS_Epigenetics/EDEN/RH4_SEvTE/RH4_SEvTE.TAD.codegene_fpkm0_nearest-genes.txt", sep="\t", header=T)
    EDEN1[,7:11] <- lapply(EDEN1[,7:11],as.factor)
    EDEN1 <- subset(EDEN1,!(EDEN1$FPKM %in% c(".")))
    EDEN1$FPKM <- as.numeric(as.character(EDEN1$FPKM))+0.0001 #gets rid of infinite values of log
    EDEN1$log2FPKM <- log(EDEN1$FPKM,2)
    EDEN1$eMTFtype <- as.factor(do.call(paste, c(EDEN1[c("enhancer","SUM")], sep = " ")))
    EDEN1$eP3Ftype <- as.factor(do.call(paste, c(EDEN1[c("enhancer","SEwP3F")], sep = " ")))
    EDEN1$eP3Ftype <- factor(EDEN1$eP3Ftype, levels=c("typical 0", "typical 1", "super 0", "super 1")) ##just setting order for plotting
    EDEN1$enhancer <- factor(EDEN1$enhancer, levels=c("typical", "super")) ##just setting order for plotting

###T-testing
EDEN1.pval.SEvTE <- signif(t.test(FPKM ~ enhancer,EDEN1)$p.value,3)
EDEN1.pval.EwP3F <- signif(t.test(FPKM ~ SEwP3F,EDEN1)$p.value,3)
EDEN1.pval.EwMTFSUM <- signif(t.test(EDEN1$FPKM[EDEN1$SUM==1], EDEN1$FPKM[EDEN1$SUM==4])$p.value,3)
EDEN1.pval.TE.wowP3F <- signif(t.test(EDEN1$FPKM[EDEN1$eP3Ftype=="typical 0"], EDEN1$FPKM[EDEN1$eP3Ftype=="typical 1"])$p.value,3)
EDEN1.pval.SE.wowP3F <- signif(t.test(EDEN1$FPKM[EDEN1$eP3Ftype=="super 0"], EDEN1$FPKM[EDEN1$eP3Ftype=="super 1"])$p.value,3)

library(ggplot2)
library(gridExtra)
EDEN1.plot.ylim=c(-2,62)
EDEN1.plotSEvTE <- ggplot(EDEN1, aes(x=enhancer,y=FPKM)) + theme_bw() + stat_boxplot(geom ='errorbar',stat_params = list(width = 0.2)) + geom_boxplot(outlier.shape = NA, width=0.7) + coord_cartesian(ylim=EDEN1.plot.ylim) + annotate("text",x=1.5,y=58,label=EDEN1.pval.SEvTE) +annotate("segment", x = 1, xend = 2, y = 55, yend = 55)
EDEN1.plotEwP3F <- ggplot(EDEN1, aes(x=SEwP3F,y=FPKM)) + theme_bw() + stat_boxplot(geom ='errorbar',stat_params = list(width = 0.2)) + geom_boxplot(outlier.shape = NA,width=0.7) + coord_cartesian(ylim=EDEN1.plot.ylim) + annotate("text",x=1.5,y=58,label=EDEN1.pval.EwP3F) + annotate("segment", x = 1, xend = 2, y = 55, yend = 55) 
EDEN1.plotEwMTFSUM <- ggplot(EDEN1, aes(x=SUM,y=FPKM)) + theme_bw() + stat_boxplot(geom ='errorbar',stat_params = list(width = 0.2)) + geom_boxplot(outlier.shape = NA,width=0.7) + coord_cartesian(ylim=EDEN1.plot.ylim) 
EDEN1.ploteP3Ftype <- ggplot(EDEN1, aes(x=eP3Ftype,y=FPKM)) + theme_bw() + stat_boxplot(geom ='errorbar',stat_params = list(width = 0.2)) + geom_boxplot(outlier.shape = NA,width=0.7) + coord_cartesian(ylim=EDEN1.plot.ylim) 

grid.arrange(EDEN1.plotSEvTE, EDEN1.ploteP3Ftype, EDEN1.plotEwP3F, EDEN1.plotEwMTFSUM, widths=c(0.7, 1, 0.7, 1), ncol=2)

ggplot(EDEN1, aes(x=eMTFtype,y=FPKM)) + theme_bw() + stat_boxplot(geom ='errorbar',stat_params = list(width = 0.2)) + geom_boxplot(outlier.shape = NA, width=0.7) + coord_cartesian(ylim=EDEN1.plot.ylim) 


###EDEN2 output with 2nd parameter###
EDEN2 <- read.table("RMS_Epigenetics/EDEN/RH4_SEvTE/RH4_SEvTE.noTAD.allgene_fpkm0_nearest-genes.txt", sep="\t", header=T)
EDEN2[,7:11] <- lapply(EDEN2[,7:11],as.factor)
EDEN2 <- subset(EDEN2,!(EDEN2$FPKM %in% c(".")))
EDEN2$FPKM <- as.numeric(as.character(EDEN2$FPKM))+0.0001
EDEN2$log2FPKM <- log(EDEN2$FPKM,2)
EDEN2$eMTFtype <- as.factor(do.call(paste, c(EDEN2[c("enhancer","SUM")], sep = " ")))
EDEN2$eP3Ftype <- as.factor(do.call(paste, c(EDEN2[c("enhancer","SEwP3F")], sep = " ")))
EDEN2$eP3Ftype <- factor(EDEN2$eP3Ftype, levels=c("typical 0", "typical 1", "super 0", "super 1")) ##just setting order for plotting
EDEN2$enhancer <- factor(EDEN2$enhancer, levels=c("typical", "super")) ##just setting order for plotting

###T-testing
EDEN2.pval.SEvTE <- signif(t.test(FPKM ~ enhancer,EDEN2)$p.value,3)
EDEN2.pval.EwP3F <- signif(t.test(FPKM ~ SEwP3F,EDEN2)$p.value,3)
EDEN2.pval.EwMTFSUM <- signif(t.test(EDEN2$FPKM[EDEN2$SUM==3], EDEN2$FPKM[EDEN2$SUM<3])$p.value,3)
EDEN2.pval.TE.wowP3F <- signif(t.test(EDEN2$FPKM[EDEN2$eP3Ftype=="typical 0"], EDEN2$FPKM[EDEN2$eP3Ftype=="typical 1"])$p.value,3)
EDEN2.pval.SE.wowP3F <- signif(t.test(EDEN2$FPKM[EDEN2$eP3Ftype=="super 0"], EDEN2$FPKM[EDEN2$eP3Ftype=="super 1"])$p.value,3)

library(ggplot2)
library(gridExtra)
EDEN2.plot.ylim=c(-2,62)
EDEN2.plotSEvTE <- ggplot(EDEN2, aes(x=enhancer,y=FPKM)) + theme_bw() + stat_boxplot(geom ='errorbar',stat_params = list(width = 0.2)) + geom_boxplot(outlier.shape = NA, width=0.7) + coord_cartesian(ylim=EDEN2.plot.ylim) + annotate("text",x=1.5,y=58,label=EDEN2.pval.SEvTE) +annotate("segment", x = 1, xend = 2, y = 55, yend = 55)
EDEN2.plotEwP3F <- ggplot(EDEN2, aes(x=SEwP3F,y=FPKM)) + theme_bw() + stat_boxplot(geom ='errorbar',stat_params = list(width = 0.2)) + geom_boxplot(outlier.shape = NA,width=0.7) + coord_cartesian(ylim=EDEN2.plot.ylim) + annotate("text",x=1.5,y=58,label=EDEN2.pval.EwP3F) + annotate("segment", x = 1, xend = 2, y = 55, yend = 55) 
EDEN2.plotEwMTFSUM <- ggplot(EDEN2, aes(x=SUM,y=FPKM)) + theme_bw() + stat_boxplot(geom ='errorbar',stat_params = list(width = 0.2)) + geom_boxplot(outlier.shape = NA,width=0.7) + coord_cartesian(ylim=EDEN2.plot.ylim) 
EDEN2.ploteP3Ftype <- ggplot(EDEN2, aes(x=eP3Ftype,y=FPKM)) + theme_bw() + stat_boxplot(geom ='errorbar',stat_params = list(width = 0.2)) + geom_boxplot(outlier.shape = NA,width=0.7) + coord_cartesian(ylim=EDEN2.plot.ylim) 

grid.arrange(EDEN2.plotSEvTE, EDEN2.ploteP3Ftype, EDEN2.plotEwP3F, EDEN2.plotEwMTFSUM, widths=c(0.7, 1, 0.7, 1), ncol=2)

ggplot(EDEN2, aes(x=eMTFtype,y=FPKM)) + theme_bw() + stat_boxplot(geom ='errorbar',stat_params = list(width = 0.2)) + geom_boxplot(outlier.shape = NA, width=0.7) + coord_cartesian(ylim=EDEN2.plot.ylim) 



######################################
###EDEN parameter comparison plots!###
######################################

grid.arrange(EDEN1.plotSEvTE, EDEN1.plotEwP3F, EDEN2.plotSEvTE, EDEN2.plotEwP3F, widths=c(1, 1, 1, 1), ncol=2)
