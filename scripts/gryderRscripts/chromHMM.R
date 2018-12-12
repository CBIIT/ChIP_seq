#######################################################################
###make 4 column BED (with first row removed) from original HMM file###
#######################################################################

outdir = "RMS_Epigenetics/RH4/RH4chromHMM/cleanbams15/output_16/"
cellstates = "RH4_16"

###using command line: tail -n +2 RH4_16_dense.bed > RH4_16_dense.nh.bed ###
chromHMM.orig <- read.table(paste(outdir,cellstates,"_dense.nh.bed",sep=""),sep="\t", header=F)
colnames(chromHMM.orig)<-c("chrom","Start","End","emission","zero","dot","start2","end2",'RGB')
###define new state numbers
##lookup table
emission=c(4,5,3,2,1,12,15,13,14,16,10,8,7,6,11,9) #original states, ordered by new assigments
element=c('Active Promoter','Poised Promoter','Genic Enhancer','Enhancer','Weak Enhancer','Bivalent Enhancer','Strong Transcription','Transcription Transition','Weak Transcription','Insulator','ZNFgenes&repeats','Strong PC Repressed','PC Repressed','Weak PC Repressed','Heterochromatin','LowSignal') #new states
RGBnew=c('238,52,43','160,75,157','255,157,0','255,157,44','200,157,144','179,226,31','249,216,44','0,149,51','109,199,57','41,199,170','89,89,255','111,146,200','71,126,181','111,164,216','147,149,152','188,190,192') #new colors
swap.states=data.frame(emission,element,RGBnew)
library(plyr)
chromHMM.new <- join(chromHMM.orig, swap.states, by = "emission")
chromHMM <- chromHMM.new[,c(1,2,3,10)]

##write bed files
chromHMM.IGV.bed<-chromHMM.new[,c(1,2,3,10,5,6,7,8,11)]
write.table(chromHMM.IGV.bed, file=paste(outdir,cellstates,".IGV.bed",sep=""), sep="\t", row.names=FALSE, col.names=FALSE, quote=FALSE)
write.table(chromHMM, file=paste(outdir,cellstates,".bed",sep=""), sep="\t", row.names=FALSE, col.names=FALSE, quote=FALSE)
######

##if available, instead load chromHMM from prepared file
chromHMM <- read.table(paste(outdir,cellstates,".bed",sep=""),sep="\t", header=F)
colnames(chromHMM)<-c("chrom","Start","End","element")

#load chromHMM.inter dataframe (a subset result from BEDtools bedintersect)
chromHMM.inter <- read.table(paste(outdir,cellstates,"_w_RH4_EZH2_p-7.bed",sep=""),sep="\t", header=F)
colnames(chromHMM.inter)<-c("chrom","Start","End","element")

#plot it
library(ggplot2)
library(gridExtra)
p1<-ggplot(chromHMM, aes(factor(element))) + geom_bar() + coord_flip() + theme_bw()
p2<-ggplot(chromHMM.inter, aes(factor(element))) + geom_bar() + coord_flip() + theme_bw()
grid.arrange(arrangeGrob(p1,p2,nrow=1))

###initiate required dataframes
allstates = chromHMM[!duplicated(chromHMM[,c('element')]),]  #make buffer dataframe with 1 of each state to buffer samples lacking a state
files <- list.files(path=outdir, pattern="*_w_*", full.names=T, recursive=FALSE)
files <- files[grep(".bed", files)]

lapply(files, function(x) {
  
  #load chromHMM.inter dataframe (a subset result from BEDtools bedintersect)
  chromHMM.inter <- read.table(x,sep="\t", header=F)
  colnames(chromHMM.inter)<-c("chrom","Start","End","element")
  chromHMM.inter = rbind(chromHMM.inter,allstates) #add 1 of each state, so that elements tables still match if a state is missing
  interPeaknumber = nrow(chromHMM.inter)
  
  ###build counted elements dataframe from tables
  elements.chromHMM<-table(chromHMM$element)
  elements.inter<-table(chromHMM.inter$element)
   
  elements.df<-as.data.frame(elements.chromHMM)
  elements.df.inter<-as.data.frame(elements.inter)
  elements.df.inter$Freq = elements.df.inter$Freq-1 #remove 1 each from the added allstates
  elements.df$inter<-elements.df.inter[,2]
  elements.df$perc<-elements.df$inter/elements.df$Freq*100 ##calc percentages
  
  ##order the elements, from alphabetical to chromatical
  #elements.df<- within(elements.df,Var1<-factor(Var1,levels=c("PromActive","PromWeak","PromPoised","EnhStrong","EnhWeak","EnhPoised","Low","ReprPC","ReprPCWeak","Hetero")))
  elements.df$neworder<-c(1,6,4,3,15,10,16,13,2,12,7,8,5,14,9,11)
  library(plyr)
  elements.df<-ddply(elements.df, c('neworder'))
  
  ###load genomics overlay stats
  chromHMM.stats <- read.table(paste(outdir,cellstates,"_overlap.txt",sep=""),sep="\t", header=T)
  ###add new states and elements and clean off bottom row
  chromHMM.stats$emission=chromHMM.stats$state..Emission.order.
  chromHMM.stats <- join(chromHMM.stats, swap.states, by = "emission")
  chromHMM.stats=chromHMM.stats[1:nrow(elements.df),]
  
  #plot it
  #ggplot(chromHMM.stats, aes(x = factor(element), y = Genome..)) + geom_bar(stat = "identity") + coord_flip() + theme_bw()
  
  ##add columns (ordered by 'neworder' so they match up) for peaks per bp
  colnames(elements.df)<-c("element","Freq","inter","perc","neworder")
  elements.df <- join(elements.df,chromHMM.stats, by = "element")
  elements.df$chromHMM.Gb<-(elements.df$Genome..)*3129054144/1000000000 #Gb of element
  elements.df$inter.per.chromHMM.Gb<-elements.df$inter/elements.df$chromHMM.Gb ##calc percentages
  
  ##order the elements accoring to new order
  elements.df$element <- factor(elements.df$element, levels = elements.df$element[order(elements.df$neworder, decreasing = TRUE)])
  
  ###plot it
  output = paste(unlist(strsplit(x, split='bed', fixed=T)),"perGbstate.pdf",sep = '')
  output2 = paste(unlist(strsplit(x, split='bed', fixed=T)),"pdf",sep = '')
  bedname = gsub(".bed","",basename(x))
  library(ggplot2)
    interplot=ggplot(elements.df, aes(x = factor(element), y = inter.per.chromHMM.Gb)) + geom_bar(stat = "identity") + coord_flip() + theme_bw() +labs(title=paste(bedname,", n = ",interPeaknumber,sep="")) 
    ggsave(output,plot = interplot,width = 5, height = 7)
    interplot2=ggplot(elements.df, aes(x = factor(element), y = inter)) + geom_bar(stat = "identity") + coord_flip() + theme_bw() +labs(title=paste(bedname,", n = ",interPeaknumber,sep="")) 
    ggsave(output2,plot = interplot2,width = 5, height = 7)
   
})

#testing plot
ggplot(elements.df, aes(x = factor(element), y = inter.per.chromHMM.Gb)) + geom_bar(stat = "identity") + coord_flip() + theme_bw()

