###Make NGS configs for a list of Gene sets

genesetpath = "RMS_Epigenetics/EpiDrugs/NGSplots/HDACi_Pol2_K27ac/GeneSetNGS/Spike/"

config = read.table(paste(genesetpath,"NGSgenes.config.txt",sep=""), sep="\t", header=T)

files <- list.files(path=genesetpath, pattern=".genelist", full.names=T, recursive=FALSE)

lapply(files, function(x) {
  config$RegionsGenes <<- basename(x)
  configfilename <<- gsub(".genelist",".config.txt",basename(x))
  write.table(config,file=paste(genesetpath,configfilename,sep=""),sep="\t", row.names=F,col.names=F)
})


###Make NGS configs for a list of BED files
genesetpath = "RMS_Epigenetics/EpiDrugs/NGSplots/HDACi_Pol2_K27ac/TES/Spike/"

config = read.table(paste(genesetpath,"NGSconfig.setup.txt",sep=""), sep="\t", header=T)

files <- list.files(path=genesetpath, pattern="*.bed", full.names=T, recursive=FALSE)

lapply(files, function(x) {
  config$RegionsGenes <<- basename(x)
  configfilename <<- gsub(".bed",".config.txt",basename(x))
  write.table(config,file=paste(genesetpath,configfilename,sep=""),sep="\t", row.names=F,col.names=F)
})