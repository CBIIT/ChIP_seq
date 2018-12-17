#############
#  Venn plotting for runBCHN.sh
#  Berkley Gryder
#  2018
#############

##load data. 

args = commandArgs(TRUE)
bed1.file = args[1]
bed2.file = args[2]
bedoverlap.file = args[3]
bed1.name = args[4]
bed2.name = args[5]
output = paste(bedoverlap.file,".VENN.pdf",sep="")

bed1 <- read.table(bed1.file, sep="\t", header=F)  #example: "AML/PeakCompare/AE_RAD21/BCHN/KG1a_D4_RAD21_p-11.bed"
bed2 <- read.table(bed2.file, sep="\t", header=F) #example: "AML/PeakCompare/AE_RAD21/BCHN/KG1a_A4_RAD21_p-11.bed"
bed1overlap2 <- read.table(bedoverlap.file, sep="\t", header=F) #example: "AML/PeakCompare/AE_RAD21/BCHN/DvA/KG1a_D4_RAD21_p-11_inter_KG1a_A4_RAD21_p-11.800.bed"

set1 <- nrow(bed1)
set2 <- nrow(bed2)
overlap <- nrow(bed1overlap2)

#############################################################################################
### new way to Venn!!! http://matticklab.com/index.php?title=Weighted_Venn_diagrams_in_R  ###
### install.packages(c("graph", "RBGL"), dependencies=TRUE)                               ###
### install.packages("Vennerable", repos="http://R-Forge.R-project.org")                  ###
#############################################################################################
library(Vennerable)

bedVenn = Venn(SetNames = c(basename(bed1.name),basename(bed2.name)), Weight = c(0,set1-overlap,set2-overlap,overlap))
#plot(bedVenn, show = list(Faces = FALSE))

w <- compute.Venn(bedVenn, type = "circles")
gp <- VennThemes(w)
gp[["Face"]][["11"]]$fill <-  "#CCCC99"
gp[["Face"]][["01"]]$fill <-  "#BFA157"
gp[["Face"]][["10"]]$fill <-  "#82C6FF"

gp[["SetText"]][["Set2"]]$col <-  "goldenrod3"
gp[["SetText"]][["Set1"]]$col <-  "steelblue3"

gp[["Set"]][["Set2"]]$lty <-  0
gp[["Set"]][["Set1"]]$lty <-  0


library(grid)

pdf(output)
  plot(w, gp = gp)
  dev.off()
