###Circle plot

data = read.table("RMS_Epigenetics/ROSE/RMS_Tumors/ChosenPbed/RH4_SEoverlap_summaryTable.txt",sep="\t",header=T)
library(ggplot2)

ggplot(data, aes(x=reorder(Group, -Count),y=Count,fill=Group)) + geom_bar(stat = "identity") + coord_polar(theta = "y")

