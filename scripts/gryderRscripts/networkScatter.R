#####
#
# 140521 
# Federation
#
#####

library(ggplot2)

#args = commandArgs(TRUE)
degFile = "RMS_Epigenetics/coltron/coltron_out/coltron_test_DEGREE_TABLE.txt"
sizeTable = "RMS_Epigenetics/coltron/coltron_out/coltron_test_CANIDATE_TF_AND_SUPER_TABLE.txt"
output = "RMS_Epigenetics/coltron/coltron_out/coltron_test_NETWORK_SCATTER.pdf"

degreeData = read.delim(degFile)
sortedDegreeData = degreeData[order(degreeData$Total_Connections, decreasing=T),]
sortedDegreeData = sortedDegreeData[1:50,]

TFnames = as.character(sortedDegreeData$Tf)
TFandSuperTable = read.delim(sizeTable)
  
# Make the scatter with points based on enhancer size

m = matrix(nrow=length(sortedDegreeData$Tf),ncol=3)
rownames(m) = sortedDegreeData$Tf
colnames(m) = c('In','Out','Length')
m[,1] = sortedDegreeData$In_Degree
m[,2] = sortedDegreeData$Out_Degree
m[,3] = rep(0,length(TFnames))

for(i in seq(length(row.names(TFandSuperTable)))) {
  gene = as.character(TFandSuperTable[i,2])
  if(is.element(gene, TFnames)) {
    m[gene,3] =+ (TFandSuperTable[i,7])
  }
}

pdf(output)
qplot(m[,1],m[,2], size=m[,3],label=rownames(m)) +
  geom_text(parse=T, size=m[,3]/15000, hjust=0.5, vjust=0.2, position=position_jitter(h=1,w=1)) +
  labs(x = 'Inward Binding', y='Outward Binding') +
  theme_bw() +
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank())

dev.off()