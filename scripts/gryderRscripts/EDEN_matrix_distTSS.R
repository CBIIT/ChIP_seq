##Load EDEN summary xlsx worked up for gene expression
require(XLConnect)
wb = loadWorkbook("RMS_Epigenetics/EDEN/RH4_P3F_shscr/RH4_NuRD_sh_near.xlsx")
EDEN.matrix <- readWorksheet(wb, sheet = "Matrix", header = TRUE)
EDEN.matrix <- subset(EDEN.matrix,EDEN.matrix$Class %in% c("UP","DOWN","NOCHANGE"))
EDEN.matrix$DistanceTSS = as.numeric(EDEN.matrix$DistanceTSS)
library(ggplot2)
ggplot(EDEN.matrix, aes(x=abs(DistanceTSS),colour=factor(column_5),fill=factor(column_5))) +
  geom_histogram(alpha=0.5) +theme_bw() +coord_flip() +facet_wrap(~Class)

ggplot(EDEN.matrix, aes(y=log10(abs(DistanceTSS)),x=Class,fill=Class)) +
  geom_violin(alpha=0.5) +theme_bw()  +facet_wrap(~column_5)