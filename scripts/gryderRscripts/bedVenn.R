##load data.  each has 2 bam files coverage
bed1 <- read.table("NB_Thiele/ENCODE/KCNR_C_enhancers.merged.bed", sep="\t", header=F)
bed2 <- read.table("NB_Thiele/ENCODE/ENCODE.enhancers.merged.bed", sep="\t", header=F)
bed1overlap2 <- read.table("NB_Thiele/ENCODE/KCNR_and_ENCODE.bed", sep="\t", header=F)

#math sanity check with "no overlap" file (bedtools intersect -v)
bed1no2 <- read.table("NB_Thiele/ENCODE/KCNR_notENCODE.bed", sep="\t", header=F)
bedmathcheck <- nrow(bed1no2) + nrow(bed1overlap2)
bedmathcheck2 <- nrow(bed1) - nrow(bed1overlap2)

require(venneuler)

set1 <- nrow(bed1)
set2 <- nrow(bed2)
overlap <- nrow(bed1overlap2)

v <- venneuler(c(A=set1-overlap, B=set2-overlap, 'A&B'=overlap))

v$labels<- c(
  paste("KCNR\n",set1),
  paste("ENCODE\n",set2)
)

plot(v)
