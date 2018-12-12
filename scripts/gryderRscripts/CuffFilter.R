###load protein coding list
coding <- read.table("RNAseq_DATA/Protein_coding_GRCh37.71.txt", sep="\t", header=T)

###load cufflinks (.fpkm) output
EXP <- read.table("RNAseq_DATA/RH4Seq_T_D21KAACXX.UCSC.fpkm", sep="\t", header=T)

###remove non-coding RNA entries
EXP$coding = EXP$gene_id %in% coding$symbol
EXP.coding = subset(EXP, EXP$coding %in% c("TRUE"))

##write files
write.table(EXP.coding, file="RNAseq_DATA/RH4Seq_T_D21KAACXX.UCSC.fpkm.coding", sep="\t", row.names=FALSE, col.names=TRUE, quote=FALSE)
