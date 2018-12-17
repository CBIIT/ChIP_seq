#!/usr/bin/env bash
#PBS -m be
#
cd $PBS_O_WORKDIR

module load bedtools

#bedtools bamtobed -i /home/chouh/ChIP_seq/DATA/Sample_RH4_H3K27ac_001_C_H5TLGBGXX/Sample_RH4_H3K27ac_001_C_H5TLGBGXX.dd.bam > Sample_RH4_H3K27ac_001_C_H5TLGBGXX.dd.bam.bed
#/home/chouh/ChIP_seq/scripts/bedExtention.pl -i Sample_RH4_H3K27ac_001_C_H5TLGBGXX.dd.bam.bed -e 421 -o  Sample_RH4_H3K27ac_001_C_H5TLGBGXX.dd.bam.extended.bed
#bedtools genomecov -i Sample_RH4_H3K27ac_001_C_H5TLGBGXX.dd.bam.extended.bed -g hg19.genome -bg > Sample_RH4_H3K27ac_001_C_H5TLGBGXX.dd.bam.extended.bedGraph
bedtools bedtobam -i Sample_RH4_H3K27ac_001_C_H5TLGBGXX.dd.bam.extended.bed -g hg19.genome > Sample_RH4_H3K27ac_001_C_H5TLGBGXX.dd.extend.bam
