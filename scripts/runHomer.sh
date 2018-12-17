#!/usr/bin/env bash

module load homer

input_file=/data/khanlab/projects/ChIP_seq/DATA/Sample_RH4_H3K27ac_001_C_H5TLGBGXX/Sample_RH4_H3K27ac_001_C_H5TLGBGXX.dd.bam_7X_peaks.narrowPeak.bed
output_file=/data/khanlab/projects/ChIP_seq/DATA/Sample_RH4_H3K27ac_001_C_H5TLGBGXX/Sample_RH4_H3K27ac_001_C_H5TLGBGXX.dd.bam_7X_peaks.narrowPeak.bed.annotation.txt

annotatePeaks.pl $input_file hg19 > $output_file
