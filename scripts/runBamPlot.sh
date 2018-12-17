#!/bin/bash
#
# this file is myjob.sh
#
#PBS -N RH4_bamplot
#PBS -m be
#
cd $PBS_O_WORKDIR

module load bamliquidator
module load R

export PATH=$PATH:/data/khanlab/projects/ChIP_seq/scripts
cd /usr/local/apps/bamliquidator/0.9/pipeline/
#cd /home/chouh/pipeline-master
input_file=/home/chouh/ChIP_seq/DATA/Sample_RH4_H3K27ac_001_C_H5TLGBGXX/Sample_RH4_H3K27ac_001_C_H5TLGBGXX.dd.bam_7X_peaks.narrowPeak.bed

#input_file=/home/chouh/ChIP_seq/DATA/Sample_RH4_H3K27ac_001_C_H5TLGBGXX/ROSE_OUT_7fold/Sample_RH4_H3K27ac_001_C_H5TLGBGXX_Enhancers_withSuper.bed
input_region='chr7:.:129989618-130015223'
bamPlot_turbo.py -g HG19 -e 0 -i $input_region -b /data/khanlab/projects/ChIP_seq/DATA/Sample_RH4_H3K27ac_001_C_H5TLGBGXX/Sample_RH4_H3K27ac_001_C_H5TLGBGXX.dd.bam -o /data/khanlab/projects/ChIP_seq/DATA/Sample_RH4_H3K27ac_001_C_H5TLGBGXX/ROSE_OUT_7fold/

