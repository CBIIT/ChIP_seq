#!/usr/bin/env bash

module load bamliquidator
module load R
module load python

export PATH=$PATH:$7
#cd /usr/local/apps/bamliquidator/0.9/pipeline
cd /usr/local/apps/bamliquidator/pipeline
ulimit -c unlimited
echo "./ROSE2_main.py -i $1 -g $2 -r $3 -t $4 -s $5 -o $6"
./ROSE2_main.py -i $1 -g $2 -r $3 -t $4 -s $5 -o $6
#/usr/local/apps/bamliquidator/0.9/pipeline/python ROSE2_geneMapper.py -g HG19 -i /data/khanlab/projects/ChIP_seq/DATA/Sample_RH4_H3K27ac_001_C_H5TLGBGXX/ROSE_OUT_7fold/Sample_RH4_H3K27ac_001_C_H5TLGBGXX_SuperEnhancers.table.txt
