#!/usr/bin/env bash

module load bamliquidator
module load R
module load python

export PATH=$PATH:$8
cd /usr/local/apps/bamliquidator/pipeline
python -V
ulimit -c unlimited
echo "./ROSE2_main.py -i $1 -g $2 -r $3 -c $4 -t $5 -s $6 -o $7"
./ROSE2_main.py -i $1 -g $2 -r $3 -c $4 -t $5 -s $6 -o $7
#/usr/local/apps/bamliquidator/0.9/pipeline/python ROSE2_geneMapper.py -g HG19 -i /data/khanlab/projects/ChIP_seq/DATA/Sample_RH4_H3K27ac_001_C_H5TLGBGXX/ROSE_OUT_7fold/Sample_RH4_H3K27ac_001_C_H5TLGBGXX_SuperEnhancers.table.txt
