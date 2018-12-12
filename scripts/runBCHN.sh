#!/usr/bin/env bash

module load bedtools
module load R
module load ngsplot
module load homer

#qsub -N BCHN -l nodes=1:gpfs,mem=250gb,ncpus=32 -v bed1=a -v config=BCHN/CTR_DvT.BCHNconfig -v s=5 -v output_dir=/data/khanlab/projects/ChIP_seq/data_by_file_type/test/scripts/output runBCHN.sh

SCRIPT_HOME=/data/khanlab/projects/ChIP_seq/data_by_file_type/test/scripts

$SCRIPT_HOME/BCHN.pl -a $bed1 -b $bed2 -g $g -c $config -s $s -n -t '-GO total –SC 0,5' -o $output_dir
