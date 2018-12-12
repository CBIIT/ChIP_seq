#!/bin/bash

module load homer

#example: 
# sbatch -J corr --partition=ccr --ntasks=8 --mem=32g --time=24:00:00 --export=in_file=/data/khanlab/projects/ChIP_seq/data_by_file_type/test/scripts/bedsplitmotif_test/t.bed,config_file=/data/khanlab/projects/ChIP_seq/data_by_file_type/test/scripts/bedsplitmotif_test/bedsplit.conf /data/khanlab/projects/ChIP_seq/scripts/BEDSplitMotif.sh

/data/khanlab/projects/ChIP_seq/scripts/BEDSplitMotif.pl -i $in_file -c $config_file
