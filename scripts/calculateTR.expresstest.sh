#!/bin/bash

module load bedtools

#example: sbatch --partition=ccr --time=24:00:00 --export=sample_list=/data/khanlab/projects/ChIP_seq/data_by_file_type/test/scripts/corr/Mari/sample_list.txt /data/khanlab/projects/ChIP_seq/data_by_file_type/test/scripts/calculateTR.sh

data_home='/data/khanlab/projects/ChIP_seq/DATA'
script_home='/data/khanlab/projects/ChIP_seq/data_by_file_type/test/scripts'
#sample_list='/data/khanlab/projects/ChIP_seq/data_by_file_type/test/scripts/corr/Mari/sample_list.txt'

IFS=$'\n' read -d '' -r -a samples < $sample_list

for sample in "${samples[@]}"
do
	echo "$script_home/calculateTR.pl -b $data_home/$sample/$sample.dd.bam -e CTR_Tram_48h_T_H3YCHBGXX.gene.TPM.txt -o $sample.trv"
	$script_home/calculateTR.pl -b $data_home/$sample/$sample.dd.bam -e CTR_Tram_48h_T_H3YCHBGXX.gene.TPM.txt -o $sample.trv
done	
