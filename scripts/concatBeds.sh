#!/usr/bin/env bash

data_home='/data/khanlab/projects/ChIP_seq/DATA/'

sample_list_file=$1
evalue=$2
rose_dist=$3
out_file=$4

while read -r sample
do
	cat "${data_home}${sample}/MACS_Out_p_${evalue}/ROSE_out_${rose_dist}/${sample}_peaks_AllEnhancers.table.super.GREAT.bed" >> $out_file
	
done < $sample_list_file
