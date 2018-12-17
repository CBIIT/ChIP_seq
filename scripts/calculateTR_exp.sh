#!/bin/bash

module load bedtools

#example: sbatch --partition=ccr --time=24:00:00 --export=sample_list=/data/khanlab/projects/ChIP_seq/data_by_file_type/test/scripts/corr/Mari/sample_list.txt /data/khanlab/projects/ChIP_seq/data_by_file_type/test/scripts/calculateTR.sh

data_home='/data/khanlab/projects/ChIP_seq/DATA'
script_home='/data/khanlab/projects/ChIP_seq/data_by_file_type/test/scripts'
#sample_list='/data/khanlab/projects/ChIP_seq/data_by_file_type/test/scripts/corr/Mari/sample_list.txt'

while IFS=$'\t' read -r -a cols
do
	echo "$script_home/calculateTR.pl -b $data_home/${cols[0]}/${cols[0]}.dd.bam -o ${cols[0]}.trv -e ${cols[1]}"
        $script_home/calculateTR.pl -b $data_home/${cols[0]}/${cols[0]}.dd.bam -o ${cols[0]}.trv -e ${cols[1]}
done < $sample_list

