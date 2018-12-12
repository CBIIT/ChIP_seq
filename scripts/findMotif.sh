#!/usr/bin/env bash

module load homer

#input_file=/data/khanlab/projects/ChIP_seq/data_by_file_type/test/DATA/Sample_RH4_H3K27ac_001_C_H5TLGBGXX/ROSE_out_narrow_1000/trimmed_AllEnhancers.sorted.table.txt_fpkm3.regular.bed
#output_dir=/data/khanlab/projects/ChIP_seq/data_by_file_type/test/DATA/Sample_RH4_H3K27ac_001_C_H5TLGBGXX/ROSE_out_narrow_1000/motif_regular
#preparse_dir=/data/khanlab/projects/ChIP_seq/data_by_file_type/test/scripts/preparse

findMotifsGenome.pl $input_file $genome $output_dir -size $size -len $len -preparsedDir $preparsedDir -p $p
chgrp -R khanlab $output_dir
chmod -R 775 $output_dir
