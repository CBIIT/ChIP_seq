#!/usr/bin/env bash

module load java
module load bedtools

jar_file=/data/khanlab/apps/ChromHMM/ChromHMM2015.jar

home_dir=/data/khanlab/projects/ChIP_seq/data_by_file_type/test

#bedtools bamtobed -i /data/khanlab/projects/ChIP_seq/DATA/Sample_RH4_30X_D_C6FEFANXX/Sample_RH4_30X_D_C6FEFANXX.dd.bam > /data/khanlab/projects/ChIP_seq/data_by_file_type/bed/bam2bed/Sample_RH4_30X_D_C6FEFANXX.dd.bed
#bedtools bamtobed -i /data/khanlab/projects/ChIP_seq/DATA/Sample_RH4_CTCF_008_C_HHC7JBGXX/Sample_RH4_CTCF_008_C_HHC7JBGXX.dd.bam > /data/khanlab/projects/ChIP_seq/data_by_file_type/bed/bam2bed/Sample_RH4_CTCF_008_C_HHC7JBGXX.dd.bed
#bedtools bamtobed -i /data/khanlab/projects/ChIP_seq/DATA/Sample_RH4_RAD21_008_C_HHC7JBGXX/Sample_RH4_RAD21_008_C_HHC7JBGXX.dd.bam > /data/khanlab/projects/ChIP_seq/data_by_file_type/bed/bam2bed/Sample_RH4_RAD21_008_C_HHC7JBGXX.dd.bed
#bedtools bamtobed -i /data/khanlab/projects/ChIP_seq/DATA/Sample_RH4_t6_input_002_C_HCK7HBGXX/Sample_RH4_t6_input_002_C_HCK7HBGXX.dd.bam > /data/khanlab/projects/ChIP_seq/data_by_file_type/bed/bam2bed/Sample_RH4_t6_input_002_C_HCK7HBGXX.dd.bed

#java -mx32G -jar $jar_file BinarizeBed $home_dir/scripts/hg19.genome /data/khanlab/projects/ChIP_seq/data_by_file_type/bed/bam2bed /data/khanlab/projects/ChIP_seq/data_by_file_type/chromHMM/berkley/run3/cellmarkfiletable.txt /data/khanlab/projects/ChIP_seq/data_by_file_type/chromHMM/berkley/run3/binary

num_states=(25)
for num_state in "${num_states[@]}"
do
	java -mx32G -jar $jar_file LearnModel /data/khanlab/projects/ChIP_seq/data_by_file_type/chromHMM/berkley/run3/binary /data/khanlab/projects/ChIP_seq/data_by_file_type/chromHMM/berkley/run3/output_$num_state $num_state hg19
done
