#!/usr/bin/env bash

module load java
module load bedtools

jar_file=/data/khanlab/apps/ChromHMM/ChromHMM2015.jar

home_dir=/data/khanlab/projects/ChIP_seq/data_by_file_type/test

files=("Sample_CTR_D48_input_002_C_HCK7HBGXX" "Sample_CTR_H3K27ac_001_C_H5TLGBGXX" "Sample_CTR_H3K27me3_001_C_H5TLGBGXX" "Sample_CTR_D48_H3K36me3_004_C_HH7C5BGXX" "Sample_CTR_H3K4me1_002_C_HCK7HBGXX" "Sample_CTR_H3K4me2_002_C_HCK7HBGXX" "Sample_CTR_D48_H3K4me3_004_C_HH7C5BGXX")

for file in "${files[@]}"
do
if [ ! -f "/data/khanlab/projects/ChIP_seq/data_by_file_type/bed/bam2bed/mari/${file}.dd.bed" ]
then
    bedtools bamtobed -i /data/khanlab/projects/ChIP_seq/DATA/${file}/${file}.dd.bam > /data/khanlab/projects/ChIP_seq/data_by_file_type/bed/bam2bed/mari/${file}.dd.bed
fi
done

java -mx32G -jar $jar_file BinarizeBed $home_dir/scripts/hg19.genome /data/khanlab/projects/ChIP_seq/data_by_file_type/bed/bam2bed/mari $home_dir/scripts/cellmarkfiletable_mari_0907.txt /data/khanlab/projects/ChIP_seq/data_by_file_type/chromHMM/mari_0907/binary
num_states=(10 11 12 13 14 15)
for num_state in "${num_states[@]}"
do
  java -mx32G -jar $jar_file LearnModel /data/khanlab/projects/ChIP_seq/data_by_file_type/chromHMM/mari_0907/binary /data/khanlab/projects/ChIP_seq/data_by_file_type/chromHMM/mari_0907/output_${num_state} ${num_state} hg19
done
