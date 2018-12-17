#!/usr/bin/env bash

module load java
module load bedtools

#example:
#sbatch -J xxx --partition=ibfdr --ntasks=32 --mem=32g --export=sample_list=sample_list.txt,output_dir=/data/khanlab/projects/ChIP_seq/data_by_file_type/test/scripts/chromHMM_Joana,input_sample=Sample_CTR_D48_input_002_C_HCK7HBGXX runChromHMM.sh
#if cellmarkfiletable.txt already exists:
#sbatch -J xxx --partition=ibfdr --ntasks=32 --mem=32g --export=sample_list=sample_list.txt,output_dir=/data/khanlab/projects/ChIP_seq/data_by_file_type/test/scripts/chromHMM_Joana runChromHMM.sh

jar_file=/data/khanlab/apps/ChromHMM/ChromHMM2015.jar

home_dir=/data/khanlab/projects/ChIP_seq/

if [ ! -f "$output_dir/cellmarkfiletable.txt" ]
then
   perl -sne 'chomp;($type, $name) = $_ =~ /Sample_(.*?)_(.*?)_.*/;print "$type\t$name\t$_.dd.bed\t$input_sample.dd.bed\n"' -- -input_sample=$input_sample $sample_list > $output_dir/cellmarkfiletable.txt
fi

while read file || [[ -n "$file" ]]
do
if [ ! -f "${home_dir}/data_by_file_type/bed/bam2bed/${file}.dd.bed" ]
then
	echo "bedtools bamtobed -i ${home_dir}/DATA/${file}/${file}.dd.bam > /${home_dir}/data_by_file_type/bed/bam2bed/${file}.dd.bed"
    bedtools bamtobed -i ${home_dir}/DATA/${file}/${file}.dd.bam > /${home_dir}/data_by_file_type/bed/bam2bed/${file}.dd.bed
fi
done <${sample_list}

if [ ! -f "${home_dir}/data_by_file_type/bed/bam2bed/${input_sample}.dd.bed" ]
then
    bedtools bamtobed -i ${home_dir}/DATA/${file}/${file}.dd.bam > /${home_dir}/data_by_file_type/bed/bam2bed/${file}.dd.bed
fi

java -mx32G -jar $jar_file BinarizeBed ${home_dir}/scripts/hg19.genome ${home_dir}/data_by_file_type/bed/bam2bed/ $output_dir/cellmarkfiletable.txt $output_dir/binary
num_states=(14 15 16)
for num_state in "${num_states[@]}"
do
  java -mx32G -jar $jar_file LearnModel $output_dir/binary $output_dir/output_${num_state} ${num_state} hg19
done
