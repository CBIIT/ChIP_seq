#!/bin/bash

module load bedtools
module load samtools

out_file_count=${out_file}.raw.bed
	
bedtools multicov -bams $bam_1 $bam_2 $bam_3 $bam_4 -bed $bed_file > $out_file_count

total_reads1=`samtools view -bh -L ${bed_file} ${bam_1} |samtools flagstat - |head -1|sed -e 's/\s/\\t/g' |cut -f 1`
total_reads2=`samtools view -bh -L ${bed_file} ${bam_2} |samtools flagstat - |head -1|sed -e 's/\s/\\t/g' |cut -f 1`
total_reads3=`samtools view -bh -L ${bed_file} ${bam_3} |samtools flagstat - |head -1|sed -e 's/\s/\\t/g' |cut -f 1`
total_reads4=`samtools view -bh -L ${bed_file} ${bam_4} |samtools flagstat - |head -1|sed -e 's/\s/\\t/g' |cut -f 1`

awk -F "\t" -v total_reads1="$total_reads1" -v total_reads2="$total_reads2" -v total_reads3="$total_reads3" -v total_reads4="$total_reads4" 'BEGIN{OFS="\t"}{print $1,$2,$3,$4*1000000/total_reads1,$5*1000000/total_reads2,$6*1000000/total_reads3,$7*1000000/total_reads4}' $out_file_count > $out_file

chgrp khanlab $bam_1
chgrp khanlab $bam_2
chgrp khanlab $bam_3
chgrp khanlab $bam_4
chgrp khanlab $bed_file
chgrp khanlab $out_file_count
chgrp khanlab $out_file