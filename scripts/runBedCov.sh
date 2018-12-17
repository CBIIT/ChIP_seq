#!/bin/bash

module load bedtools
module load samtools

##sort the first bam file, if needed
	bam_1_sorted=${bam_1%.*}.sorted
		if [ ! -f ${bam_1_sorted}.bam ]
		then
		samtools sort $bam_1 $bam_1_sorted
		samtools index ${bam_1_sorted}.bam
		fi
	bam_1_sorted=${bam_1_sorted}.bam

##sort the second bam file, if needed
	bam_2_sorted=${bam_2%.*}.sorted
		if [ ! -f ${bam_2_sorted}.bam ]
		then
		samtools sort $bam_2 $bam_2_sorted
		samtools index ${bam_2_sorted}.bam
		fi
	bam_2_sorted=${bam_2_sorted}.bam

bedtools multicov -bams $bam_1_sorted $bam_2_sorted -bed $bed_file > $out_file
chgrp khanlab $bam_1_sorted
chgrp khanlab $bam_2_sorted
chgrp khanlab $bed_file
chgrp khanlab $out_file
