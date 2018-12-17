#!/bin/bash

module load igvtools
module load bedtools

#inputs: bed_file, bam_file, tdf_file, out_file

#make sure to use a TDF which is not RPM normalized (doesn't have RPM in the name, just *.50.tdf for example)

in_file_bg=${tdf_file}.bedgraph
out_file_count=${out_file}.raw.bed
out_bg=${out_file%.*}.bedgraph
out_file_TDF=${out_file%.*}.RPMPR.tdf

#run bedtools to count reads within all BED regions
bedtools multicov -bams $bam_file -bed $bed_file > $out_file_count
	
#sum all reads counted within regions of $out_file_count
readsum=`awk -F"\t" '{sum+=$4} END {print sum}' $out_file_count` 

#make a bedgraph if there isn't already one for the corresponding TDF 
if [ ! -f $in_file_bg ]
then
	igvtools tdftobedgraph $in_file_tdf $in_file_bg
fi

#convert input bedgraph to RPMP by dividing the 4th column by the total reads in the peaks, then multiplying by a million
if [ ! -f $out_bg ]
then
	awk -F $'\t' -v readsum=$readsum 'function abs(x){return ((x < 0.0) ? -x : x)} BEGIN{OFS=FS}{print $1,$2,$3,($4/readsum*1000000)}' $in_file_bg > $out_bg
fi
 
chmod 775 $out_bg
chgrp khanlab $out_bg
igvtools toTDF $out_bg $out_file_TDF hg19

chmod 775 $out_file_TDF
chgrp khanlab $out_file_TDF

rm $out_file_count