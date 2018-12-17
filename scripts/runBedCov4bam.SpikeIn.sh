#!/bin/bash

module load bedtools

out_file_count=${out_file}.cov.raw.bed
out_file_spikescaled=${out_file}.cov.RRPM.bed
	
bedtools multicov -bams $bam_1 $bam_2 $bam_3 $bam_4 -bed $bed_file > $out_file_count

bam_path1=`dirname $bam_1` 
bam_path2=`dirname $bam_2`
bam_path3=`dirname $bam_3`
bam_path4=`dirname $bam_4`

total_reads1=`grep -n 2 $bam_path1/SpikeIn/spike_map_summary | cut -f 3`
total_reads2=`grep -n 2 $bam_path2/SpikeIn/spike_map_summary | cut -f 3`
total_reads3=`grep -n 2 $bam_path3/SpikeIn/spike_map_summary | cut -f 3`
total_reads4=`grep -n 2 $bam_path4/SpikeIn/spike_map_summary | cut -f 3`

#echos for debugging:
echo `grep -n 2 $bam_path1/SpikeIn/spike_map_summary | cut -f 3`
echo total_reads1
echo `grep -n 2 $bam_path2/SpikeIn/spike_map_summary | cut -f 3`
echo total_reads2
echo `grep -n 2 $bam_path3/SpikeIn/spike_map_summary | cut -f 3`
echo total_reads3
echo `grep -n 2 $bam_path4/SpikeIn/spike_map_summary | cut -f 3`
echo total_reads4


awk -F "\t" -v total_reads1="$total_reads1" -v total_reads2="$total_reads2" -v total_reads3="$total_reads3" -v total_reads4="$total_reads4" 'BEGIN{OFS="\t"}{print $1,$2,$3,$4*1000000/total_reads1,$5*1000000/total_reads2,$6*1000000/total_reads3,$7*1000000/total_reads4}' $out_file_count > $out_file_spikescaled

chgrp khanlab $bam_1
chgrp khanlab $bam_2
chgrp khanlab $bam_3
chgrp khanlab $bam_4
chgrp khanlab $bed_file
chgrp khanlab $out_file_count
chgrp khanlab $out_file_spikescaled