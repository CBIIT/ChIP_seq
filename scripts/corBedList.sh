#!/bin/bash

module load bedtools
module load R

#example: sbatch -J corr --partition=ccr --ntasks=8 --mem=32g --time=24:00:00 --export=sample_list=/data/khanlab/projects/ChIP_seq/data_by_file_type/test/scripts/corr/Mari/sample_list.txt,bed_file_list=/data/khanlab/projects/ChIP_seq/data_by_file_type/test/scripts/corr/Mari/bed_list.txt,bam_file_list=/data/khanlab/projects/ChIP_seq/data_by_file_type/test/scripts/corr/Mari/bam_list.txt,out_file=/data/khanlab/projects/ChIP_seq/data_by_file_type/test/scripts/corr/Mari/out.matrix /data/khanlab/projects/ChIP_seq/data_by_file_type/test/scripts/corBedList.sh
#command to generate bed_list.txt (may need to modify)
#perl -sne 'chomp;$bed_file="$data_home/$_/MACS_Out_p_1e-07/$_"."_peaks.narrowPeak.nobl.bed";$bed_file2="$data_home/$_/MACS_Out_p_1e-07_/$_"."_peaks.narrowPeak.nobl.bed";$bed_file3="$data_home/$_/MACS_Out_p_1e-05/$_"."_peaks.narrowPeak.nobl.bed";if (-e $bed_file) {print "$bed_file\n"} elsif (-e $bed_file2) {print "$bed_file2\n";} elsif (-e $bed_file3) {print "$bed_file3\n";} else{print("not found - $bed_file2\n")}' -- -data_home=$data_home sample_list.txt


data_home='/data/khanlab/projects/ChIP_seq/DATA'
script_home='/data/khanlab/projects/ChIP_seq/data_by_file_type/test/scripts'

if [ ! -f "$bed_file_list" ]
then
	perl -sne 'chomp;$bed_file="$data_home/$_/MACS_Out_p_1e-07/$_"."_peaks.narrowPeak.nobl.bed";$bed_file2="$data_home/$_/MACS_Out_p_1e-07_/$_"."_peaks.narrowPeak.nobl.bed";$bed_file3="$data_home/$_/MACS_Out_p_1e-05/$_"."_peaks.narrowPeak.nobl.bed";if (-e $bed_file) {print "$bed_file\n"} elsif (-e $bed_file2) {print "$bed_file2\n";} elsif (-e $bed_file3) {print "$bed_file3\n";} else{print("not found - $bed_file2\n")}' -- -data_home=$data_home $sample_list > $bed_file_list
fi

if [ ! -f "$bam_file_list" ]
then
	perl -sne 'chomp;($type, $name) = $_ =~ /Sample_(.*?)_(.*?)_.*/;$total_reads=readpipe("$script_home/getTotalReads.sh $_");print "$name\t$total_reads\t$data_home/$_/$_.bam\n"' -- -script_home=$script_home -data_home=$data_home $sample_list > $bam_file_list
fi

bed_list="$(tr '\n' ' ' < $bed_file_list)"
bam_list="$(cut -f3 $bam_file_list | tr '\n' ' ')"
echo $bam_list

if [ ! -f combined.bed ]
then
	cat $bed_list > combined.bed
fi
if [ ! -f merged.bed ]
then
	sort -k1,1 -k2,2n combined.bed | cut -f 1-3 > combined.sorted.bed
	bedtools merge -i combined.sorted.bed > merged.bed
fi
if [ ! -f merged_cov.bed ]
then
	echo bedtools multicov -bams $bam_list -bed merged.bed > merged_cov.bed
	bedtools multicov -bams $bam_list -bed merged.bed > merged_cov.bed
fi
cut -f1 $bam_file_list | perl -ne 'chomp;$s.=$_."\t";END{$s=~s/\t$//;$s;print "$s\n"}' > $out_file
cut -f4- merged_cov.bed >> $out_file
Rscript ${script_home}/plot_corr.r merged_cov.bed $bam_file_list
