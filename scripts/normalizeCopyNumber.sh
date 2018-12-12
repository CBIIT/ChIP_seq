#!/usr/bin/env bash

module load bamliquidator
module load R
module load bedtools

out_folder=$1
map_file=$2
cp_file=$3

enhancer_bed=${map_file%.*}_enhancer.bed
input_bed=${map_file%.*}_input.bed
enhancer_normalized_bed=${map_file%.*}_enhancer_normalized.bed
input_normalized_bed=${map_file%.*}_input_normalized.bed
normalized_map=${map_file%.*}_normalized.txt
awk 'BEGIN{OFS="\t"}$1 != "REGION_ID"{print $2,$3,$4,$7,$1,$5,$6}' $map_file | grep -v CHROM | sort -k1,1 -k2,2n > $enhancer_bed
awk 'BEGIN{OFS="\t"}$1 != "REGION_ID"{print $2,$3,$4,$8,$1,$5,$6}' $map_file | grep -v CHROM | sort -k1,1 -k2,2n > $input_bed
bedtools intersect -a $enhancer_bed -b $cp_file -wao | awk -F'\t' 'BEGIN{OFS="\t"}{if (o1 != "" && ($1 != o1 || $2 != o2 || $3 != o3)) {cn = (total==0)? 1 : (total+o3-o2-total_length)/(o3-o2); normalized = o4/cn;print o1,o2,o3,normalized,cn;total=0;total_length=0;}o1=$1;o2=$2;o3=$3;o4=$4;total+=2^$11*$12;total_length+=$12}END{cn = (total==0)? 1 : (total+o3-o2-total_length)/(o3-o2); normalized = o4/cn;print o1,o2,o3,normalized,cn;}' > $enhancer_normalized_bed
bedtools intersect -a $input_bed -b $cp_file -wao | awk -F'\t' 'BEGIN{OFS="\t"}{if (o1 != "" && ($1 != o1 || $2 != o2 || $3 != o3)) {cn = (total==0)? 1 : (total+o3-o2-total_length)/(o3-o2); normalized = o4/cn;print o1,o2,o3,normalized,cn;total=0;total_length=0;}o1=$1;o2=$2;o3=$3;o4=$4;total+=2^$11*$12;total_length+=$12}END{cn = (total==0)? 1 : (total+o3-o2-total_length)/(o3-o2); normalized = o4/cn;print o1,o2,o3,normalized,cn;}' > $input_normalized_bed
#bedtools intersect -a $input_bed -b $cp_file -wao | awk -F'\t' 'BEGIN{OFS="\t"}{if (o1 != "" && ($1 != o1 || $2 != o2 || $3 != o3)) {cn = (total==0)? 1 : total/(o3-o2); normalized = o4/cn;print o1,o2,o3,normalized,cn;total=0;}o1=$1;o2=$2;o3=$3;o4=$4;total+=2^$11*$12;}END{cn = (total==0)? 1 : total/(o3-o2); normalized = o4/cn;print o1,o2,o3,normalized,cn;}' > $input_normalized_bed
head -n 1 $map_file > $normalized_map
paste <(cut -f5 $enhancer_bed) <(cut -f1-3 $enhancer_normalized_bed) <(cut -f6,7 $enhancer_bed) <(cut -f4 $enhancer_normalized_bed) <(cut -f4 $input_normalized_bed) >> $normalized_map
R --no-save $out_folder/ $normalized_map 'cn_normalized' 'input' < /usr/local/apps/bamliquidator/pipeline/ROSE2_callSuper.R


