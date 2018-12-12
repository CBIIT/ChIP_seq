#!/bin/bash

module load igvtools
module load bedtools

out_file1=${in_file1%.*}.bedgraph
out_file2=${in_file2%.*}.bedgraph

if [ -z $method ]
then
	method="subtraction"
fi

if [ -z $pseudo_count ]
then
	pseudo_count=0
fi
		
if [ ! -f $out_file1 ]
then
	igvtools tdftobedgraph $in_file1 $out_file1
fi
if [ ! -f $out_file2 ]
then
	igvtools tdftobedgraph $in_file2 $out_file2
fi
out_union_file=${out_file%.*}.union.bedgraph

if [ ! -f $out_union_file ]
then
	bedtools unionbedg -i $out_file1 $out_file2 > $out_union_file
fi
out_bg=${out_file%.*}.bedgraph
if [ ! -f $out_file ]
then
	if [ "$method" == "division" ]
	then
		echo "awk -F $'\t' -v pseudo_count="$pseudo_count" 'BEGIN{OFS=FS}{print $1,$2,$3,($4+pseudo_count)/($5+pseudo_count)}' $out_union_file > $out_bg"
		awk -F $'\t' -v pseudo_count="$pseudo_count" 'BEGIN{OFS=FS}{print $1,$2,$3,($4+pseudo_count)/($5+pseudo_count)}' $out_union_file > $out_bg		
	else
		awk -F $'\t' 'function abs(x){return ((x < 0.0) ? -x : x)} BEGIN{OFS=FS}{print $1,$2,$3,($4-$5)}' $out_union_file > $out_bg
	fi
	
fi
#rm $out_file1 $out_file2 $out_union_file
chmod 775 $out_bg
chgrp khanlab $out_bg
igvtools toTDF $out_bg $out_file hg19
chmod 775 $out_file
chgrp khanlab $out_file

