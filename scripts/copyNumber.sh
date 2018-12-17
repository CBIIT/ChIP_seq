#!/usr/bin/env bash
#
# Author: Rajesh Patidar rajbtpatidar@gmail.com
# Takes a bed file and calculate sum, RPKM and log2(RPKM) for every bed entry
#

module load samtools/0.1.19

mil=1000000
constant="0.01"
bedFile=$1
bamFile=$2
outFile=$3

total_reads=`samtools view -bh -L ${bedFile} ${bamFile} |samtools flagstat - |head -1|sed -e 's/\s/\\t/g' |cut -f 1`

if [ -f "$outFile" ];then
	rm $outFile
fi

#total_reads=`samtools flagstat $bamFile |head -1 | sed 's/\s/\t/g' | cut -f1`
#total_reads=407932018

while read -r line
do
	elements=( $line )
	mpileup_cord="${elements[0]}:${elements[1]}-${elements[2]}";
	sum=`samtools depth -Q 10 -r $mpileup_cord $bamFile | awk '{sum += $3} END {if (NR > 0) print sum  / NR;}'`; 
	if [ -z "$sum" ];then
		sum="0"
		echo -e "${elements[0]}\\t${elements[1]}\\t${elements[2]}\\t$sum" >> $outFile
	else
		RPKM=`echo "($sum * $mil)/ $total_reads "|bc -l`
		echo -e "${elements[0]}\\t${elements[1]}\\t${elements[2]}\\t$RPKM" >> $outFile
	fi
done < $bedFile
