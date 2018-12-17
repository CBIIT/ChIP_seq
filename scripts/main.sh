#!/bin/sh
module load samtools

interval=$1
bam=$2
out_dir=$3
out_file=$4

TotalReads=`samtools view -bh -L ${interval} ${bam} |samtools flagstat - |head -1|sed -e 's/\s/\\t/g' |cut -f 1`
sh ../copyNumber.sh ${TotalReads} ${interval} ${bam} ${out_file}

