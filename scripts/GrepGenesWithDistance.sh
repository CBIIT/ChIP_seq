#!/bin/sh
#module load bedtools
#
# Align temp to afixed place rather then shifting it based on number of genes on either side.
#
#

if [ "$1" == "-h" ]; then
  echo "Usage: $0 Chr Start End"
  exit 0
fi
echo -e "$1\t$2\t$3\ttemp\tstrand" >temp
grep -P "^$1\t" /data/khanlab/ref/annovar/gene_coordinates >>temp
/usr/local/bedtools/bin/sortBed -i temp >temp1
START=$2
END=$3
echo -ne "$1\t$2\t$3\t"
grep -B5 temp temp1 |cut -f 2,3,4,5 |awk -v start=$START -v end=$END -F "\t" '{OFS="\t"};{if($4 == "+") print $3,start-$1; else if ($4 == "-")print $3, start-$2; else print $3}' |perl -pe 's/\n/\t/g'
grep -A5 temp temp1 |cut -f 2,3,4,5 |awk -v start=$START -v end=$END -F "\t" '{OFS="\t"};{if($4 == "+") print $3,$1-end; else if ($4 == "-")print $3, $2-end; }'|perl -pe 's/\n/\t/g'
#grep -b5 temp temp1 |cut -f 2,3,4,5 |awk -F "\t" '{OFS="\t"};{if($4 == "+") print $3,$1"##"; else if ($4 == "-")print $3,$2"##"; else print $3}'|perl -pe 's/\n/\t/g'
#grep -b5 temp temp1 |cut -f 2,3,4,5 |perl -pe 's/\n/\t/g'
echo -ne "\n"
#awk 'BEGIN{RS="\t"};{print $1,$2,$3,$4,$5,$6,$7,$8,$9,$10,$11}' temp
rm -rf temp temp1
