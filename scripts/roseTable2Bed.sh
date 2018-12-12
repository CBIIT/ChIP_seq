#!/usr/bin/env bash

if [ $# -ne 3 ]
then
  echo "Usage: `basename $0` [Input ROSE TABLE] [Output Regular BED FILE] [Ouput super BED file]"
  exit 65
fi

TEMP=$(mktemp /tmp/temporary-file.XXXXXXXX)
grep -v '^[#|REGION]' $1 | awk -v OFS='\t' -F'\t' '$NF==0 {for(i=2; i<=NF; i++) {printf $i;printf (i<NF?"\t":"\n")}}' > ${TEMP}
bedtools sort -i ${TEMP} > $2
grep -v '^[#|REGION]' $1 | awk -v OFS='\t' -F'\t' '$NF==1 {for(i=2; i<=NF; i++) {printf $i;printf (i<NF?"\t":"\n")}}' > ${TEMP}
bedtools sort -i ${TEMP} > $3

