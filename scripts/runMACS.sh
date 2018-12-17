#!/bin/bash
#This file is runMacs
#
#

# load the latest (default) version of MACS
#module load macs  
#this line didn't work for some reason, so we had to put the full location in the callpeak line
#cd /data/khanlab/projects/ChIP_seq/data/bam
/usr/local/apps/macs/2.0.10-20130731/bin/macs2 callpeak -t /data/khanlab/projects/ChIP_seq/data/bam/E-MTAB-1565_ERK/ERR248781/ERR248781_sort_dd.bam -c /data/khanlab/projects/ChIP_seq/data/bam/E-MTAB-1565_ERK/ERR248779/ERR248779_sort_dd.bam -f BAM -q 0.01 --name ERK2_ESC_ERR248781_sort_dd
#the function --outdir should write the files to the specified directory, but I get an error message: ...macs2: error: unrecognized arguments: --outdir /data/khanlab/projects/ChIP_seq/data/macs
