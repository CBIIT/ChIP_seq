#!/usr/bin/env bash

module load bedtools
module load R
module load ngsplot
module load homer

#qsub -N BCHN -l nodes=1:gpfs,mem=250gb,ncpus=32 -v bed1=a -v config=BCHN/CTR_DvT.BCHNconfig -v s=5 -v output_dir=/data/khanlab/projects/ChIP_seq/data_by_file_type/test/scripts/output runBCHN.sh

SCRIPT_HOME=/data/khanlab/projects/ChIP_seq/data_by_file_type/test/scripts

$SCRIPT_HOME/BCHN.pl -a $bed1 -b $bed2 -g $g -c $config -e $e -n -m -t '-GO total -P 0 -IN 0 -L 3200 ' -o $output_dir #-s $s -YAS 0,0.5

bedoverlap=`$output_dir/$bed1_inter_$bed2`
echo bedoverlap

Rscript $SCRIPT_HOME/plotBCHNVenn.R $bed1 $bed2 $bedoverlap

#EDENbed=$output_dir/combined.$e.bed
#outputfilePrefix=combined.$e
#echo $SCRIPT_HOME/EDEN.pl -b $EDENbed -e $exp -d $d -o $output_dir -x $outputfilePrefix -t $t -n $n -a $a
#$SCRIPT_HOME/EDEN.pl -b $EDENbed -e $exp -d $d -o $output_dir -x $outputfilePrefix -c -t $t -n $n -a $a

chgrp khanlab -R $output_dir

#-t '-GO total –SC 0,5' -s $s
#Usage:	
#	
#$0 [options]	
#	
#required options:	
#	
#  -a	BED1
#  -b	BED2
#  -c	Config file
#  -o	Output directory
#  	
#optional options:    	
#  -e	BED extention length (default: $bed_ext_len)  
#  -s	The column to be sorted (e.g -s 5 : sort by 5th column)
#  -u	sbatch option (default: $sbatch_option)
#  -m	No NGSPlot
#  -t	NGSPlot option (default: $ngsplot_option)
# HOMER optional options:	
#  -n	Do not run HOMER
#  -g	Genome (default: $homer_genome)
#  -r	Preparsed dir (default: $homer_preparsedDir)
#  -p	Number of threads (default: $homer_p)
#  -z	Motif search size (default: $homer_size)
#  -l	Motif length (default: $homer_len) 
