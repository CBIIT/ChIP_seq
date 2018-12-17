#!/usr/bin/env bash

#command example
#grep -v -P '^[#|REGION]' trimmed_AllEnhancers_SuperStretch.table.txt | sort -k2,2 -k3n > sample_AllEnhancers.sorted.table.txt
#qsub -N findEnhancerGenes -l nodes=1:gpfs,mem=250gb,ncpus=32 -v home=/data/khanlab/projects/ChIP_seq/data_by_file_type/test -v db_file=/data/khanlab/projects/ChIP_seq/data_by_file_type/test/scripts/chipseq.db -v exp_file=/data/khanlab/projects/working_DATA/Sample_RH4Seq_T_D21KAACXX/Sample_RH4Seq_T_D21KAACXX.Cufflinks_UCSC/genes.fpkm_tracking -v fpkm_cutoff=3 -v enhancer_table=/data/khanlab/projects/ChIP_seq/data_by_file_type/test/DATA/Sample_RH4_H3K27ac_001_C_H5TLGBGXX/ROSE_out_narrow_1000/trimmed_AllEnhancers_SuperStretch.sorted.table.txt -v sample_id=Sample_RH4Seq_T_D21KAACXX findEnhancerGenes.sh
#/data/khanlab/projects/ChIP_seq/data_by_file_type/test/scripts/findEnhancerGenes.pl -d /data/khanlab/projects/ChIP_seq/data_by_file_type/genome_bins/TADs/IMR90_domains_hg19.bed -e /data/khanlab/projects/working_DATA/Sample_RH4Seq_T_D21KAACXX/Sample_RH4Seq_T_D21KAACXX.Cufflinks_UCSC/genes.fpkm_tracking -t /data/khanlab/projects/ChIP_seq/data_by_file_type/test/DATA/Sample_RH4_H3K27ac_001_C_H5TLGBGXX/ROSE_out_narrow_1000/trimmed_AllEnhancers.sorted.table.txt -f 2 -s 200000

echo "$home/scripts/findEnhancerGenes.pl -d $tad_file -e $exp_file -t $enhancer_table -f $fpkm_cutoff -s $dist_cutoff"
$home/scripts/findEnhancerGenes.pl -d $tad_file -e $exp_file -t $enhancer_table -f $fpkm_cutoff -s $dist_cutoff
