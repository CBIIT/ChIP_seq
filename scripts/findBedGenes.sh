#!/usr/bin/env bash

#command example
#qsub -N findBedGenes -l nodes=1:gpfs,mem=250gb,ncpus=32 -v home=/data/khanlab/projects/ChIP_seq/data_by_file_type/test -v db_file=/data/khanlab/projects/ChIP_seq/data_by_file_type/test/scripts/chipseq.db -v bed_file=/data/khanlab/projects/ChIP_seq/data_by_file_type/test/DATA/Sample_RH4_H3K27ac_001_C_H5TLGBGXX/ROSE_out_narrow_1000/trimmed_Enhancers_withSuperStretch.sorted.dd.bed -v exp_file=/data/khanlab/projects/working_DATA/Sample_RH4Seq_T_D21KAACXX/Sample_RH4Seq_T_D21KAACXX.Cufflinks_UCSC/genes.fpkm_tracking -v out_file=/data/khanlab/projects/ChIP_seq/data_by_file_type/test/DATA/Sample_RH4_H3K27ac_001_C_H5TLGBGXX/ROSE_out_narrow_1000/trimmed_Enhancers_withSuperStretch_nearest_expressed_gene.fpkm3.bed -v fpkm_cutoff=3 -v super_table=/data/khanlab/projects/ChIP_seq/data_by_file_type/test/DATA/Sample_RH4_H3K27ac_001_C_H5TLGBGXX/ROSE_out_narrow_1000/trimmed_AllEnhancers_SuperStretch.table.txt -v sample_id=Sample_RH4Seq_T_D21KAACXX findBedGenes.sh

echo "$home/scripts/findBedGenes.pl -d $db_file -b $bed_file -e $exp_file -o $out_file -t $super_table -f $fpkm_cutoff -s $sample_id"
$home/scripts/findBedGenes.pl -d $db_file -b $bed_file -e $exp_file -o $out_file -t $super_table -f $fpkm_cutoff -s $sample_id
