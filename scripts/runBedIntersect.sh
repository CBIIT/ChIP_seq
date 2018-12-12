#!/usr/bin/env bash

module load bedtools

export PATH=$PATH:/data/khanlab/projects/ChIP_seq/scripts
bedtools intersect -a $bed_1 -b $bed_2 > $intersect
bedtools intersect -a $bed_1 -b $bed_2 -v > $bed1_only
bedtools intersect -a $bed_2 -b $bed_1 -v > $bed2_only


