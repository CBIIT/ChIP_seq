#!/bin/bash

module load bedtools

sort -k1,1 -k2,2n $bed_file | bedtools coverage -a stdin -b $bam_1 -counts > $out_file


