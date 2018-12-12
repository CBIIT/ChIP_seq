#!/bin/bash

module load macs
module load samtools
module load homer
module load bedtools
module load IGVTools

/data/khanlab/projects/ChIP_seq/scripts/runChipseqPipeline.pl $config_file
