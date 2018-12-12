#!/bin/bash

module load macs/2.1.0.20150420
module load homer
module load bedtools/2.22.0
module load igvtools

/data/khanlab/projects/ChIP_seq/data_by_file_type/test/scripts/runChipseqPipeline.pl $config_file
