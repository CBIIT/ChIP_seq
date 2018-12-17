#!/bin/bash

cd $PBS_O_WORKDIR

module load R
module load ngsplot

ngs.plot.r -G hg19 -R bed -C $config_file -O $output_name -GO none