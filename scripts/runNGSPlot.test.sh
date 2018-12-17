#!/bin/bash

cd $PBS_O_WORKDIR

module load R
module load ngsplot

ngs.plot.r -G hg19 -R bed -C /data/khanlab/projects/ChIP_seq/data_by_file_type/test/scripts/BCHN/test.conf -O test
