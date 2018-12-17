#!/bin/bash
#
# this file is myjob.sh
#
#PBS -N ngsplot
#PBS -m be
#
cd $SLURM_SUBMIT_DIR

###code to use config files to make plots
module load R
module load ngsplot

ngs.plot.r -G hg19 -R bed -C $left_config -O $left_output_name -L 1600 -P 0 -IN 0 -GO total -WD 7 -RR 10  #use bed with 400 bp flank around CTCF motif center

ngs.plot.r -G hg19 -R bed -C $right_config -O $right_output_name -L 1600 -P 0 -IN 0 -GO total -WD 7 -RR 10

ngs.plot.r -G hg19 -R bed -C $super_config -O $super_output_name -L 3000 -P 0 -IN 1 -GO total -WD 9 -RR 10 -MW 2
