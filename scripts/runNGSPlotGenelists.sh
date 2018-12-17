#!/bin/bash
#
# this file is myjob.sh
#
#PBS -N ngsplot
#PBS -m be
#
cd $SLURM_SUBMIT_DIR

module load R
module load ngsplot

outputname_genebody=${output_name}.genebody
outputname_TSS=${output_name}.TSS
outputname_TES=${output_name}.TES

echo ngs.plot.r -G hg19 -R genebody -C $sample_config -O $outputname_genebody     -T $output_name     -L 5000 -P 0 -RR 5 
echo ngs.plot.r -G hg19 -R tss      -C $sample_config -O $outputname_TSS -T $output_name -L 3000 -P 0 -RR 5 #-YAS 0,150 #-SC 0,3 -D ensembl -L 3200 -GO none
echo ngs.plot.r -G hg19 -R tes      -C $sample_config -O $outputname_TES -T $output_name -L 3000 -P 0 -RR 5 #-YAS 0,80  #-SC 0,3 -D ensembl -L 3200 -GO none

ngs.plot.r -G hg19 -R genebody -C $sample_config -O $outputname_genebody     -T $output_name     -L 5000 -P 0 -RR 5 #-YAS 0,150 #-SC 0,3 -D ensembl -L 3200 -GO none 
ngs.plot.r -G hg19 -R tss      -C $sample_config -O $outputname_TSS -T $output_name -L 3000 -P 0 -RR 5 #-YAS 0,150 #-SC 0,3 -D ensembl -L 3200 -GO none
ngs.plot.r -G hg19 -R tes      -C $sample_config -O $outputname_TES -T $output_name -L 3000 -P 0 -RR 5 #-YAS 0,80  #-SC 0,3 -D ensembl -L 3200 -GO none
