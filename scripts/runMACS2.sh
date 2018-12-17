!/usr/bin/env bash

module load macs

macs2 callpeak -t $sample_file -c $control_file -f BAM --name $sample_file --outdir $outdir --extsize $extsize

