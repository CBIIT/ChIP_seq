#!/usr/bin/env bash

module load sicer

cd $PBS_O_WORKDIR
echo ${sicer_path}SICER.sh $input_dir $input_bed $control_bed $output_dir $genome $redundancy $window $frag_size $genome_fraction $gap $qvalue_cutoff
${sicer_path}SICER.sh $input_dir $input_bed $control_bed $output_dir $genome $redundancy $window $frag_size $genome_fraction $gap $qvalue_cutoff
chgrp -R khanlab $output_dir
chmod -R 775 $output_dir


