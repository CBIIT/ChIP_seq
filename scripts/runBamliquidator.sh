#!/usr/bin/env bash

module load bamliquidator

export PATH=$PATH:/usr/local/apps/bamliquidator/0.9/pipeline/bamliquidator_internal/
cd /usr/local/apps/bamliquidator/0.9/pipeline/bamliquidator_internal/
#/usr/local/apps/bamliquidator/0.9/pipeline/bamliquidator_internal/bamliquidatorbatch/bamliquidator_batch.py -n 32 -r $bed_file -f -e $extend_length -o $output_dir $input_file

python /usr/local/apps/bamliquidator/0.9/pipeline/bamToGFF_turbo.py -b $input_file -i $bed_file -e $extend_length -c $bin_size -o $input_file.bin.txt

