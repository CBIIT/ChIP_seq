#!/usr/bin/env bash

module load igvtools 

igvtools count -w $bin_size -e $extend_length $input_file $output_file hg19


