#!/bin/bash

perl -ne 'if ($.==1) {($total_num_reads)=$_=~/(.*?)\s.*/;print $total_num_reads;}' /data/khanlab/projects/ChIP_seq/DATA/${1}/${1}.flagstat.txt
