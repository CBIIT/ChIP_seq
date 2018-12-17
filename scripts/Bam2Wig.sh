#!/bin/sh

fileName=$1
N=`basename $fileName .bam`
samtools mpileup $fileName | perl -ne 'BEGIN{print "track type=wiggle_0 name=fileName description=fileName\n"};($c, $start, undef, $depth) = split; if ($c ne $lastC) { print "variableStep chrom=$c\n"; };$lastC=$c;next unless $. % 10 ==0;print "$start\t$depth\n" unless $depth<3;' >/data/khanlab/projects/ChIP_seq/data/wig/$N.wig
