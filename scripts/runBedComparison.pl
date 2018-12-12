#!/usr/bin/perl -w
use strict;
use Getopt::Long;

my $bed_list_file;
my $input_dir_prefix = "/data/khanlab/projects/ChIP_seq/DATA/";
my $pbs_log_dir = "/data/khanlab/projects/ChIP_seq/data_by_file_type/macs/pbs_log/";
my $output_dir_prefix = "/data/khanlab/projects/ChIP_seq/data_by_file_type/bed/comparison/";

my $usage = <<__EOUSAGE__;

Usage:

$0 <bed file list> 

__EOUSAGE__

$bed_list_file=$ARGV[0];

unless($bed_list_file) {
  die($usage);
}

my @bed_list = ();

open(FILE, $bed_list_file);
while(<FILE>){
    chomp;
    push @bed_list, $_;
}
close(FILE);
#foreach my $bed_file (@bed_list) {
    #my $filtered = "$input_dir_prefix$bed_file.7x.bed";
    #my $cmd = "awk '\$7>7' $input_dir_prefix$bed_file > $filtered";
    #system($cmd);
#}
for (my $i=0;$i<=$#bed_list-1;$i++) {
     for (my $j=$i+1;$j<=$#bed_list;$j++) {
          my $file_a = $bed_list[$i];
          my $file_b = $bed_list[$j];
          my ($file_a_short) = ($file_a =~ /.*Sample_(.*)\.dd\..*/);
          my ($file_b_short) = ($file_b =~ /.*Sample_(.*)\.dd\..*/);
          my $output = $output_dir_prefix.$file_a_short."_vs_".$file_b_short;
          my $cmd = "mkdir -p $output";
          my $intersect = $output."/$file_a_short"."_inter_"."$file_b_short.bed";
          my $bed1_only = $output."/$file_a_short.only.bed";
          my $bed2_only = $output."/$file_b_short.only.bed";
          system($cmd);
          my $jobname = "comp_bed";
          my $script_file = "runBedIntersect.sh";
          my $qsub_cmd = "qsub -N $jobname -l nodes=1:gpfs,mem=250gb,ncpus=32 -v bed_1=$input_dir_prefix$file_a -v bed_2=$input_dir_prefix$file_b -v intersect=$intersect -v bed1_only=$bed1_only -v bed2_only=$bed2_only $script_file";      
          #print $qsub_cmd."\n";
          system($qsub_cmd);
          
     }     
}

