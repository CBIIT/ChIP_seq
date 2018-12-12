#!/usr/bin/perl -w
use strict;

my $bin_size = 25;
$bin_size = $ARGV[0] if ($ARGV[0]);
my $read_length = 75;

my $script_file = "/data/khanlab/projects/ChIP_seq/scripts/runBamliquidator.sh";
my $input_dir_prefix = "/data/khanlab/projects/ChIP_seq/DATA/";
my $output_dir_prefix = "/data/khanlab/projects/ChIP_seq/data_by_file_type/bamliquidator/bin_$bin_size/";
my $pbs_log_dir = "/data/khanlab/projects/ChIP_seq/data_by_file_type/bamliquidator/bin_$bin_size/pbs_log/";

my %samples = (
            'Sample_CTR_H3K27ac_001_C_H5TLGBGXX' => 421,
            'Sample_CTR_H3K27me3_001_C_H5TLGBGXX' => 400,
            'Sample_CTR_Input_001_C_H5TLGBGXX' => 396,
            'Sample_CTR_MYC_001_C_H5TLGBGXX' => 412,
            'Sample_CTR_MYOD1_001_C_H5TLGBGXX' => 406,
            'Sample_RH4_H3K27ac_001_C_H5TLGBGXX' => 380,
            'Sample_RH4_H3K27me3_001_C_H5TLGBGXX' => 385,
            'Sample_RH4_H3K4me2_001_C_H5TLGBGXX' => 373,
            'Sample_RH4_H3K4me3_001_C_H5TLGBGXX' => 394,
            'Sample_RH4_Input_001_C_H5TLGBGXX' => 382,
            'Sample_RH4_MYOD1_001_C_H5TLGBGXX' => 396,
            );

while(my($sample, $library_length) = each %samples) {
    my $jobname = $sample;
    if ($sample =~ /Sample_(.*)_001/) {
        $jobname = $1; 
    }
    my $output_dir = $output_dir_prefix.$sample;
    my $input_file = $input_dir_prefix.$sample."/".$sample.".dd.bam";
    system("mkdir -p $output_dir");
    system("mkdir -p $pbs_log_dir");
    my $extend_length = $library_length - $read_length;
    my $pbs_out = $pbs_log_dir.$sample.".o";
    my $pbs_err = $pbs_log_dir.$sample.".e";
    my $qsub_cmd = "qsub -N $jobname -l nodes=1:gpfs,mem=250gb,ncpus=32 -v bin_size=$bin_size -v extend_length=$extend_length -v output_dir=$output_dir -v input_file=$input_file -e $pbs_err -o $pbs_out $script_file";      
    system($qsub_cmd);
    #print "bamliquidator_batch.py -b $bin_size -e $extend_length -o $output_dir $input_file\n";
}
	
