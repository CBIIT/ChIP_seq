#!/usr/bin/perl -w
use strict;
use Getopt::Long qw(GetOptions);
use Time::Piece;
use File::Basename;
use Cwd;

my $input_file;
my $output_file;
my $scale_factor;

my $usage = <<__EOUSAGE__;

Usage:

$0 [options]

required options:

	-i <string> Input TDF file  
	-o <string> Output TDF file
	-f <float>  Scale factor
	
__EOUSAGE__



GetOptions (
	'i=s' => \$input_file,
	'o=s' => \$output_file,
	'f=f' => \$scale_factor
);

unless ($input_file && $output_file && $scale_factor) {
	die "Please input, output file and scale factor\n$usage";
}

&main();

sub main {

	my $bg_in_file = $input_file.".bedgraph";	
	my $bg_out_file = $output_file.".bedgraph";
	#if (!-e $bg_in_file){
	system("igvtools tdftobedgraph $input_file $bg_in_file");
	#}
	open(IN_FILE, $bg_in_file) or die "Cannot open file $bg_in_file";
	open(OUT_FILE, ">$bg_out_file") or die "Cannot open file $bg_out_file";
	while(<IN_FILE>) {
		next if (/^#/);
		chomp;
		my ($chr, $start_pos, $end_pos, $cnt) = split(/\t/);
		$cnt = $cnt * $scale_factor;
		print OUT_FILE "$chr\t$start_pos\t$end_pos\t$cnt\n";		
	}
	close(IN_FILE);
	system("igvtools toTDF $bg_out_file $output_file hg19");
	system("rm $bg_out_file $bg_in_file");
	system("chmod 775 $output_file");
	system("chgrp khanlab $output_file");
}
