#!/usr/bin/perl -w
use strict;
use Getopt::Long qw(GetOptions);
use Time::Piece;
use DBI;
use File::Basename;
use Cwd;

my $input_file;
my $output_file;
my $num_window = 5;
my $exclude_region;
#my $output_dir = cwd();

my $usage = <<__EOUSAGE__;

Usage:

$0 [options]

required options:

	-i <string> Input Bedgraph file  
	-o <string> Output Bedgraph file

optional options:    
	-w <integer> Number of smoothed windows (default: $num_window)
	-e <string>  Region to be excluded (example: chr2:1000-2000)
	
__EOUSAGE__



GetOptions (
	'i=s' => \$input_file,
	'o=s' => \$output_file,
	'e=s' => \$exclude_region,
	'w=i' => \$num_window  
);

unless ($input_file && $output_file) {
	die "Please input and output file\n$usage";
}

&main();

sub main {
	
	my $bg_in_file = $input_file;
	my $bg_ex_file = $input_file.".excluded";
	my $bg_out_file = $output_file;
	if ($exclude_region) {
		if ($exclude_region =~ /(.*):(.*)-(.*)/) {
			my $region = "$1\t$2\t$3";
			system("echo '$region' | bedtools intersect -a $bg_in_file -b stdin -v > $bg_ex_file");
			$bg_in_file = $bg_ex_file;
		} else {
			print "invalid -e option. No region will be removed";
		}		
	}	 
	open(IN_FILE, $bg_in_file) or die "Cannot open file $bg_in_file";
	open( my $out_fh, ">$bg_out_file") or die "Cannot open file $bg_out_file";
	my $pre_chr = "";
	my $bin_size = 0;
	my @smoothed_pos = ();
	my @smoothed_cnt = ();
	my $min_coord = 0;
	my $max_coord = 0;
	my %raw_cnt = ();
	while(<IN_FILE>) {
		next if (/^#/);
		chomp;
		my ($chr, $start_pos, $end_pos, $cnt) = split(/\t/);
		if ($bin_size == 0) {
			$bin_size = $end_pos - $start_pos;
			print "bin size: $bin_size\n";
		}
		if ($chr ne $pre_chr) {			
			if ($max_coord != 0) {
				&doSmooth($pre_chr, $out_fh, $min_coord, $max_coord, $bin_size, $num_window, \%raw_cnt);
				%raw_cnt = ();				
			}
			$min_coord = $start_pos;
		}		
		$pre_chr = $chr;
		$max_coord = $start_pos;
		$raw_cnt{$start_pos} = $cnt;
	}
	&doSmooth($pre_chr, $out_fh, $min_coord, $max_coord, $bin_size, $num_window, \%raw_cnt);	
	close(IN_FILE);
	#system("rm $bg_in_file $bg_out_file");
	system("chmod 775 $output_file");
	system("chgrp khanlab $output_file");
}

sub doSmooth {
	my ($chr, $out_fn, $min_coord, $max_coord, $bin_size, $num_window, $raw_cnt) = @_;
	my @win = ();
	my $win_sum = 0;
	#initial window
	for (my $i=0;$i<=$num_window;$i++) {
		push @win, 0;
	}
	
	for (my $i=$min_coord;$i<=$min_coord + ($num_window - 1) * $bin_size;$i+=$bin_size) {
		if ($$raw_cnt{$i}) {
			push @win, $$raw_cnt{$i};
			$win_sum += $$raw_cnt{$i};
		} else {
			push @win, 0;
		}
	}
	
	
	
	for (my $i=$min_coord;$i<=$max_coord;$i+=$bin_size) {
		my $removed = shift @win;
		$win_sum = $win_sum - $removed;
		if ($$raw_cnt{$i + $num_window * $bin_size}) {
			push @win, $$raw_cnt{$i + $num_window * $bin_size};
			$win_sum = $win_sum + $$raw_cnt{$i + $num_window * $bin_size};
		} else {
			push @win, 0;			
		}
		if (abs($win_sum) > 10e-10) {
			print $out_fn join("\t", $chr, $i, $i + $bin_size, $win_sum / (2*$num_window+1))."\n";
		}		
	}
	
}