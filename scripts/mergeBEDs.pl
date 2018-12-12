#!/usr/bin/perl -w
use strict;
use Getopt::Long qw(GetOptions);
use File::Basename;
use Cwd 'abs_path';


my $bed1_file;
my $bed2_file;
my $append_bed1_file_name = 0;
my $append_bed2_file_name = 0;
my $usage = <<__EOUSAGE__;

Usage:

$0 [options]

required options:

  -a	<string>	BED1
  -b	<string>	BED2
  -c			Append BED1 file name
  -d			Apeend BED2 file name
 
__EOUSAGE__



GetOptions (
  'a=s' => \$bed1_file,
  'b=s' => \$bed2_file,
  'c'	=> \$append_bed1_file_name,
  'd'	=> \$append_bed2_file_name
);

unless ($bed1_file && $bed2_file) {
    die "Please input BEDs\n$usage";
}

open BED1, "$bed1_file" || print STDERR "cannot read file: $bed1_file";
open BED2, "$bed2_file" || print STDERR "cannot read file: $bed2_file";

my $first_line1 = <BED1>;
my $first_line2 = <BED2>;
chomp $first_line1;
chomp $first_line2;
my @bed_cols1 = split(/\t/, $first_line1);
my @bed_cols2 = split(/\t/, $first_line2);

my $bed_col_len1 = $#bed_cols1 + $append_bed1_file_name;
my $bed_col_len2 = $#bed_cols2 + $append_bed2_file_name;

seek BED1, 0, 0;
seek BED2, 0, 0;

while(<BED1>) {
	chomp;
	print;	
	if ($bed_col_len1 < $bed_col_len2) {
		for (my $i=0; $i<$bed_col_len2 - $bed_col_len1;$i++) {
			print "\t-";
		}
	}
	if ($append_bed1_file_name) {
		print "\t$bed1_file";
	}
	print "\n";
}

while(<BED2>) {
	chomp;
	print;	
	if ($bed_col_len1 > $bed_col_len2) {
		for (my $i=0; $i<$bed_col_len1 - $bed_col_len2;$i++) {
			print "\t-";
		}
	}
	if ($append_bed2_file_name) {
		print "\t$bed2_file";
	}
	print "\n";
}

close(BED1);
close(BED2);

