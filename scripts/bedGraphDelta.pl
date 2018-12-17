#!/usr/bin/perl -w
use strict;
use Getopt::Long qw(GetOptions);
use Time::Piece;
use DBI;
use File::Basename;
use Cwd;


my $in_file1;
my $in_file2;
my $out_file;
my $out_dir = cwd();
my $usage = <<__EOUSAGE__;

Usage:

$0 [options]

required options:

  -1	first bedgraph file
  -2	second bedgraph file
  -o    output bedgraph file
    
__EOUSAGE__



GetOptions (
  '1=s' => \$in_file1,
  '2=s' => \$in_file2,
  'o=s' => \$out_file
);

unless ($in_file1 && $in_file1 && $out_file) {
    die "Please input required files\n$usage";
}

&main();

sub main {
	open(INPUT_FILE_1, "$in_file1") or die "Cannot open file $in_file1";
	open(INPUT_FILE_2, "$in_file2") or die "Cannot open file $in_file2";
	open(OUTPUT_FILE, ">$out_file") or die "Cannot open file $out_file";
	while (<INPUT_FILE_1>) {
	     chomp;
		 my $line = <INPUT_FILE_2>;
		 chomp $line;
	     my @field1 = split(/\t/);
		 my @field2 = split(/\t/, $line);
		 if ($field1[0] ne $field2[0] || $field1[1] ne $field2[1] || $field1[2] ne $field2[2]) {
			close(OUTPUT_FILE);
			die "Two input bedgraph files have different bin regions";
		 }
		 print OUTPUT_FILE $field1[0]."\t".$field1[1]."\t".$field1[2]."\t".abs($field1[3]-$field2[3])."\n";		 
	}
}
close(INPUT_FILE_1);
close(INPUT_FILE_2);
close(OUTPUT_FILE);