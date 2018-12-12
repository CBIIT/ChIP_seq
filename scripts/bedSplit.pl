#!/usr/bin/perl -w
use strict;
use Getopt::Long qw(GetOptions);
use Time::Piece;
use DBI;
use File::Basename;
use Cwd;


my $bed_file;
my $col;
my $no_header = 0;
my $out_dir = cwd();
my $usage = <<__EOUSAGE__;

Usage:

$0 [options]

required options:

	-i	Input BED file
	-j	The column to be splited (e.g -j 5 : split by 5th column)  
  
optional options:
	-n  No header
	-o	Output directory (default: $out_dir)
  
__EOUSAGE__



GetOptions (
  'i=s' => \$bed_file,
  'j=i' => \$col,
  'o=s' => \$out_dir,  
  'n' => \$no_header
);

unless ($bed_file && $col) {
    die "Please input BED and split column\n$usage";
}

$out_dir = &formatDir($out_dir);


&main();

sub main {
	open(INPUT_FILE, "$bed_file") or die "Cannot open file $bed_file";
	my %split_data = ();
	my $header;
	if (!$no_header) {
		$header = <INPUT_FILE>;
	}
	
	while (<INPUT_FILE>) {
	     chomp;
	     my @fields = split(/\t/);
		 $split_data{$fields[$col-1]} .= $_."\n";		 
	}
	
	while (my ($key, $value) = each(%split_data)) {
		my $out_file = $out_dir.basename($bed_file,".bed")."_".$key.".bed";
		open(OUTPUT_FILE, ">$out_file") or die "Cannot open file $out_file";
		if (!$no_header) {
			print OUTPUT_FILE $header;
		}
		print OUTPUT_FILE $value;
	}
}

sub formatDir {
    my ($dir) = @_;
    if ($dir !~ /\/$/) {
        $dir = $dir."/";
    }
    return $dir;
}

