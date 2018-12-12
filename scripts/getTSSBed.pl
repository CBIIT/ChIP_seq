#!/usr/bin/perl -w
use strict;
use Getopt::Long qw(GetOptions);
use Time::Piece;
use DBI;
use File::Basename;
use Cwd;


my $in_file;
my $length;
my $out_dir = cwd();
my $usage = <<__EOUSAGE__;

Usage:

$0 [options]

required options:

  -i	Input coordinate file (tab seperated:chr->start->end->gene->strand)
  -l	TSS length  
    
__EOUSAGE__



GetOptions (
  'i=s' => \$in_file,
  'l=i' => \$length
);

unless ($in_file && $length) {
    die "Please input coordinate file and TSS length\n$usage";
}

&main();

sub main {
	open(INPUT_FILE, "$in_file") or die "Cannot open file $in_file";
	while (<INPUT_FILE>) {
	     chomp;
	     my @fields = split(/\t/);
		 my $tss_start;
		 my $tss_end;
		 if ($fields[4] eq "-") {
			$tss_start = $fields[2] - $length;
			$tss_end = $tss_start + $length;
		 }
		 if ($fields[4] eq "+") {
			$tss_end = $fields[1] + $length;
			$tss_start = $fields[1] - $length;			
		 }
		 print $fields[0]."\t".$tss_start."\t".$tss_end."\t".$fields[3]."\t".$fields[4]."\n";
	}
}