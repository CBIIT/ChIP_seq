#!/usr/bin/perl -w
use strict;
use Getopt::Long qw(GetOptions);
use Time::Piece;
use DBI;
use File::Basename;
use Cwd;

my $trv_file;
my $out_file;
my $value_file = "/data/khanlab/projects/ChIP_seq/data_by_file_type/test/data_by_file_type/ref/prcnt_gene_TR_values.list.txt";
my $value_idx = 9;
my $usage = <<__EOUSAGE__;

Usage:

$0 [options]

required options:

  -i  <string>  Input TRV file
  -o  <string>  Output file

optional options:

  -v  <string>  Value list file (default: $value_file)
  -c  <integer> Value column index (default: $value_idx)
  
__EOUSAGE__



GetOptions (
  'v=s' => \$value_file,
  'i=s' => \$trv_file,
  'o=s' => \$out_file  
);

unless ($trv_file && $out_file) {
    die "Please input and output file\n$usage";
}

&main();

sub main {
	
	my @break_points = ();
	my @break_point_values = ();
	open(VALUE_FILE, $value_file) or die "Cannot open file $value_file";
	while (<VALUE_FILE>) {
		chomp;
		push @break_points, $_;
		push @break_point_values, 0;
	}
	close(VALUE_FILE);
	open(TRV_FILE, $trv_file) or die "Cannot open file $trv_file";
	$value_idx--;
	my $total_genes = 0;
	while (<TRV_FILE>) {
		chomp;
		my @fields = split(/\t/);
		next if ($#fields < $value_idx);		
		my $value = $fields[$value_idx];
		next if ($value == 0);
		for (my $i=0;$i<=$#break_points;$i++) {
			if ($break_points[$i] >= $value) {
				$break_point_values[$i]++;
			}
		}
		$total_genes++;
	}
	close(TRV_FILE);
	open(OUT_FILE, ">$out_file") or die "Cannot open file $out_file";
	for (my $i=0;$i<=$#break_points;$i++) {
		print OUT_FILE $break_points[$i]."\t".$break_point_values[$i]."\t".($break_point_values[$i]/$total_genes)."\n";		
	}	
	close(OUT_FILE);
}



sub formatDir {
    my ($dir) = @_;
    if ($dir !~ /\/$/) {
        $dir = $dir."/";
    }
    return $dir;
}

sub runCommand {
    my ($cmd) = @_;
    system($cmd);
}

