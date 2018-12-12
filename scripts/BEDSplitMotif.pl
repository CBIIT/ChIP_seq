#!/usr/bin/perl -w
use strict;
use Getopt::Long qw(GetOptions);
use Time::Piece;
use File::Basename;
use Config::Simple;
use Cwd;

my $bed_file;
my $config_file;
my $output_dir = cwd();
$output_dir = &formatDir($output_dir)."output";
my $usage = <<__EOUSAGE__;

Usage:

$0 [options]

required options:

	-i	Input BED file
	-c	The config file
  
optional options:
	-o	Output directory (default: $output_dir)
  
__EOUSAGE__



GetOptions (
  'i=s' => \$bed_file,
  'c=s' => \$config_file,  
  'o=s' => \$output_dir
);

unless ($bed_file && $config_file) {
    die "Please input BED and config file\n$usage";
}

#global variables
my $config = new Config::Simple($config_file);

my $script_home = &formatDir(dirname($0));
$output_dir = &formatDir($output_dir);
#system("rm -rf $output_dir");
system("mkdir -p $output_dir");
system($script_home."bedSplit.pl -i $bed_file -j 4 -n -o $output_dir");
my @files = grep { -f } glob $output_dir."*.bed";
foreach my $file (@files) {	
	my $out_dir = &formatDir($output_dir.basename($file,".bed"));
	my $stdout = $out_dir."output.txt";
	my $cmd = "findMotifsGenome.pl $file ".$config->param("HOMER.genome")." $out_dir -size ".$config->param("HOMER.size")." -len ". $config->param("HOMER.len")." -preparsedDir ".$config->param("HOMER.preparsedDir")." -p ".$config->param("HOMER.p");
	print "find motif in ".basename($file)."...";
	#print "\n$cmd\n";
	system($cmd);
	print "done\n";
}

my %id_list = ();
my %motif_list = ();
my %pvalue_matrix = ();

foreach my $file (@files) {	
	my $out_dir = &formatDir($output_dir.basename($file,".bed"));
	my $homer_file = $out_dir."knownResults.txt";
	if (-e $homer_file) {						
		my $id = basename($file,".bed");
		$id_list{$id} = '';
		&readMotifFile($homer_file, $id);
	}
}

my $motif_file = $output_dir."motif.txt";
my @ids = sort {$a cmp $b} keys %id_list;
my @motifs = sort {$a cmp $b} keys %motif_list;
open (OUT_FILE, ">$motif_file");
my $header = "Motif";
foreach my $id(@ids) {
	$header .= "\t$id";
}
print OUT_FILE $header."\n";
foreach my $motif(@motifs) {
	print OUT_FILE $motif;
	foreach my $id(@ids) {
		if ($pvalue_matrix{$motif}{$id}) {
			print OUT_FILE "\t".$pvalue_matrix{$motif}{$id};
		}
		else {
			print OUT_FILE "\t0";
		}
	}
	print OUT_FILE "\n";
}

sub readMotifFile {
	my ($f, $id) = @_;
	open (MOTIF_FILE, $f);
	<MOTIF_FILE>;
    while (<MOTIF_FILE>) {
         my @fields = split(/\t/);
		 my $motif_id = $fields[0];
		 my $pvalue = $fields[3] * -1;
		 $motif_list{$motif_id} = '';
		 $pvalue_matrix{$motif_id}{$id} = $pvalue;             
    }
    close(MOTIF_FILE);
	
}

	
sub formatDir {
    my ($dir) = @_;
    if ($dir !~ /\/$/) {
        $dir = $dir."/";
    }
    return $dir;
}
