#!/usr/bin/perl -w
use strict;
use Getopt::Long qw(GetOptions);
use Time::Piece;
use File::Basename;
use Config::Simple;
use Scalar::Util qw(looks_like_number);
use Cwd;
	
my $data_home = "/data/khanlab/projects/ChIP_seq/DATA";
my $sample_file;
my $usage = <<__EOUSAGE__;

Usage:

$0 [options]

  -d  <string> Data home dir (default: $data_home)
  -s  <string> Sample description file  
  
__EOUSAGE__



GetOptions (
  'd=s' => \$data_home,
  's=s' => \$sample_file
);


unless ($sample_file) {
    die $usage;
}

my %id_list = ();
my %motif_list = ();
my %pvalue_matrix = ();

&main();

sub main {
	$data_home = &formatDir($data_home);
	my @sample_list = &readSampleDescriptionFile($sample_file);
	my %sum_macs = ();
	my %sum_rose = ();
	my %rose_stitch_lens = ();
	
	foreach my $sample(@sample_list) {
		my $input_dir = $data_home.$sample."/";
		opendir(DIR, $input_dir);
		my @macs_dirs = grep(/MACS_Out_/,readdir(DIR));
		closedir(DIR);
		foreach my $macs_dir (@macs_dirs) {				
			my ($p_value)=$macs_dir=~/MACS_Out_p_(.*)/;			
			if ($p_value) {
					$macs_dir = $input_dir.$macs_dir."/";
					my $homer_file = $macs_dir."motif/knownResults.txt";
					if (-e $homer_file) {						
						my $id = "$sample.$p_value";
						$id_list{$id} = '';
						&readMotifFile($homer_file, $id);
					}
					opendir(DIR, $macs_dir);
					my @rose_dirs = grep(/ROSE_out_/,readdir(DIR));
					closedir(DIR);
					foreach my $rose_dir (@rose_dirs) {
						my ($rose_stitch_len)=$rose_dir=~/ROSE_out_(.*)/;			
						if ($rose_stitch_len) {
							$rose_dir = $macs_dir.$rose_dir."/";
							$homer_file = $rose_dir."motif_super/knownResults.txt";
							if (-e $homer_file) {
								my $id = "$sample.$p_value.$rose_stitch_len";
								$id_list{$id} = '';
								&readMotifFile($homer_file, $id);
							}					
						}
					}				
			}
		}		
	}
	my @ids = sort {$a cmp $b} keys %id_list;
	my @motifs = sort {$a cmp $b} keys %motif_list;
	my $header = "Motif";
	foreach my $id(@ids) {
			$header .= "\t$id";
	}
	print $header."\n";
	foreach my $motif(@motifs) {
		print $motif;
		foreach my $id(@ids) {
			if ($pvalue_matrix{$motif}{$id}) {
				print "\t".$pvalue_matrix{$motif}{$id};
			}
			else {
				print "\t0";
			}
		}
		print "\n";
	}
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

sub runCommand {
    my ($cmd) = @_;
    system($cmd);
}

sub runCommandWithReturn {
    my ($cmd) = @_;
    my $result = readpipe($cmd);
	chomp $result;
    return $result;
}

sub runCommandToSTDOUT {
    my ($cmd, $stdout_file, $stderr_file) = @_;
    &print_log($cmd);
    open OUT_FILE, ">$stdout_file" || print STDERR "cannot create file: $stdout_file";
    open ERR_FILE, ">$stderr_file" || print STDERR "cannot create file: $stderr_file";
    open (CMD, "($cmd | sed 's/^/STDOUT:/') 2>&1 |");
    while (<CMD>) {
         if (s/^STDOUT://)  {
             print OUT_FILE $_;
         } else {
             print ERR_FILE $_;
         }
    }
    close(OUT_FILE);
    close(ERR_FILE);

}

sub readSampleDescriptionFile {
    my ($sample_desc_file) = @_;
    my @sample_list;
    open(FILE, $sample_desc_file) or die "Cannot open file $sample_desc_file";
    my $header_str = <FILE>;
    chomp $header_str;
    my @headers = split(/\t/, $header_str);
    my $sample_file_index;
    for (my $i=0;$i<=$#headers;$i++) {
         $sample_file_index = $i if ($headers[$i] eq 'SampleFiles');
    }
    while(<FILE>){    
         chomp;
		 next if (/^#/);
         my @fields = split(/\t/);
         next if ($#fields != $#headers);
         my $sample_id = $fields[$sample_file_index];
         push @sample_list, $sample_id;
    }
    close(FILE);
    return @sample_list;
}

sub formatDir {
    my ($dir) = @_;
    if ($dir !~ /\/$/) {
        $dir = $dir."/";
    }
    return $dir;
}
sub getFileLineNum {
	my ($file) = @_;
	my $result = &runCommandWithReturn("wc -l $file");
	return (split(/\s/,$result))[0];
}

