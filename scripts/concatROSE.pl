#!/usr/bin/perl -w
use strict;
use Getopt::Long qw(GetOptions);
use Time::Piece;
use File::Basename;
use Config::Simple;
use Scalar::Util qw(looks_like_number);
use Cwd;
	
my $data_home = "/data/khanlab/projects/ChIP_seq/DATA";
my $output_dir = getcwd;
my $sample_file;
my $usage = <<__EOUSAGE__;

Usage:

$0 [options]

  -d  <string> Data home dir (default: $data_home)
  -s  <string> Sample description file  
  -o  <string> Output dir (default: $output_dir)
  
__EOUSAGE__



GetOptions (
  'd=s' => \$data_home,
  's=s' => \$sample_file,
  'o=s' => \$output_dir
);


unless ($sample_file) {
    die $usage;
}

&main();

sub main {
	$data_home = &formatDir($data_home);
	$output_dir = &formatDir($output_dir);
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
					opendir(DIR, $macs_dir);
					my @rose_dirs = grep(/ROSE_out_/,readdir(DIR));
					closedir(DIR);
					foreach my $rose_dir (@rose_dirs) {
						my ($rose_stitch_len)=$rose_dir=~/ROSE_out_(.*)/;			
						if ($rose_stitch_len) {
							$rose_stitch_lens{$rose_stitch_len} = '';
							$rose_dir = $macs_dir.$rose_dir."/";							
							my $re_file = $rose_dir.$sample."_peaks_AllEnhancers.table.txt";
							my $se_file = $rose_dir.$sample."_peaks_SuperEnhancers_ENHANCER_TO_GENE.txt";
							my $out_file = $output_dir."$sample.$p_value.$rose_stitch_len.combined.txt";
							print $re_file."\n";
							print $se_file."\n";
							print $out_file."\n";
							if (-e $re_file && -e $se_file) {
									my $cmd = "cut -f1-12,14 $se_file > $out_file;awk -F'\\t' '{if (\$NF==0) {for (i=1;i<NF;i++){printf \"%s\\t\",\$i;}printf \"\\t\\t\\t%s\\n\",\$NF;}}' $re_file >> $out_file";
									print $cmd."\n";
									&runCommand($cmd);
							}							
						}
					}				
			}
		}       		
	}
	
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

