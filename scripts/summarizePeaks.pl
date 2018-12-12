#!/usr/bin/perl -w
use strict;
use Getopt::Long qw(GetOptions);
use Time::Piece;
use File::Basename;
use Config::Simple;
use Scalar::Util qw(looks_like_number);

my $data_home = "/data/khanlab/projects/ChIP_seq/DATA";
my $sample_file;
my $usage = <<__EOUSAGE__;

Usage:

$0 [options]

  -d  <string>  Data home directory (default: $data_home)
  -s  <string>  Sample description file  
  
__EOUSAGE__



GetOptions (
  'd=s' => \$data_home,
  's=s' => \$sample_file
);


unless ($sample_file) {
    die $usage;
}

&main();

sub main {
	my $data_home = &formatDir($data_home);	
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
				my $macs_out_file = $macs_dir.$sample."_peaks.narrowPeak.nobl.bed";
				if (-e $macs_out_file) {
					my $num_peak = &getFileLineNum($macs_out_file);
					$sum_macs{$sample}{$p_value} = $num_peak;
					opendir(DIR, $macs_dir);
					my @rose_dirs = grep(/ROSE_out_/,readdir(DIR));
					closedir(DIR);
					foreach my $rose_dir (@rose_dirs) {
						my ($rose_stitch_len)=$rose_dir=~/ROSE_out_(.*)/;			
						if ($rose_stitch_len) {
							$rose_stitch_lens{$rose_stitch_len} = '';
							$rose_dir = $macs_dir.$rose_dir."/";							
							my $re_file = $rose_dir.$sample."_peaks_AllEnhancers.table.regular.bed";
							my $se_file = $rose_dir.$sample."_peaks_AllEnhancers.table.super.bed";
							my $total = 0;
							my $re = 0;
							my $se = 0;
							if (-e $re_file) {
									$re = &getFileLineNum($re_file);
									$total += $re;
									#print $rose_stitch_len."\t".$num_peak."\n";
							}
							if (-e $se_file) {
									$se = &getFileLineNum($se_file);
									$total += $se;
									#print $rose_stitch_len."\t".$num_peak."\n";
							}							
							$sum_rose{$sample}{$p_value}{$rose_stitch_len} = "$total\t$re\t$se";
						}
					}
				}
			}
		}       		
	}
	my @rose_stitch_len_list = sort {$a <=> $b} keys %rose_stitch_lens;
	my $header = "Sample\tp-value\tMACS";
	foreach my $rose_stitch_len (@rose_stitch_len_list) {
		$header .= "\tTotal-$rose_stitch_len\tRegular-$rose_stitch_len\tSuper-$rose_stitch_len";
	}
	$header .= "\n";
	print $header;
	while (my ($sample, $macs) = each %sum_macs) {
		while (my ($pvalue, $num_peak) = each %{$macs}) {
			my $rose_results = "";
			foreach my $rose_stitch_len (@rose_stitch_len_list) {
				my $rose_result = $sum_rose{$sample}{$pvalue}{$rose_stitch_len};
				if ($rose_result) {
					$rose_results .= "\t".$rose_result;
				} else {
					$rose_results .= "\t0\t0\t0";
				}
				if ($sum_rose{$sample}{$pvalue}{$rose_stitch_len}) {
					$rose_result .= "\t";
				}
			}
			print "$sample\t$pvalue\t$num_peak$rose_results\n";
		}
	}
	
}




sub runCommand {
    my ($cmd) = @_;
    &print_log($cmd);
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

