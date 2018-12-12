#!/usr/bin/perl -w
use strict;
use Getopt::Long qw(GetOptions);
use Time::Piece;
use DBI;
use File::Basename;
use Cwd;

my $input_file;
my $input_dir = "/data/khanlab/projects/ChIP_seq/DATA/";
my $output_dir = "/data/khanlab/projects/ChIP_seq/DATA/";
my $mismatch_allowed = 1;
my $check_both = 0;
#my $output_dir = cwd();

my $usage = <<__EOUSAGE__;

Usage:

$0 [options]

required options:

	-i <string> Input index file  

optional options:    
	-d <string>  Input directory (default: $input_dir)
	-o <string>  Output directory (default: $output_dir)
	-m <integer> Mismatches allowed (default: $mismatch_allowed)
	-b           Check both direction
  
__EOUSAGE__



GetOptions (
	'i=s' => \$input_file,
	'd=s' => \$input_dir,
	'o=s' => \$output_dir,
	'm=i' => \$mismatch_allowed,
	'b' => \$check_both
);

unless ($input_file) {
	die "Please input index file\n$usage";
}

&main();

sub main {

	$output_dir = &formatDir($output_dir);
	$input_dir = &formatDir($input_dir);

	open(IN_FILE, $input_file) or die "Cannot open file $input_file";	
	
	while(<IN_FILE>) {
		next if (/^#/);
		chomp;
		my ($sample, $out_name, $first_cut_seq, $first_strand, $second_cut_seq, $second_strand) = split(/\t/);
		my $fastq_file = $input_dir.$sample."/".$sample."_R1.fastq.gz";
		my $output_first_dir = $output_dir.$out_name."_FirstCut";
		my $output_second_dir = $output_dir.$out_name."_SecondCut";
		my $output_both_dir = $output_dir.$out_name."_Both";
		my $output_neither_dir = $output_dir.$out_name."_Neither";
		my $output_first_file = $output_first_dir."/".$out_name."_FirstCut_R1.fastq";
		my $output_second_file = $output_second_dir."/".$out_name."_SecondCut_R1.fastq";
		my $output_both_file = $output_both_dir."/".$out_name."_Both_R1.fastq";
		my $output_neither_file = $output_neither_dir."/".$out_name."_Neither_R1.fastq";
		system("mkdir -p $output_first_dir");
		system("mkdir -p $output_first_dir"."/pbs_log");
		system("mkdir -p $output_second_dir");
		system("mkdir -p $output_second_dir"."/pbs_log");
		system("mkdir -p $output_both_dir");
		system("mkdir -p $output_both_dir"."/pbs_log");
		system("mkdir -p $output_neither_dir");
		system("mkdir -p $output_neither_dir"."/pbs_log");
		my $first_cut_rc_seq = &reverse_complement($first_cut_seq);
		my $second_cut_rc_seq = &reverse_complement($second_cut_seq);
		($first_cut_seq, $first_cut_rc_seq) = ($first_cut_rc_seq, $first_cut_seq) if ($first_strand eq 'R');
		($second_cut_seq, $second_cut_rc_seq) = ($second_cut_rc_seq, $second_cut_seq) if ($second_strand eq 'R');
		#print "second strand: $second_strand\n";
		#print "second_cut_seq: $second_cut_seq\n";
		#print "second_cut_rc_seq: $second_cut_rc_seq\n";
		open(FASTQ, "gunzip -c $fastq_file |") or die "cannot open pipe to $fastq_file: $~\n";
		open(FIRST, ">$output_first_file") or die "cannot open $output_first_file\n";
		open(SECOND, ">$output_second_file") or die "cannot open $output_second_file\n";
		open(BOTH, ">$output_both_file") or die "cannot open $output_both_file\n";
		open(NEITHER, ">$output_neither_file") or die "cannot open $output_neither_file\n";
		print "processing sample: $sample\n";
		my $total_reads = 0;
		my $first_reads = 0;
		my $second_reads = 0;
		my $both_reads = 0;
		while (<FASTQ>) {
			my $line1 = $_;
			my $line2 = <FASTQ>;
			my $line3 = <FASTQ>;
			my $line4 = <FASTQ>;
			chomp $line2;
			chomp $line4;
			my ($trimmed_line2, $trimmed_line4) = &remove_viewpoint($first_cut_seq, $line2, $line4);
			
			if ($trimmed_line2 eq "" && $check_both) {
				($trimmed_line2, $trimmed_line4) = &remove_viewpoint($first_cut_rc_seq, $line2, $line4);
			}
			my $found_first = 0;
			my $found_second = 0;
			if ($trimmed_line2 ne "") {
				print FIRST $line1.$trimmed_line2."\n".$line3.$trimmed_line4."\n";
				print BOTH $line1.$trimmed_line2."\n".$line3.$trimmed_line4."\n";
				$first_reads++;
				$both_reads++;
				$found_first = 1;
			}
			($trimmed_line2, $trimmed_line4) = &remove_viewpoint($second_cut_seq, $line2, $line4);
			if ($trimmed_line2 eq "" && $check_both) {
				($trimmed_line2, $trimmed_line4) = &remove_viewpoint($second_cut_rc_seq, $line2, $line4);
			}
			if ($trimmed_line2 ne "") {
				print SECOND $line1.$trimmed_line2."\n".$line3.$trimmed_line4."\n";
				print BOTH $line1.$trimmed_line2."\n".$line3.$trimmed_line4."\n";
				$second_reads++;
				$both_reads++;
				$found_second = 1;
			}
			if (!$found_first && !$found_second) {
				print NEITHER $line1.$line2."\n".$line3.$line4."\n";
			}
			$total_reads++;			
		}
		close(FASTQ);
		close(FIRST);
		close(SECOND);
		close(BOTH);
		print "total reads: $total_reads\n";
		print "first cut reads: $first_reads (".sprintf("%.1f",($first_reads/$total_reads)*100)."%)\n";
		print "second cut reads: $second_reads (".sprintf("%.1f",($second_reads/$total_reads)*100)."%)\n";
		print "both reads: $both_reads (".sprintf("%.1f",($both_reads/$total_reads)*100)."%)\n";
		system("gzip -f $output_first_file");
		system("gzip -f $output_second_file");
		system("gzip -f $output_both_file");
		system("gzip -f $output_neither_file");
		system("chmod -R 775 $output_first_dir");
		system("chgrp -R khanlab $output_first_dir");
		system("chmod -R 775 $output_second_dir");
		system("chgrp -R khanlab $output_second_dir");
		system("chmod -R 775 $output_both_dir");
		system("chgrp -R khanlab $output_both_dir");
		system("chmod -R 775 $output_neither_dir");
		system("chgrp -R khanlab $output_neither_dir");
	}
	close(IN_FILE);
}

sub reverse_complement {
	my $dna = shift;
	my $revcomp = reverse($dna);
	$revcomp =~ tr/ACGTacgt/TGCAtgca/;
	return $revcomp;
}

sub remove_viewpoint_old {
	my ($viewpoint, $viewpoint_rc, $seq, $qual) = @_;
	my $idx_f = index($seq, $viewpoint);
	my $idx_r = index($seq, $viewpoint_rc);	
	if ($idx_f != -1) {		
		return (substr($seq, $idx_f+length($viewpoint)), substr($qual, $idx_f+length($viewpoint)));
	}
	if ($idx_r != -1) {
		return (substr($seq, 0, $idx_r+length($viewpoint_rc)), substr($qual, $idx_r+length($viewpoint_rc)));
	}	
	return ("","");
}

sub remove_viewpoint {
	my ($viewpoint, $seq, $qual) = @_;
	#check forward strand of the read
	if (&strMatch($seq, $viewpoint)) {		
		return (substr($seq, length($viewpoint)), substr($qual, length($viewpoint)));
	} else {
		#check reverse strand of the read
		my $seq_rc = &reverse_complement($seq);
		if (&strMatch($seq_rc, $viewpoint)) {
			return (substr($seq, 0, length($seq) - length($viewpoint) + 1), substr($qual, 0, length($seq) - length($viewpoint) + 1));
		}	
	}
	return ("","");
}

sub formatDir {
    my ($dir) = @_;
    if ($dir !~ /\/$/) {
        $dir = $dir."/";
    }
    return $dir;
}

sub strMatch {
	my ($seq, $viewpoint) = @_;
	my @seq_arr = split(//, $seq);
	my @viewpoint_arr = split(//, $viewpoint);
	my $total_mismatch = 0;
	for (my $i=0;$i<length($viewpoint);$i++) {
		if ($seq_arr[$i] ne $viewpoint_arr[$i]) {
			$total_mismatch++;
		}
		if ($total_mismatch > $mismatch_allowed) {
			return 0;
		}
	}
	return 1;
}
