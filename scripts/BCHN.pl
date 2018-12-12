#!/usr/bin/perl -w
use strict;
use Getopt::Long qw(GetOptions);
use Time::Piece;
use DBI;
use File::Basename;
use Cwd 'abs_path';


my $bed1_file;
my $bed2_file;
my $config_file;
my $output_dir;
my $bed_ext_len = 4000;
my $genome = "hg19";
my $homer_preparsedDir="/data/khanlab/projects/ChIP_seq/data_by_file_type/preparse";
my $homer_p = 32;
my $homer_size = 1000;
my $homer_len = 8;
my $sbatch_option = '--partition=ccr --cpus-per-task=32 --ntasks=2 --mem=32g';
my $ngsplot_option = '-GO total –SC 0,5';
my $sort_column = -1;
my $no_homer = 0;
my $no_ngsplot = 0;
my $script_home = &formatDir(dirname(abs_path($0)));
my $usage = <<__EOUSAGE__;

Usage:

$0 [options]

required options:

  -a	BED1
  -b	BED2
  -c	Config file
  -o	Output directory
  
optional options:    
  -e	BED extention length (default: $bed_ext_len, 0: no extention)  
  -s	The column to be sorted (e.g -s 5 : sort by 5th column)
  -u	sbatch option (default: $sbatch_option)
  -m	No NGSPlot
  -t	NGSPlot option (default: $ngsplot_option)
  -g	Genome (default: $genome)
 HOMER optional options:
  -n	Do not run HOMER  
  -r	Preparsed dir (default: $homer_preparsedDir)
  -p	Number of threads (default: $homer_p)
  -z	Motif search size (default: $homer_size)
  -l	Motif length (default: $homer_len)  
  
__EOUSAGE__



GetOptions (
  'a=s' => \$bed1_file,
  'b=s' => \$bed2_file,
  'c=s' => \$config_file,
  'o=s' => \$output_dir,
  'e=i' => \$bed_ext_len,  
  's=i' => \$sort_column,
  'n' => \$no_homer,
  'm' => \$no_ngsplot,
  'g=s' => \$genome,
  'r=s' => \$homer_preparsedDir,
  'p=i' => \$homer_p,
  'z=i' => \$homer_size,
  'l=i' => \$homer_len,
  'u=s' => \$sbatch_option,
  't=s' => \$ngsplot_option   
);


unless ($bed1_file && $bed2_file && $config_file) {
    die "Please input BED and config files\n$usage";
}

unless ($output_dir) {
    die "Please assign output dir\n";
}

#Global variables

$output_dir = &formatDir($output_dir);


&main();

sub main {
    &runCommand("mkdir -p $output_dir");
    my $bed1_base = basename($bed1_file,".bed");
	my $bed2_base = basename($bed2_file,".bed");
	my $intersect_file = $output_dir.$bed1_base."_inter_".$bed2_base.".bed";	
	my $bed1_only_file = $output_dir.$bed1_base."_only.bed";
	my $bed2_only_file = $output_dir.$bed2_base."_only.bed";
	&runCommand("bedtools intersect -a $bed1_file -b $bed2_file > $intersect_file") if (!-e $intersect_file);
	&runCommand("bedtools intersect -a $bed1_file -b $bed2_file -v > $bed1_only_file") if (!-e $bed1_only_file);
	&runCommand("bedtools intersect -a $bed2_file -b $bed1_file -v > $bed2_only_file") if (!-e $bed2_only_file);
	my $intersect_anno_file = $output_dir.$bed1_base."_inter_".$bed2_base.".annotation.txt";	
	my $bed1_only_anno_file = $output_dir.$bed1_base."_only.annotation.txt";
	my $bed2_only_anno_file = $output_dir.$bed2_base."_only.annotation.txt";
	my $intersect_anno_sum_file = $output_dir.$bed1_base."_inter_".$bed2_base.".summary.txt";	
	my $bed1_only_anno_sum_file = $output_dir.$bed1_base."_only.summary.txt";
	my $bed2_only_anno_sum_file = $output_dir.$bed2_base."_only.summary.txt";
	&runPeakAnnotation($intersect_file, $intersect_anno_file, $intersect_anno_sum_file);
	&runPeakAnnotation($bed1_only_file, $bed1_only_anno_file, $bed1_only_anno_sum_file);
	&runPeakAnnotation($bed2_only_file, $bed2_only_anno_file, $bed2_only_anno_sum_file);	
	
	if (!$no_homer) {
	    &runFindMotifQsub($intersect_file);
	    &runFindMotifQsub($bed1_only_file);
	    &runFindMotifQsub($bed2_only_file);
	}
	if ($sort_column != -1) {
	    my $bed1_only_sorted_file = dirname($bed1_only_file)."/".basename($bed1_only_file,".bed").".sorted.bed";
		my $bed2_only_sorted_file = dirname($bed2_only_file)."/".basename($bed2_only_file,".bed").".sorted.bed";
		&runCommand("sort -k$sort_column,$sort_column -nr $bed1_only_file > $bed1_only_sorted_file");
		&runCommand("sort -k$sort_column,$sort_column -nr $bed2_only_file > $bed2_only_sorted_file");
		$bed1_only_file = $bed1_only_sorted_file;
		$bed2_only_file = $bed2_only_sorted_file;
	}
	my $bed1_only_flanked_file = dirname($bed1_only_file)."/".basename($bed1_only_file,".bed").".$bed_ext_len.bed";
	my $bed2_only_flanked_file = dirname($bed2_only_file)."/".basename($bed2_only_file,".bed").".$bed_ext_len.bed";
	my $intersect_flanked_file = dirname($intersect_file)."/".basename($intersect_file,".bed").".$bed_ext_len.bed";
	my $combined_flanked_file = $output_dir."combined.".$bed_ext_len.".bed";
	if ($bed_ext_len == 0) {
		$bed1_only_flanked_file = $bed1_only_file;
		$bed2_only_flanked_file = $bed2_only_file;
		$intersect_flanked_file = $intersect_file;
	} else {
		&adjustBedFlank($bed1_only_file, $bed_ext_len, $bed1_only_flanked_file);
		&adjustBedFlank($bed2_only_file, $bed_ext_len, $bed2_only_flanked_file);
		&adjustBedFlank($intersect_file, $bed_ext_len, $intersect_flanked_file);
	}
	#&runCommand("cat $bed1_only_flanked_file $bed2_only_flanked_file $intersect_flanked_file > $combined_flanked_file");
	my $tmp_combined_file = $output_dir."combined.".$bed_ext_len.".tmp.bed";
	&runCommand($script_home."mergeBEDs.pl -a $bed1_only_flanked_file -b $intersect_flanked_file -c -d > $tmp_combined_file");
	&runCommand($script_home."mergeBEDs.pl -a $tmp_combined_file -b $bed2_only_flanked_file -d > $combined_flanked_file");
	my $bed1_only_conf_file = dirname($bed1_only_file)."/".basename($bed1_only_file,".bed").".$bed_ext_len.txt";
	my $bed2_only_conf_file = dirname($bed2_only_file)."/".basename($bed2_only_file,".bed").".$bed_ext_len.txt";
	my $intersect_conf_file = dirname($intersect_file)."/".basename($intersect_file,".bed").".$bed_ext_len.txt";
	my $combined_conf_file = dirname($intersect_file)."/"."combined".".$bed_ext_len.txt";
	if (!$no_ngsplot) {
		&prepareNGSconfigFile($config_file, $bed1_only_flanked_file, $bed1_only_conf_file);
		&prepareNGSconfigFile($config_file, $bed2_only_flanked_file, $bed2_only_conf_file);
		&prepareNGSconfigFile($config_file, $intersect_flanked_file, $intersect_conf_file);
		&prepareNGSconfigFile($config_file, $combined_flanked_file, $combined_conf_file);
		&runNGSPlot($bed1_only_conf_file, $bed1_base);
		&runNGSPlot($bed2_only_conf_file, $bed2_base);
		&runNGSPlot($intersect_conf_file, basename($intersect_file,".bed"));
		&runNGSPlot($combined_conf_file, "combined");
	}
}

sub formatDir {
    my ($dir) = @_;
    if ($dir !~ /\/$/) {
        $dir = $dir."/";
    }
    return $dir;
}

sub	runFindMotifQsub {
	my ($bed_file) = @_;
	my $homer_dir = $output_dir."HOMER_".basename($bed_file,".bed")."/";
	system("mkdir -p $homer_dir");
	my $log_dir = $homer_dir."pbs_log/";
	system("mkdir -p $log_dir");
	my $pbs_out = $log_dir."findMotif.".basename($bed_file,".bed").".o";
    my $pbs_err = $log_dir."findMotif.".basename($bed_file,".bed").".e";
    my $jobname = "findMotif";
    my $script_file = $script_home."findMotif.sh";
	my $sbatch_cmd = "sbatch $sbatch_option -J $jobname --output $pbs_out --error $pbs_err --export=input_file=$bed_file,genome=$genome,output_dir=$homer_dir,size=$homer_size,len=$homer_len,preparsedDir=$homer_preparsedDir,p=$homer_p $script_file";
	#my $qsub_cmd = "sbatch -N $jobname -o $pbs_out -e $pbs_err $qsub_option -v input_file=$bed_file -v genome=$genome -v output_dir=$homer_dir -v size=$homer_size -v len=$homer_len -v preparsedDir=$homer_preparsedDir -v p=$homer_p $script_file";
	my $output_file = $homer_dir."/homerResults.html";
    if (!-s $output_file) {
       &runCommand($sbatch_cmd);
    }	
}

sub runPeakAnnotation {
	my ($bed_file, $annotation_file, $annotation_summary_file) = @_;
	my $cmd = "annotatePeaks.pl $bed_file ".$genome;
    &runCommandToSTDOUT($cmd, $annotation_file, $annotation_summary_file);
}

sub adjustBedFlank {
    my ($bed_file, $flank_length, $out_file) = @_;	
    my $script_file = $script_home."bedFlank.pl";
    if (!-s $out_file) {
       &runCommand("$script_file -i $bed_file -f $flank_length -o $out_file");	   
    }
}

sub runNGSPlot {
    my ($config_file, $output_name) = @_;
	#my $script_file = $script_home."runNGSPlot.sh";
	#my $qsub_cmd = "qsub -N NGSPlot $qsub_option -v config_file=$config_file -v output_name=$output_name $script_file";
	chdir $output_dir;
	&runCommand("ngs.plot.r -G $genome -R bed -C $config_file -O $output_name $ngsplot_option");
}

sub prepareNGSconfigFile {
	my ($config_input_file, $bed_file, $config_output_file) = @_;
	open(INPUT_FILE, "$config_input_file") or die "Cannot open file $config_input_file";
	open(OUTPUT_FILE, ">$config_output_file") or die "Cannot open file $config_output_file";
	<INPUT_FILE>;
	while (<INPUT_FILE>) {
	     chomp;
	     my @fields = split(/\t/);
		 $fields[1] = $bed_file;		 
		 print OUTPUT_FILE join("\t",@fields)."\n";
	}
}

sub runCommand {
    my ($cmd) = @_;
	print $cmd."\n";
    system($cmd);
}

sub runCommandToSTDOUT {
    my ($cmd, $stdout_file, $stderr_file) = @_;
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
