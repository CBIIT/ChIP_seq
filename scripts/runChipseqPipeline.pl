#!/usr/bin/perl -w
use strict;
use Getopt::Long qw(GetOptions);
use Time::Piece;
use File::Basename;
use Config::Simple;
use Data::Dumper;
use Scalar::Util qw(looks_like_number);

my $usage = <<__EOUSAGE__;

Usage:

$0 <ChIPSeq config file> 
  
__EOUSAGE__

my $config_file = $ARGV[0];

unless ($config_file) {
    die $usage;
}

#global variables 
my $default_coltron_exp_cutoff = 2;
my $config = new Config::Simple($config_file);
my $override = $config->param("General.override");
my $home = &formatDir($config->param("General.home"));
my $sample_desc_file = $home.$config->param("General.sample_desc_file");
my $data_home = &formatDir($home.$config->param("General.data_dir"));
my $script_home = &formatDir($home.$config->param("General.script_dir"));
my $project_home = &formatDir($home.$config->param("General.project_dir"));
my $output_home = &formatDir($home.$config->param("General.output_dir"));
my $log_dir = &formatDir($home.$config->param("General.log_dir"));
my $bam_ext = $config->param("General.bam_extension");
my $current_time = localtime->strftime('%Y-%m-%d_%H:%M:%S');
my $cmd_log_file = $log_dir."cmd_$current_time.log";
my $incomplete_log_file = $log_dir."incomplete_$current_time.log";
my $default_read_length = $config->param("General.default_read_length");
my $default_lib_size = $config->param("General.default_lib_size");
my @genome_features = qw/3UTR miRNA ncRNA TTS pseudo Exon Intron Intergenic Promoter 5UTR snoRNA snRNA rRNA/;
my @incomplete_lists = (); 
#check we should use p-value or q-value cutoff. If p-value assigned, the q-value will be ignored
my @pvalue_cutoffs = $config->param("MACS.pvalue_cutoff");
my @qvalue_cutoffs = $config->param("MACS.qvalue_cutoff");
my @nfr_types = $config->param("HOMER.nfr");
my $keep_dup = $config->param("MACS.keep-dup");
if (!defined $keep_dup) {
	$keep_dup = "all";
}
my $macs_suffix = "";
if (defined $config->param("MACS.suffix")) {
	$macs_suffix = $config->param("MACS.suffix");
}
my @bin_sizes = $config->param("IGVTools.bin_size");
my @cutoffs = ();
my $isqvalue = 1;
if (@pvalue_cutoffs) {
    @cutoffs = @pvalue_cutoffs;
	$isqvalue = 0;
} else {
    @cutoffs = @qvalue_cutoffs;
}
my %exclude_chr_list = ();
my $exclude_chr_file = $config->param("General.exclude_chr");
print "exclude_chr_file: $exclude_chr_file\n";
if ($exclude_chr_file && -s $home.$exclude_chr_file) {	
	&getExcludeChrList($home.$config->param("General.exclude_chr"));
}
print "Running... for detailed log, please see: $cmd_log_file\n";
# Sample list. Array of Sample objects
my @sample_list = ();
# projects is the projects-samples mapping hash. 
my %projects = ();
&runCommand("mkdir -p $log_dir");
&main();

=nd

Main pipeline function:

1. read sample meta file
2. run IGV tools
3. call peaks (narrow and broad) using MACS and SICER
4. run enhancer pipeline
5. output peak summary
6. peak file comparison

=cut
sub main {
	@sample_list = &readSampleDescriptionFile($sample_desc_file);	
	# anno_summary is a multi-dim hash used to store the peak summary for each genome feature
	# The structure of anno_summary is:
	# $anno_summary{$project_id}{$sample_id}{$cutoff_log}{$type} where 
	#   project_id: ChIPseq project
	#   sample_id: sample ID
	#   cutoff_log: p-value or q-value cutoff (one MACS run might have multiple cutoffs)
	#   type: broad or narrow peaks calling
	my %anno_summary = ();
	foreach my $sample(@sample_list) {
        &print_log("running: ".$sample->sample_id);
		push @{$projects{$sample->project}}, $sample;
		if (!-e $sample->sample_file) {
		   print STDERR "Sample bam file: ".$sample->sample_file." not found!\n";
		   next;
		}
		# create tdf file for each sample with read extension
		&runIGVTools($config, $sample);
		#next;
		if (lc($sample->chip_target) ne 'input') {
			foreach my $nfr_type(@nfr_types) {
				if ($sample->chip_target eq $nfr_type) {
					&runHOMERfindPeak($config, $sample);
					last;
				}
			}
			if (!-e $sample->sample_file) {
			   print STDERR "Control bam file: ".$sample->control_file." not found!\n";
			   next;
			}        
			#my @bed_type = ("narrow", "broad");
			my $type = $sample->peak_type;		
			if ($type eq "narrow" || $type eq "broad") {
				&print_log("type: $type");
				my @sample_pvalues = @cutoffs;
				if ($sample->chosenP ne "" && $sample->chosenP ne ".") {
					@sample_pvalues = split(/,/, $sample->chosenP);
				}
				foreach my $cutoff (@sample_pvalues) {
					&print_log("cutoff: $cutoff");
					my $cutoff_log = sprintf("%1.0e",$cutoff);			   
					my ($bed_file, $summit_file) = &runMACS2($config, $sample, $type, $cutoff, $isqvalue);
					my $bed_anno = &annotatePeaks($config, $bed_file);
					$anno_summary{$sample->project}{$sample->sample_id}{$cutoff_log} = $bed_anno;									
					if ($sample->enhance_pipe && $sample->enhance_pipe eq "yes") {
						my $rose_out_file = $sample->output_dir."ROSE_out_$type/".$sample->sample_id."_Enhancers_withSuper.bed";
						&runROSE2_with_bedmerge($config, $sample, $bed_file, $summit_file);
						#&runROSE2($config, $sample, $bed_file, $summit_file);
					} 
					if (($type eq "broad") && ($config->param("SICER.enable") eq "true")){
						&runSICER($config, $sample);						
					}
				}
				my $cutoff_type = ($isqvalue)? "q" : "p";
				my $summary_file = dirname($sample->sample_file)."/".basename($sample->sample_file).".$cutoff_type.$type.summary";
				&outputSamplePeaksSummary($summary_file, $sample->project, $sample->sample_id, %anno_summary);
			}
			
		}
	}	
	&outputProjectPeaksSummary(%anno_summary);
    &runBedComparison($config, @sample_list);
	open IN_COMPLETE, ">$incomplete_log_file";
	print IN_COMPLETE "The following output files do not exists or have zero size:\n";
	foreach my $incomplete_file (@incomplete_lists) {
		print IN_COMPLETE $incomplete_file."\n";
	}
	close(IN_COMPLETE);
	&runCommand("chgrp -f -R khanlab $home");
	&runCommand("chmod -f -R 775 $home");
	print "Done!\n";
}

sub outputSamplePeaksSummary {
	my ($summary_file, $project_id, $sample_id, %anno_summary) = @_;
	open SUM_FILE, ">$summary_file";
	#print header
	my $total_peaks = "total_peaks";
	my $avg_peaks_length = "avg_peaks_length";
	
	my @sample_cutoffs = keys %{$anno_summary{$project_id}{$sample_id}};
	#foreach my $cutoff (@cutoffs) {
	#print Dumper(%anno_summary);
	#print Dumper($anno_summary{$project_id}{$sample_id});
	#print Dumper(keys $anno_summary{$project_id}{$sample_id});
	#print Dumper(@sample_cutoffs);
	print "sample_cutoff count: ".$#sample_cutoffs."\n";
	foreach my $cutoff (@sample_cutoffs) {
		my $cutoff_log = sprintf("%1.0e",$cutoff);
		print "$summary_file cutoff: $cutoff\n";
		print SUM_FILE "\t".$cutoff;
		$total_peaks = $total_peaks."\t".$anno_summary{$project_id}{$sample_id}{$cutoff}->total_peaks;
		$avg_peaks_length = $avg_peaks_length."\t".$anno_summary{$project_id}{$sample_id}{$cutoff}->avg_peaks_length;
	}
	print SUM_FILE "\n$total_peaks\n$avg_peaks_length\n";
	foreach my $genome_feature(@genome_features) {
		print SUM_FILE $genome_feature;		
		#foreach my $cutoff (@cutoffs) {
		foreach my $cutoff (@sample_cutoffs) {
		    my $cutoff_log = sprintf("%1.0e",$cutoff);
		    my %dist = %{$anno_summary{$project_id}{$sample_id}{$cutoff}->genome_feature_dist};
			my $dist_feature = $dist{$genome_feature};
			if (!$dist_feature) {
				$dist_feature = 0;
			}
			print SUM_FILE "\t".sprintf("%2.2f",($dist_feature/$anno_summary{$project_id}{$sample_id}{$cutoff}->total_peaks)*100)."\%";
		}
		print SUM_FILE "\n";
	}
	close(SUM_FILE);
	
}

sub outputProjectPeaksSummary {
	my (%anno_summary) = @_;	
	#print header	
	foreach my $project_id (keys %projects) {
	    my $project_dir = $project_home.$project_id;
		my $total_peaks = "total_peaks";
	    my $avg_peaks_length = "avg_peaks_length";
		&runCommand("mkdir -p $project_dir");
	    my @samples = @{$projects{$project_id}};		
		my $summary_file = $project_dir."/Peaks.summary";
		open SUM_FILE, ">$summary_file";		
		foreach my $sample (@samples) {
				my $type = $sample->peak_type;
			    next if (!$anno_summary{$project_id}{$sample->sample_id});
				foreach my $cutoff (@cutoffs) {
				    my $cutoff_log = sprintf("%1.0e",$cutoff);
				    if ($anno_summary{$project_id}{$sample->sample_id}{$cutoff_log}) {
						print SUM_FILE "\t".$sample->sample_id."(".$cutoff_log.")";
						$total_peaks = $total_peaks."\t".$anno_summary{$project_id}{$sample->sample_id}{$cutoff_log}->total_peaks;
						$avg_peaks_length = $avg_peaks_length."\t".$anno_summary{$project_id}{$sample->sample_id}{$cutoff_log}->avg_peaks_length;
					}
				}
		}	
		print SUM_FILE "\n$total_peaks\n$avg_peaks_length\n";
		foreach my $genome_feature(@genome_features) {
			    print SUM_FILE $genome_feature;
				foreach my $sample (@samples) {
                    next if (!$anno_summary{$project_id}{$sample->sample_id});
					foreach my $cutoff (@cutoffs) {
					    my $cutoff_log = sprintf("%1.0e",$cutoff);
					    if ($anno_summary{$project_id}{$sample->sample_id}{$cutoff_log}) {
							my %dist = %{$anno_summary{$project_id}{$sample->sample_id}{$cutoff_log}->genome_feature_dist};
							my $dist_feature = $dist{$genome_feature};
							if (!$dist_feature) {
								$dist_feature = 0;
							}
							print SUM_FILE "\t".sprintf("%2.2f",($dist_feature/$anno_summary{$project_id}{$sample->sample_id}{$cutoff_log}->total_peaks)*100)."\%";
						}
					}					
				}
				print SUM_FILE "\n";
		}
		close(SUM_FILE);		
	}	
}


sub runBedComparison {
    &print_log("Start bed comparison");
    my ($config, @sample_list) = @_;
	my $bed_comparison_home = $project_home;	
	my @adjust_flank = $config->param("BedComparison.adjust_flank");
	foreach my $project (keys %projects) {
	    my @samples = @{$projects{$project}};
		my $bed_comp_project_home = $bed_comparison_home.$project."/comp/";
		&runCommand("mkdir -p $bed_comp_project_home");
		foreach my $cutoff (@cutoffs) {
			for (my $i=0;$i<=$#samples-1;$i++) {
				for (my $j=$i+1;$j<=$#samples;$j++) {
					my $sample_a = $samples[$i];
					my $sample_b = $samples[$j];
					my $cutoff_log = sprintf("%1.0e",$cutoff);
					my $cutoff_type = ($isqvalue)? "q" : "p";
					my $out_dir = $bed_comp_project_home."MACS_Out_$cutoff_type"."_".$cutoff_log."/".$sample_a->sample_id."_vs_".$sample_b->sample_id;					 
					my $type_a = $sample_a->peak_type;
					my $type_b = $sample_b->peak_type;
					my $file_a = &getMACSOutputDir($sample_a, $cutoff, $isqvalue).&getMACSNoBlackListBedFileName($sample_a->sample_id, ($type_a eq "broad")); 
					my $file_b = &getMACSOutputDir($sample_b, $cutoff, $isqvalue).&getMACSNoBlackListBedFileName($sample_b->sample_id, ($type_b eq "broad"));
					print_log("Comparing: $file_a and $file_b");
					if (-e $file_a && -e $file_b) {
								    if (!-e $out_dir) {
						                my $cmd = "mkdir -p $out_dir";
						                &runCommand($cmd);
									}	
									my $intersect = $out_dir."/".$sample_a->sample_id.".".$type_a."_inter_".$sample_b->sample_id.".$type_b.bed";
									my $bed1_only = $out_dir."/".$sample_a->sample_id.".".$type_a.".vs.$type_b.only.bed";
									my $bed2_only = $out_dir."/".$sample_b->sample_id.".".$type_b."vs.$type_a.only.bed";
									my $script_file = $script_home."runBedIntersect.sh";                              
									if (!-s $bed2_only || $override) {
										&print_log("Running bed intersection:");
										&runBedIntersect($file_a, $file_b, $intersect, $bed1_only, $bed2_only);
										push (@incomplete_lists, $intersect) if (!-s $intersect);
										push (@incomplete_lists, $bed1_only) if (!-s $bed1_only);
										push (@incomplete_lists, $bed2_only) if (!-s $bed2_only);
										foreach my $flank_length (@adjust_flank) {
										   &adjustBedFlank($intersect, $flank_length, "$intersect.$flank_length.bed");
										   &adjustBedFlank($bed1_only, $flank_length, "$bed1_only.$flank_length.bed");
										   &adjustBedFlank($bed2_only, $flank_length, "$bed2_only.$flank_length.bed");
										}
									}
					}
				}
			}
		}
	}     	
}

sub print_log {
    my ($msg) = @_;
    open CMD_FILE, ">>$cmd_log_file" || print "cannot create command log file";
    print CMD_FILE "[".localtime->strftime('%Y-%m-%d %H:%M:%S')."] $msg\n";
    close(CMD_FILE);
}

sub annotatePeaks {
    my ($config, $bed_file) = @_;
    my $filename = basename($bed_file);
    my $pbs_out = $log_dir."homer.$filename.o";
    my $homer_detail_file = $bed_file.".annotation.txt";
    my $homer_summary_file = $bed_file.".annotation.summary";
    my $cmd = "annotatePeaks.pl $bed_file ".$config->param('HOMER.genome');
    if (!-s $homer_detail_file || $override) {
        &print_log("Running homer: Outputfile => $homer_detail_file");
        &runCommandToSTDOUT($cmd, $homer_detail_file, $homer_summary_file);
		push (@incomplete_lists, $homer_detail_file) if (!-s $homer_detail_file);
    }
	open FILE, "$homer_summary_file";
	my $found = 0;
	my %genome_feature_dist = ();
	while (<FILE>) {
	   chomp;
	   my @fields = split('\t');
	   if ($#fields == 5) {
		   if ($fields[2] eq "Annotation" && $fields[3] eq "Number of peaks") {
			   $found = 1;
			   next;
		   } 
	   } else {
		   last if ($found);
	   }   
	   if ($found) {
		   $genome_feature_dist{$fields[2]} = $fields[3]; 
	   }
	}
	my $total_peaks = &getFileLineNum($bed_file);
	my $avg_peaks_length = &runCommandWithReturn("awk -F \$'\t' '{sum+=(\$3-\$2);line++}END {print sum/line}' $bed_file");
	my $bed_anno = new BEDAnnotation(total_peaks=>$total_peaks, avg_peaks_length=>$avg_peaks_length, genome_feature_dist=>\%genome_feature_dist);
    return $bed_anno;	
}

sub runROSE2_with_bedmerge {
    my ($config, $sample, $bed_file, $summit_file) = @_;
    my @stitch_dists = $config->param('ROSE.stitch_distance');
	my $bed_file_no_tss = dirname($bed_file)."/".basename($bed_file,".bed")."_no_TSS.bed";
	&removePeaks($bed_file, $home.$config->param("General.tss_file"),$bed_file_no_tss);
	$bed_file = $bed_file_no_tss;
    foreach my $stitch_dist(@stitch_dists) {
	   my $merged_bed_file = dirname($bed_file)."/".basename($bed_file,".bed")."_$stitch_dist.bed";
	   my $cmd = "bedtools merge -i $bed_file -d $stitch_dist -c 4,5,6 -o distinct,sum,distinct > $merged_bed_file";
	   &runCommand($cmd);
       my $rose_out_dir = dirname($merged_bed_file)."/ROSE_out_$stitch_dist";
       my $rose_out_file = "$rose_out_dir/".$sample->sample_id."_peaks_AllEnhancers.table.txt";
	   $stitch_dist = 0;
       $cmd = $script_home."runROSE2.sh $merged_bed_file ".$config->param('ROSE.genome')." ".$sample->sample_file." ".$sample->control_file." ".$config->param('ROSE.tss_distance')." $stitch_dist $rose_out_dir $script_home";
       if ($sample->control_id eq ".") {
           $cmd = $script_home."runROSE2-nc.sh $merged_bed_file ".$config->param('ROSE.genome')." ".$sample->sample_file." ".$config->param('ROSE.tss_distance')." $stitch_dist $rose_out_dir $script_home";
       }
       if (!-s $rose_out_file || $override) {
          &print_log("Running ROSE2: Outputfile => $rose_out_file");
          &runCommand($cmd);
		  push (@incomplete_lists, $rose_out_file) if (!-s $rose_out_file);
       } 
	   my $se_bed_file = dirname($rose_out_file)."/".basename($rose_out_file,".txt").".super.bed";
	   my $re_bed_file = dirname($rose_out_file)."/".basename($rose_out_file,".txt").".regular.bed";	
	   $cmd = $script_home."roseTable2Bed.sh $rose_out_file $re_bed_file $se_bed_file";
	   &runCommand($cmd) if (!-s $se_bed_file || $override);
       my $great_file = dirname($rose_out_file)."/".basename($se_bed_file, ".bed").".GREAT.bed";;
       &convertGREAT($se_bed_file, $great_file);
	   $great_file = dirname($rose_out_file)."/".basename($re_bed_file, ".bed").".GREAT.bed";;
	   &convertGREAT($re_bed_file, $great_file);
	   &runColtron($config, $sample, $rose_out_file, $sample->sample_file, $config->param('ROSE.genome'), $rose_out_dir."/coltron/", $sample->sample_id);
       &runEDEN($config, $sample, $re_bed_file, basename($rose_out_file,".txt")."_re", 0);
	   &runEDEN($config, $sample, $se_bed_file, basename($rose_out_file,".txt")."_se", 1);
       my $summit_se_bed_file = dirname($rose_out_file)."/".$sample->sample_id."_summit_super.bed";
	   my $summit_re_bed_file = dirname($rose_out_file)."/".$sample->sample_id."_summit_regular.bed";
       $cmd = "bedtools intersect -wa -a $summit_file -b $se_bed_file > $summit_se_bed_file";
       if (!-s $summit_se_bed_file || $override) {
               &print_log("Running finding summits in super enhancers: Output file => $summit_se_bed_file");
               &runCommand($cmd);
			   push (@incomplete_lists, $summit_se_bed_file) if (!-s $summit_se_bed_file);
       }
	   if ($config->param('HOMER.find_motif')) {
           &runFindMotifQsub($config, $summit_se_bed_file, dirname($rose_out_file)."/motif_super");
	   }       
	   $cmd = "bedtools intersect -wa -a $summit_file -b $re_bed_file > $summit_re_bed_file";
       if (!-s $summit_re_bed_file || $override) {
               &print_log("Running finding summits in regular enhancers: Output file => $summit_re_bed_file");
               &runCommand($cmd);
			   push (@incomplete_lists, $summit_re_bed_file) if (!-s $summit_re_bed_file);
       }
	   if ($config->param('HOMER.find_motif')) {
           &runFindMotifQsub($config, $summit_re_bed_file, dirname($rose_out_file)."/motif_regular");
	   }
	   if ($config->param("MACS.do_filtering") eq "yes") {
		   my $fold_change_field = $config->param("MACS.fold_change_field");
		   my @fold_change_cutoffs = $config->param("MACS.fold_change_cutoffs");
		   foreach my $fold_change_cutoff (@fold_change_cutoffs) {
				my $bed_filtered = dirname($bed_file)."/".&getMACSBedFilteredFileName($sample->sample_id, 0, $fold_change_cutoff);
				&runBedFilter($bed_file,$bed_filtered, $fold_change_field, $fold_change_cutoff);
				my $rose_out_dir = dirname($bed_file)."/ROSE_out_$stitch_dist"."_$fold_change_cutoff"."X";
				my $rose_out_file = "$rose_out_dir/".$sample->sample_id."_peaks_AllEnhancers.table.txt";
				$cmd = $script_home."runROSE2.sh $bed_file ".$config->param('ROSE.genome')." ".$sample->sample_file." ".$sample->control_file." ".$config->param('ROSE.tss_distance')." $stitch_dist $rose_out_dir $script_home";
				if (!-s $rose_out_file || $override) {
				   &print_log("Running ROSE2: Outputfile => $rose_out_file");
				   &runCommand($cmd);
				   push (@incomplete_lists, $rose_out_file) if (!-s $rose_out_file);
				}
				my $se_bed_file = &runEDEN($config, $sample, $rose_out_file);
				my $summit_se_bed_file = dirname($rose_out_file)."/".$sample->sample_id."_summit_super.bed";
				$cmd = "bedtools intersect -wa -a $summit_file -b $se_bed_file > $summit_se_bed_file";
				if (!-s $summit_se_bed_file || $override) {
				   &print_log("Running finding summits in super enhancers: Outputfile => $summit_se_bed_file");
				   &runCommand($cmd);
				   push (@incomplete_lists, $summit_se_bed_file) if (!-s $summit_se_bed_file);
				}
				&runFindMotifQsub($config, $summit_se_bed_file, dirname($rose_out_file)."/motif_super");
		   }
	   }
    }
}

sub runROSE2 {
    my ($config, $sample, $bed_file, $summit_file) = @_;
    my @stitch_dists = $config->param('ROSE.stitch_distance');
    foreach my $stitch_dist(@stitch_dists) {
       my $rose_out_dir = dirname($bed_file)."/ROSE_out_$stitch_dist";
       my $rose_out_file = "$rose_out_dir/".$sample->sample_id."_peaks_AllEnhancers.table.txt";
       my $cmd = $script_home."runROSE2.sh $bed_file ".$config->param('ROSE.genome')." ".$sample->sample_file." ".$sample->control_file." ".$config->param('ROSE.tss_distance')." $stitch_dist $rose_out_dir $script_home";
       if ($sample->control_id eq ".") {
           $cmd = $script_home."runROSE2-nc.sh $bed_file ".$config->param('ROSE.genome')." ".$sample->sample_file." ".$config->param('ROSE.tss_distance')." $stitch_dist $rose_out_dir $script_home";
       }
       if (!-s $rose_out_file || $override) {
          &print_log("Running ROSE2: Outputfile => $rose_out_file");
          &runCommand($cmd);
		  push (@incomplete_lists, $rose_out_file) if (!-s $rose_out_file);
       } 
	   my $se_bed_file = dirname($rose_out_file)."/".basename($rose_out_file,".txt").".super.bed";
	   my $re_bed_file = dirname($rose_out_file)."/".basename($rose_out_file,".txt").".regular.bed";	
	   $cmd = $script_home."roseTable2Bed.sh $rose_out_file $re_bed_file $se_bed_file";
	   &runCommand($cmd) if (!-s $se_bed_file || $override);
	   &runColtron($config, $sample, $rose_out_file, $sample->sample_file, $config->param('ROSE.genome'), $rose_out_dir."/coltron/", $sample->sample_id);
       &runEDEN($config, $sample, $re_bed_file, basename($rose_out_file,".txt")."_re", 0);
	   &runEDEN($config, $sample, $se_bed_file, basename($rose_out_file,".txt")."_se", 1);
       my $summit_se_bed_file = dirname($rose_out_file)."/".$sample->sample_id."_summit_super.bed";
	   my $summit_re_bed_file = dirname($rose_out_file)."/".$sample->sample_id."_summit_regular.bed";
       $cmd = "bedtools intersect -wa -a $summit_file -b $se_bed_file > $summit_se_bed_file";
       if (!-s $summit_se_bed_file || $override) {
               &print_log("Running finding summits in super enhancers: Output file => $summit_se_bed_file");
               &runCommand($cmd);
			   push (@incomplete_lists, $summit_se_bed_file) if (!-s $summit_se_bed_file);
       }
	   if ($config->param('HOMER.find_motif')) {
           &runFindMotifQsub($config, $summit_se_bed_file, dirname($rose_out_file)."/motif_super");
	   }   
	   $cmd = "bedtools intersect -wa -a $summit_file -b $re_bed_file > $summit_re_bed_file";
       if (!-s $summit_re_bed_file || $override) {
               &print_log("Running finding summits in regular enhancers: Output file => $summit_re_bed_file");
               &runCommand($cmd);
			   push (@incomplete_lists, $summit_re_bed_file) if (!-s $summit_re_bed_file);
       }
	   if ($config->param('HOMER.find_motif')) {
           &runFindMotifQsub($config, $summit_re_bed_file, dirname($rose_out_file)."/motif_regular");
	   }
       if ($config->param("MACS.do_filtering") eq "yes") {
		   my $fold_change_field = $config->param("MACS.fold_change_field");
		   my @fold_change_cutoffs = $config->param("MACS.fold_change_cutoffs");
		   foreach my $fold_change_cutoff (@fold_change_cutoffs) {
				my $bed_filtered = dirname($bed_file)."/".&getMACSBedFilteredFileName($sample->sample_id, 0, $fold_change_cutoff);
				&runBedFilter($bed_file,$bed_filtered, $fold_change_field, $fold_change_cutoff);
				my $rose_out_dir = dirname($bed_file)."/ROSE_out_$stitch_dist"."_$fold_change_cutoff"."X";
				my $rose_out_file = "$rose_out_dir/".$sample->sample_id."_peaks_AllEnhancers.table.txt";
				$cmd = $script_home."runROSE2.sh $bed_file ".$config->param('ROSE.genome')." ".$sample->sample_file." ".$sample->control_file." ".$config->param('ROSE.tss_distance')." $stitch_dist $rose_out_dir";
				if (!-s $rose_out_file || $override) {
				   &print_log("Running ROSE2: Outputfile => $rose_out_file");
				   &runCommand($cmd);
				   push (@incomplete_lists, $rose_out_file) if (!-s $rose_out_file);
				}
				my $se_bed_file = &runEDEN($config, $sample, $rose_out_file);
				my $summit_se_bed_file = dirname($rose_out_file)."/".$sample->sample_id."_summit_super.bed";
				$cmd = "bedtools intersect -wa -a $summit_file -b $se_bed_file > $summit_se_bed_file";
				if (!-s $summit_se_bed_file || $override) {
				   &print_log("Running finding summits in super enhancers: Outputfile => $summit_se_bed_file");
				   &runCommand($cmd);
				   push (@incomplete_lists, $summit_se_bed_file) if (!-s $summit_se_bed_file);
				}
				&runFindMotifQsub($config, $summit_se_bed_file, dirname($rose_out_file)."/motif_super");
		   }
	   }
    }
}

sub runEDEN {
    my ($config, $sample, $bed_file, $prefix, $do_stitch) = @_;
    my $TPM_cutoff = $config->param("EDEN.TPM_cutoff");
    my $exp_file = $sample->exp_dir;
	return if ($exp_file eq ".");
    if ($config->param('EDEN.exp_file')) {
    	if ($config->param('EDEN.exp_file') ne "") {
    		$exp_file = $sample->exp_dir.$config->param('EDEN.exp_file');
    	}
    }
    my $tad_file = $home.$config->param('EDEN.tad_file');
	my $output_dir = dirname($bed_file);
	my $eden_out_file = &formatDir($output_dir).$prefix."_TPM$TPM_cutoff"."_multi-genes.txt"; 
	my $s_option = "";	
	$s_option = " -s ".$config->param('EDEN.super_loci_distance_cutoff') if ($do_stitch);
	my $exp_type = "EXP";
	if ($config->param('EDEN.exp_type')) {
		$exp_type = $config->param('EDEN.exp_type');
	}
    my $cmd = $script_home."EDEN.pl -e $exp_file -o $output_dir -x $prefix -t $exp_type -c -d $tad_file -b $bed_file -f ".$config->param('EDEN.TPM_cutoff')." $s_option -n ".$config->param('EDEN.nearest_gene_distance_cutoff');
    if (!-s $eden_out_file || $override) {
        &print_log("Running EDEN: Outputfile => $eden_out_file");
        &runCommand($cmd);
		push (@incomplete_lists, $eden_out_file) if (!-s $eden_out_file);
    }    
}

sub runColtron {
    my ($config, $sample, $peak_file, $bam_file, $genome, $output_dir, $name) = @_;
	my $coltron_dir = $config->param('COLTRON.path');
	my $coltron_dir_noexp = $config->param('COLTRON.path_noexp');
	my $exp_type = $config->param('COLTRON.exp_type');
	my $exp_file = $config->param('COLTRON.exp_file');	
	my $exp_cutoff = $config->param("COLTRON.exp_cutoff");
	my $nearest_tf_distance_cutoff = $config->param("COLTRON.nearest_tf_distance_cutoff");	
	&runCommand("mkdir -p $output_dir");
	my $cmd = $script_home."runColtron.sh $coltron_dir_noexp $peak_file $bam_file $genome $output_dir $name";
	my $exp_test = $sample->exp_dir;
	&print_log("Sample paired RNAseq file prefix is $exp_test");
	if ($exp_test ne '.') {
		if ( ref($exp_file) ne 'ARRAY' ) {
			#$exp_file = &formatDir($sample->exp_dir).$exp_file;
			$exp_file = $sample->exp_dir.$exp_file;
			my $coltron_exp_file = $output_dir."exp.txt";
			my $cmd_process_exp_file = 'grep -v tracking_id '.$exp_file.' | cut -f1,7,10 | perl -ane \'($chr,$start,$end)=$F[1]=~/(.*):(.*)-(.*)/;print "$F[0]\tchr$chr(.):$start-$end\t$F[2]\n"\' | grep -v ^NR_';
			if ($exp_type && $exp_type eq "TPM") {
				$cmd_process_exp_file = 'grep -v Start '.$exp_file.' | perl -ane \'print "$F[4]\tchr$F[0](.):$F[1]-$F[2]\t$F[5]\n"\' | grep -v ^NR_';
			}
			#grep -v Start Cl0036_T1R_T.gene.TPM.txt | perl -ane 'print "$F[3]\tchr$F[0](.):$F[1]-$F[2]\t$F[4]\n"'
			&runCommand("$cmd_process_exp_file > $coltron_exp_file");
			if (!$exp_cutoff) {
				$exp_cutoff = $default_coltron_exp_cutoff;
			}
			$cmd = $script_home."runColtron.sh $coltron_dir $peak_file $bam_file $genome $output_dir $name $nearest_tf_distance_cutoff $coltron_exp_file $exp_cutoff";
		}		
	}
    
	my $coltron_out_file = &formatDir($output_dir).$name."__NETWORK_SCATTER.pdf";
    if (!-s $coltron_out_file || $override) {
        &print_log("Running Coltron: Outputfile => $coltron_out_file");
        &runCommand($cmd);
		push (@incomplete_lists, $coltron_out_file) if (!-s $coltron_out_file);
    }    
}


sub runFindMotif {
    my ($config, $bed_file, $output_dir) = @_;
    my $cmd = "findMotifsGenome.pl $bed_file ".$config->param("HOMER.genome")." $output_dir -size ".$config->param("HOMER.size")." -len ". $config->param("HOMER.len")." -preparsedDir ".$home.$config->param("HOMER.preparsedDir")." -p ".$config->param("HOMER.p");
    my $output_file = $output_dir."/homerResults.html";
    if (!-s $output_file || $override) {
        &print_log("Running Motif finding: Output file => $output_file");
        &runCommand($cmd);
    }
}	

sub	runFindMotifQsub {
	my ($config, $bed_file, $output_dir) = @_;
	my $pbs_out = $log_dir."findMotif.".basename($bed_file).".o";
    my $pbs_err = $log_dir."findMotif.".basename($bed_file).".e";
    my $jobname = "findMotif";
    my $script_file = $script_home."findMotif.sh";
	my $sbatch_option = $config->param("HOMER.sbatch_option");
	my $genome = $config->param("HOMER.genome");
	my $size = $config->param("HOMER.size");
	my $len = $config->param("HOMER.len");
	my $preparsedDir = $home.$config->param("HOMER.preparsedDir");
	my $p = $config->param("HOMER.p");
	#my $qsub_cmd = "qsub -N $jobname -o $pbs_out -e $pbs_err $qsub_q -l $qsub_l -v input_file=$bed_file -v genome=$genome -v output_dir=$output_dir -v size=$size -v len=$len -v preparsedDir=$preparsedDir -v p=$p $script_file";
	my $sbatch_cmd = "sbatch $sbatch_option -J $jobname --output $pbs_out --error $pbs_err --export=input_file=$bed_file,genome=$genome,output_dir=$output_dir,size=$size,len=$len,preparsedDir=$preparsedDir,p=$p $script_file";
	my $output_file = $output_dir."/homerResults.html";
    if (!-s $output_file || $override) {
       &print_log("Running Motif finding: Output file => $output_file");
       &runCommand($sbatch_cmd);
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

sub runBedFilter {
    my ($input_bed, $output_bed, $field_num, $cutoff) = @_;
    my $cmd = "awk '\$$field_num > $cutoff' $input_bed > $output_bed";
    if (!-s $output_bed || $override) {
       &print_log("Running filtering: Output file => $output_bed");
       &runCommand($cmd);
    }
}

sub removePeaks {
    my ($bed_file, $blacklist, $output_file) = @_;
	&runCommand("bedtools intersect -a $bed_file -b $blacklist -v > $output_file");
}

sub runBedIntersect {
    my ($bed_1, $bed_2, $intersect, $bed1_only, $bed2_only) = @_;
    &runCommand("bedtools intersect -a $bed_1 -b $bed_2 > $intersect");
    &runCommand("bedtools intersect -a $bed_1 -b $bed_2 -v > $bed1_only");
    &runCommand("bedtools intersect -a $bed_2 -b $bed_1 -v > $bed2_only");
}

sub adjustBedFlank {
    my ($bed_file, $flank_length, $out_file) = @_;    
    my $script_file = $script_home."bedFlank.pl";
    if (!-s $out_file || $override) {
       &print_log("Running bed adjustment:");
       &runCommand("$script_file -i $bed_file -f $flank_length -o $out_file");
	   push (@incomplete_lists, $out_file) if (!-s $out_file);
    }
}

sub runMACS2 {
    my ($config, $sample, $type, $cutoff, $isqvalue) = @_;   #type: 0: narrow, 1: broad, isqvalue: 0: pvalue, 1: qvalue
    my $output_dir = &getMACSOutputDir($sample, $cutoff, $isqvalue);
    &runCommand("mkdir -p $output_dir");
    my $bed_file = $output_dir.&getMACSBedFileName($sample->sample_id, ($type eq "broad"));  #macs default output file name: xxxx.narrowPeak or xxx.broadPeak
	my $broad = ($type eq "broad")? "--broad": "";
	my $format = "BAM";
	my $ext_option = "--extsize ".$sample->ext_size;	
	if (lc($sample->library_type) eq "paired end") {
		$format = "BAMPE";
		$ext_option = "";
	}
	my $cutoff_type = ($isqvalue)? "--qvalue": "--pvalue";
	my $control_option = " ";
	if ($sample->control_id ne ".") {
	    $control_option = " -c ".$sample->control_file;
	}
	my $genome = $config->param("MACS.genome");
	my $genome_option = "";
	if ($genome) {
		$genome_option = "-g $genome";
	}
	
    my $cmd = "macs2 callpeak $broad -t ".$sample->sample_file.$control_option." $genome_option -f $format --name ".$sample->sample_id." --keep-dup $keep_dup --outdir $output_dir $ext_option $cutoff_type $cutoff";
    if (!-s $bed_file || $override) {
       &print_log("Running MACS2: Output file => $bed_file");
       &runCommand($cmd);
    }
	push (@incomplete_lists, $bed_file) if (!-s $bed_file);
    my $sample_file_esc = $sample->sample_id;
    #$sample_file_esc =~ s/\//\\\//g;/
	#MACS output postprocessing: 1. shorten peak_id, 2. remove unfinished chr 3. remove blacklist
	my $output_short_file = $bed_file.".bed.tmp";
    my $output_exclude_chr_file = $bed_file.".bed";
	my $output_noblack_file = $output_dir.&getMACSNoBlackListBedFileName($sample->sample_id, ($type eq "broad"));	
	my $summit_file = $output_dir.$sample->sample_id."_summits.bed";
	my $summit_exclude_chr_file = $output_dir.$sample->sample_id."_summits.bed.tmp";
	my $summit_extended_file = $output_dir.$sample->sample_id."_summits.extended.bed";
	my $summit_noblack_file = $output_dir.&getMACSNoBlackListSummitFileName($sample->sample_id, ($type eq "broad"));
	if (!-s $output_exclude_chr_file || $override) {
        &runCommand("sed 's/$sample_file_esc"."_//g' $bed_file > $output_short_file");
		&removeExcludeChrList($output_short_file, $output_exclude_chr_file);
		&runCommand("rm -rf $output_short_file");
	}
    if (!-s $output_noblack_file || $override) {
        &removePeaks($output_exclude_chr_file, $home.$config->param("General.blacklist"),$output_noblack_file);		
		push (@incomplete_lists, $output_noblack_file) if (!-s $output_noblack_file);
	}
	my $great_file = $output_dir.basename($output_noblack_file, ".bed").".GREAT.bed";
	&convertGREAT($output_noblack_file, $great_file);
    my $flank = int($config->param("HOMER.size")/2);
	if (-s $summit_file) {
		if (!-s $summit_extended_file || $override) {
			&removeExcludeChrList($summit_file, $summit_exclude_chr_file);
			&adjustBedFlank($summit_exclude_chr_file, $flank, $summit_extended_file);
			&runCommand("rm -rf $summit_exclude_chr_file");
		}
		if (!-s $summit_noblack_file || $override) {
			&removePeaks($summit_extended_file, $home.$config->param("General.blacklist"),$summit_noblack_file);			
			push (@incomplete_lists, $summit_noblack_file) if (!-s $summit_noblack_file);
		}
		&runFindMotifQsub($config, $summit_noblack_file, $output_dir."/motif");
	}
	&runEDEN($config, $sample, $output_noblack_file, basename($output_noblack_file,".bed"), 0);	
    return ($output_noblack_file, $summit_noblack_file);
}

sub runHOMERfindPeak {
    my ($config, $sample) = @_;
	my $nfr_dir = $data_home.$sample->sample_id."/HOMER_NFR/";
	my $tag_dir = $nfr_dir."tag/";	
	
	my $tag_info_file = $tag_dir."tagInfo.txt";
	if (!-s $tag_info_file || $override) {
		&runCommand("makeTagDirectory $tag_dir ".$sample->sample_file);
	}
	my $control_option = "";
	if ($sample->control_id ne ".") {
		my $control_tag_dir = $data_home.$sample->control_id."/HOMER_NFR/";
		$tag_info_file = $control_tag_dir."tagInfo.txt";
		if (!-s $tag_info_file || $override) {
			&runCommand("makeTagDirectory $control_tag_dir ".$sample->control_file);
		}
	    $control_option = " -i $control_tag_dir ";
	}
	my $frag_option = "";
	if ($sample->library_size ne ".") {
		$frag_option = " -fragLength ".$sample->library_size;
	}
	#my $original_peak_file = $tag_dir."peaks.txt";
	my $original_peak_file = $tag_dir."regions.txt";
	my $peak_file = $tag_dir."peaks_nfr.txt";
	if (!-s $peak_file || $override) {
		&print_log("Running findPeaks: Output file => $peak_file");
		&runCommand("findPeaks $tag_dir -style histone -o auto -nfr $frag_option $control_option");
		&runCommand("mv $original_peak_file $peak_file");
	}
	my $peak_no_nfr_file = $tag_dir."peaks_no_nfr.txt";
	if (!-s $peak_no_nfr_file || $override) {
		&print_log("Running findPeaks: Output file => $peak_no_nfr_file");
		&runCommand("findPeaks $tag_dir -style histone -o auto $frag_option $control_option");
		&runCommand("mv  $original_peak_file $peak_no_nfr_file"); 
        }

	my $bed_file = $data_home.$sample->sample_id."/".$sample->sample_id.".nfr.bed";
	my $bed_no_nfr_file = $data_home.$sample->sample_id."/".$sample->sample_id.".no_nfr.bed";
	if (!-s $bed_file || $override) {
		&print_log("Running pos2bed.pl: Output file => $bed_file");
		&runCommand("pos2bed.pl $peak_file | grep -v '^#' > $bed_file");
	}
	 if (!-s $bed_no_nfr_file || $override) {
                &print_log("Running pos2bed.pl: Output file => $bed_no_nfr_file");
                &runCommand("pos2bed.pl $peak_no_nfr_file | grep -v '^#' > $bed_no_nfr_file");
        }
	
	&runFindMotifQsub($config, $peak_file, $nfr_dir."/motif");
		
}

sub runSICER {
    #my ($input_dir, $input_bed, $control_bed, $input_bam, $control_bam, $output_dir, $qvalue_cutoff, $frag_size, $output_file) = @_; 
    my ($config, $sample) = @_;
	my $input_bed = &formatDir($home.$config->param("SICER.bed_dir")).basename($sample->sample_file,".bam").".bed";
	my $control_bed = &formatDir($home.$config->param("SICER.bed_dir")).basename($sample->control_file,".bam").".bed";
	&runBamToBed($sample->sample_file, $input_bed);
	&runBamToBed($sample->control_file, $control_bed);
	my $qvalue_cutoff = $config->param("SICER.qvalue_cutoff");	
	my $cutoff_log = sprintf("%1.0e",$qvalue_cutoff);
	my $output_dir = &formatDir($sample->output_dir)."SICER_Out_q_".$cutoff_log."/";
	&runCommand("mkdir -p $output_dir");	
	my $output_file = $output_dir.basename($input_bed, ".bed")."-W".$config->param("SICER.window")."-normalized.wig";
	my $cmd = join(' ',&formatDir($config->param("SICER.path"))."SICER.sh", $home.$config->param("SICER.bed_dir"), basename($input_bed), basename($control_bed), $output_dir, $config->param("SICER.genome"), $config->param("SICER.redundancy"), $config->param("SICER.window"), $sample->library_size, $config->param("SICER.genome_fraction"), $config->param("SICER.gap"), $config->param("SICER.qvalue_cutoff"));
	if ($sample->control_id eq ".") {
	    $cmd = join(' ',&formatDir($config->param("SICER.path"))."SICER-rb.sh", $home.$config->param("SICER.bed_dir"), basename($input_bed), $output_dir, $config->param("SICER.genome"), $config->param("SICER.redundancy"), $config->param("SICER.window"), $sample->library_size, $config->param("SICER.genome_fraction"), $config->param("SICER.gap"), $config->param("SICER.evalue_cutoff"));
    }
	if (!-s $output_file || $override) {
       &print_log("Running SICER:");
       &runCommand($cmd);
	   push (@incomplete_lists, $output_file) if (!-s $output_file);
    }
}

sub runSICERQsub {
    #my ($input_dir, $input_bed, $control_bed, $input_bam, $control_bam, $output_dir, $qvalue_cutoff, $frag_size, $output_file) = @_; 
    my ($config, $sample) = @_;
	my $input_bed = &formatDir($home.$config->param("SICER.bed_dir")).basename($sample->sample_file,".bam").".bed";
	my $control_bed = &formatDir($home.$config->param("SICER.bed_dir")).basename($sample->control_file,".bam").".bed";
	&runBamToBed($sample->sample_file, $input_bed);
	&runBamToBed($sample->control_file, $control_bed);
	my $sicer_path = &formatDir($config->param("SICER.path"));	
	my $output_file = $sample->output_dir.basename($input_bed)."-W".$config->param("SICER.window")."-normalized.wig";
	my $pbs_out = $log_dir."sicer.".$sample->sample_id.".o";
    my $pbs_err = $log_dir."sicer.".$sample->sample_id.".e";
    my $jobname = "sicer";
    my $script_file = $script_home."runSICER.sh";
	my $qsub_l = $config->param("HOMER.qsub_l");
	my $input_dir = $home.$config->param("SICER.bed_dir");
	my $input_bed_file = basename($input_bed);
	my $control_bed_file = basename($control_bed);	
	my $genome = $config->param("SICER.genome");
	my $redundancy = $config->param("SICER.redundancy");
	my $window = $config->param("SICER.window");
	my $frag_size = $sample->library_size;
	my $genome_fraction = $config->param("SICER.genome_fraction");
	my $gap = $config->param("SICER.gap");
	my $qvalue_cutoff = $config->param("SICER.qvalue_cutoff");	
	my $cutoff_log = sprintf("%1.0e",$qvalue_cutoff);
	my $output_dir = &formatDir($sample->output_dir)."SICER_Out_q_".$cutoff_log."/";
	&runCommand("mkdir -p $output_dir");
	my $qsub_cmd = "qsub -N $jobname -o $pbs_out -e $pbs_err -l $qsub_l -v sicer_path=$sicer_path -v input_dir=$input_dir -v input_bed=$input_bed_file -v control_bed=$control_bed_file -v genome=$genome -v output_dir=$output_dir -v redundancy=$redundancy -v window=$window -v frag_size=$frag_size -v genome_fraction=$genome_fraction -v gap=$gap -v qvalue_cutoff=$qvalue_cutoff $script_file";
    if (!-s $output_file || $override) {
       &print_log("Running SICER:");
       &runCommand($qsub_cmd);
    }
}

sub runBamToBed {
    my ($bam_file, $bed_file) = @_;
    my $cmd = "bedtools bamtobed -i $bam_file > $bed_file";
    if (!-s $bed_file || $override) {
        &print_log("Running bamtobed:"); 
        &runCommand($cmd);
		push (@incomplete_lists, $bed_file) if (!-s $bed_file);
    }
}

sub runIGVTools {
    my ($config, $sample) = @_;
    
	my $ext_option = "-e ".$sample->ext_size;
	if (lc($sample->library_type) eq "paired end") {		
		$ext_option = "--pairs";
	}
	foreach my $bin_size (@bin_sizes) {
		my $output_file = $sample->output_dir.$sample->sample_id.".$bin_size.tdf";
		my $output_file_with_dup = $sample->output_dir.$sample->sample_id.".$bin_size.dup.tdf";		
		my $cmd = 'igvtools count -w '.$bin_size.' '.$ext_option.' '.$sample->sample_file." $output_file ".$config->param("IGVTools.genome");
		my $cmd_dup = 'igvtools count -w '.$bin_size.' '.$ext_option.' --includeDuplicates '.$sample->sample_file." $output_file_with_dup ".$config->param("IGVTools.genome");
		if (!-s $output_file || $override) {
			&print_log("Running IGVTools: Output file => $output_file");
			&runCommand($cmd);
			push (@incomplete_lists, $output_file) if (!-s $output_file);
		}

		if (!-s $output_file_with_dup || $override) {
			&print_log("Running IGVTools: Output file => $output_file_with_dup");
			&runCommand($cmd_dup);
			push (@incomplete_lists, $output_file_with_dup) if (!-s $output_file_with_dup);
		}
		#create normalized TDF file
		my $output_norm_tdf_file = $sample->output_dir.$sample->sample_id.".$bin_size.RPM.tdf";
		if (!-s $output_norm_tdf_file || $override) {
			my $flag_stat_file = $sample->output_dir.$sample->sample_id.".flagstat.txt";
			print $flag_stat_file."\n";
			my $total_mapped_reads = readpipe("grep mapped $flag_stat_file | head -n 1 | cut -f1 -d ' '");			
			chomp $total_mapped_reads;
			print "total reads: ".$total_mapped_reads."\n";
			my $scale_factor = 1000000 / $total_mapped_reads;
			my $cmd_norm_tdf = $script_home."scaleTDF.pl -i $output_file -o $output_norm_tdf_file -f $scale_factor";
			&print_log("Running Scaling RPM TDF(Scale factor: $scale_factor): Output file => $output_norm_tdf_file");
			&runCommand($cmd_norm_tdf);
			push (@incomplete_lists, $output_norm_tdf_file) if (!-s $output_norm_tdf_file);
		}		
		#check if spike in folder exists. If yes then scale the TDF file
		my $spike_in_file = $sample->output_dir."SpikeIn/spike_map_summary";
		if (-s $spike_in_file) {			
			open(SPIKE_IN_SUMMARY, $spike_in_file) or die "Cannot open file $spike_in_file";
			<SPIKE_IN_SUMMARY>;
			my $summary_line = <SPIKE_IN_SUMMARY>;
			my @summary_fields = split(/\t/, $summary_line);
			my $scale_factor_raw = $summary_fields[1]/$summary_fields[2];
			my $scale_factor  = sprintf "%.4f", $scale_factor_raw;
			my $scaled_tdf = $sample->output_dir.$sample->sample_id.".$bin_size.scaled.tdf";
			if (!-s $scaled_tdf || $override) {
				my $cmd_scale_tdf = $script_home."scaleTDF.pl -i $output_norm_tdf_file -o $scaled_tdf -f $scale_factor";
				&print_log("Running Scaling TDF(Scale factor: $scale_factor): Output file => $scaled_tdf");
				&runCommand($cmd_scale_tdf);
				push (@incomplete_lists, $scaled_tdf) if (!-s $scaled_tdf);
			}
		}
		my $smooth_window = $config->param("IGVTools.smooth_window");
		my $smooth_option = "";
		if ($smooth_window) {
			if ($smooth_window > 0) {
				$smooth_option = " -w $smooth_window";
				if ($sample->exclude_view) {
					$smooth_option .= " -e ".$sample->exclude_view;
				}	
				my $output_file_smoothed = $sample->output_dir.$sample->sample_id.".$bin_size.smoothed.w$smooth_window.tdf";		
				my $output_file_with_dup_smoothed = $sample->output_dir.$sample->sample_id.".$bin_size.dup.smoothed.w$smooth_window.tdf";		
				my $cmd_smooth_tdf = $script_home."smoothTDF.pl -i $output_file_with_dup -o $output_file_with_dup_smoothed $smooth_option";
				if (!-s $output_file_with_dup_smoothed || $override) {
					&print_log("Running Smoothing TDF: Output file => $output_file_with_dup_smoothed");
					&runCommand($cmd_smooth_tdf);
					push (@incomplete_lists, $output_file_with_dup_smoothed) if (!-s $output_file_with_dup_smoothed);
                }
				$cmd_smooth_tdf = $script_home."smoothTDF.pl -i $output_file -o $output_file_smoothed $smooth_option";
				if (!-s $output_file_smoothed || $override) {
					&print_log("Running Smoothing TDF: Output file => $output_file_smoothed");
					&runCommand($cmd_smooth_tdf);
					push (@incomplete_lists, $output_file_smoothed) if (!-s $output_file_smoothed);
                }
			}
		}
				
		
	}
}

sub getSampleBamFileName {
    my ($sample_name) = @_;
    return $data_home.$sample_name."/".$sample_name.".$bam_ext";
}

sub getMACSBedFileName {
    my ($sample_name, $call_broad) = @_;
    my $type = ($call_broad)? "broad" : "narrow";
    return $sample_name."_peaks.$type"."Peak";
}

sub getMACSNoBlackListBedFileName {
    my ($sample_name, $call_broad) = @_;
    my $type = ($call_broad)? "broad" : "narrow";
    return $sample_name."_peaks.$type"."Peak.nobl.bed";
}

sub getMACSNoBlackListSummitFileName {
    my ($sample_name, $call_broad) = @_;
    my $type = ($call_broad)? "broad" : "narrow";
    return $sample_name."_summits.$type".".extended.nobl.bed";
}

sub getMACSOutputDir {
    my ($sample, $cutoff, $isqvalue) = @_;
    my $cutoff_log = sprintf("%1.0e",$cutoff);
	my $cutoff_type = ($isqvalue)? "q" : "p";
	my $is_pe = "";
	if (lc($sample->library_type) eq "paired end") {
	   $is_pe = "_PE";
	}
	my $dup_str = "";
	if ($keep_dup ne "all") {
		$dup_str = "_dup".$keep_dup;
	}
    return &formatDir($sample->output_dir)."MACS_Out_$cutoff_type"."_".$cutoff_log.$is_pe.$dup_str.$macs_suffix."/";
}

sub getMACSBedFilteredFileName {
    my ($sample_name, $call_broad, $fold_change_cutoff) = @_;
    my $type = ($call_broad)? "broad" : "narrow";
    return $sample_name.".$bam_ext"."_peaks.$type"."Peak.$fold_change_cutoff"."x.bed";
}

sub getExcludeChrList {
    my ($exlude_chr_file) = @_;
    &print_log("Reading $exlude_chr_file ...");
    open(FILE, $exlude_chr_file) or die "Cannot open file $exlude_chr_file";
    while(<FILE>){    
         chomp;
		 $exclude_chr_list{$_} = '';
	}	
}

sub convertGREAT {
	my ($input_file, $output_file) = @_;
	&runCommand("cut -f1-3 $input_file > $output_file");
}

sub removeExcludeChrList {
    my ($in_file, $out_file) = @_;
    open(IN_FILE, $in_file) or die "Cannot open file $in_file";
	open(OUT_FILE, ">$out_file") or die "Cannot open file $out_file";
    while(<IN_FILE>){    
         chomp;
		 my @fields = split(/\t/);
		 if (!exists $exclude_chr_list{$fields[0]}) {
		 	print OUT_FILE $_."\n";
		 }
		 
	}
	close(IN_FILE);
	close(OUT_FILE);	
}

sub readSampleDescriptionFile {
    my ($sample_desc_file) = @_;
    &print_log("Reading $sample_desc_file ...");
    my @sample_list;
    open(FILE, $sample_desc_file) or die "Cannot open file $sample_desc_file";
    my $header_str = <FILE>;
    chomp $header_str;
    my @headers = split(/\t/, $header_str);
    my $sample_file_index;
    my $pair_file_index;
	my $lib_type_index;
    my $lib_size_index;
	my $chiptarget_index;
    my $read_length_index;
    my $exp_dir_index;
    my $project_index;
	my $chosenP_index;
    my $enhance_pipe_index;
	my $peak_type_index;
	my $exclude_view_index;
    for (my $i=0;$i<=$#headers;$i++) {
         $sample_file_index = $i if ($headers[$i] eq 'SampleFiles');
         $pair_file_index = $i if ($headers[$i] eq 'PairedInput');
         $lib_type_index = $i if ($headers[$i] eq 'LibraryType');
		 $lib_size_index = $i if ($headers[$i] eq 'LibrarySize');
         $read_length_index = $i if ($headers[$i] eq 'ReadLength');
         $exp_dir_index = $i if ($headers[$i] eq 'PairedExpression');
		 $chiptarget_index = $i if ($headers[$i] eq 'ChIPtarget');
		 $chosenP_index = $i if ($headers[$i] eq 'ChosenP');
         $project_index = $i if ($headers[$i] eq 'Project');
         $enhance_pipe_index = $i if ($headers[$i] eq 'EnhancePipe');
		 $peak_type_index = $i if ($headers[$i] eq 'PeakCalling');
		 $exclude_view_index = $i if ($headers[$i] eq 'ExcludeView');
    }
    while(<FILE>){    
         chomp;
		 next if (/^#/);
         my @fields = split(/\t/);
         next if ($#fields != $#headers);
         my $sample_id = $fields[$sample_file_index];
         my $control_id = $fields[$pair_file_index];
         $control_id =~ s/\*/Sample_/;
		 my $library_type = $fields[$lib_type_index];
		 my $chip_target = $fields[$chiptarget_index];
         my $library_size = $fields[$lib_size_index];
		 if (!looks_like_number($library_size)) {
		     print_log("Cannot find library_size info. Use default $default_lib_size instead.");
		     $library_size = $default_lib_size;
         }			 
         my $read_length = $fields[$read_length_index];
		 if (!looks_like_number($read_length)) {
		     print_log("Cannot find read length info. Use default $default_read_length instead.");
		     $read_length = $default_read_length;
         }
         my $ext_size = $library_size - $read_length;
         my $exp_dir = $fields[$exp_dir_index];
         my $project = $fields[$project_index];
		 my $chosenP = ".";
		 $chosenP = $fields[$chosenP_index] if ($chosenP_index);		 
         my $enhance_pipe = $fields[$enhance_pipe_index];
		 my $peak_type = $fields[$peak_type_index];
		 my $exclude_view;
		 $exclude_view = $fields[$exclude_view_index] if ($exclude_view_index);
         my $output_dir = $output_home.$sample_id."/";
         my $sample_file = &getSampleBamFileName($sample_id);
         print $sample_file."\n";
         my $sample_file_base = basename($sample_file, ".bam");
         my $control_file = &getSampleBamFileName($control_id);
         my $control_file_base = basename($control_file, ".bam");
		 push (@incomplete_lists, $sample_file) if (!-s $sample_file);
		 push (@incomplete_lists, $sample_file) if (!-s $control_file && $control_id ne '.');
         my $sample = new Sample(sample_id=>$sample_id, sample_file=>$sample_file, read_length=>$read_length, library_type=>$library_type, library_size=>$library_size, ext_size=>$ext_size, control_id=>$control_id, control_file=>$control_file, output_dir=>$output_dir, project=>$project, chosenP=> $chosenP, chip_target=>$chip_target, enhance_pipe=>$enhance_pipe, peak_type=>$peak_type, exp_dir=>$exp_dir, exclude_view=>$exclude_view);
         push @sample_list, $sample;
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

package Sample;

use Class::MethodMaker
    get_set       => [qw(sample_id sample_file read_length library_type library_size ext_size control_id control_file output_dir project chosenP chip_target enhance_pipe peak_type exp_dir exclude_view)],
    new_hash_init => 'new';

package BEDAnnotation;

use Class::MethodMaker
	get_set       => [qw(total_peaks avg_peaks_length genome_feature_dist)],
    new_hash_init => 'new';
