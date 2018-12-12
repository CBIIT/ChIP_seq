#!/usr/bin/perl -w
use strict;
use Getopt::Long qw(GetOptions);
use Time::Piece;
use DBI;
use File::Basename;
use Cwd;

my $tad_file;
my $exp_file;
my $exp_type;
my $input_file;
my $output_dir = cwd();
my $original_column = "";
my $output_prefix;
my $show_exp = 0;
my $do_stitch = 0;
my $fpkm_cutoff = 1;
my $lfc_cutoff = 1;
my $pvalue_cutoff;
my $qvalue_cutoff = 0.05;
my $super_loci_cutoff;
my $nearest_gene_loci_cutoff = 500000;
my $split_col;
my $max_dist = 249250621;              #this is the max chromosome length
my $gene_coor_file = dirname($0)."/gene_coordinates.txt";
my $matrix_col = 1;
my $matrix_header = "";
my $usage = <<__EOUSAGE__;

Usage:

$0 [options]

required options:

  -e  <string>  expression file
  -t  <string>  Type (Cufflink expression: 'EXP' , Khanlab TPM: 'TPM', Cuffdiff: 'DIFF', Matrix file: 'MATRIX')
  -b  <string>  Input BED file  

optional options:  
  -m  <integer> mth column in matrix file (start with 0. default: $matrix_col)  
  -d  <string>  TAD file
  -f  <float>   FPKM/TPM cutoff (default: $fpkm_cutoff)
  -l  <float>   Log fold change cutoff for differentially expressed genes(absolute value, default: $lfc_cutoff)
  -p  <float>   p-value cutoff for differentially expressed genes
  -q  <float>   q-value cutoff for differentially expressed genes (default: $qvalue_cutoff. Ignored if -p specified)
  -s  <integer> Super loci length cutoff
  -o  <string>  Output directory (default: $output_dir) 
  -n  <integer> The distance to nearest gene (default: $nearest_gene_loci_cutoff)
  -c            Output expression columns (cuffdiff columns if -t = 'DIFF')
  -a  <string>  Output columns in BED file (e.g.,'4,6' will output 4th and 6th column)
  -j  <integer> The column to be splited (e.g -j 5 : split by 5th column)
  -x  <string>  Output file prefix
  -r  <string>  Coordinate file (default: $gene_coor_file)
  
__EOUSAGE__



GetOptions (
  'e=s' => \$exp_file,
  'm=i' => \$matrix_col,
  't=s' => \$exp_type,
  'b=s' => \$input_file,
  'o=s' => \$output_dir,
  'x=s' => \$output_prefix,
  'd=s' => \$tad_file,
  'f=f' => \$fpkm_cutoff,
  'l=f' => \$lfc_cutoff,
  'p=f' => \$pvalue_cutoff,
  'q=f' => \$qvalue_cutoff,
  's=i' => \$super_loci_cutoff,
  'n=i' => \$nearest_gene_loci_cutoff,
  'j=i' => \$split_col,
  'a=s' => \$original_column,
  'r=s' => \$gene_coor_file,
  'c' => \$show_exp  
);

if (!$exp_type || ($exp_type ne 'EXP' && $exp_type ne 'TPM' && $exp_type ne 'DIFF' && $exp_type ne 'MATRIX')) {
    die "The t option must be either 'EXP' OR 'DIFF' or 'MATRIX'\n$usage";
}
unless ($exp_file && $input_file) {
    die "Please input expression file and BED file\n$usage";
}

#Global variables
my %exp = ();
my %tad = ();
my %tss = ();
my %coord = ();
my @beds = ();
my $cmd_log_file;
my @original_columns = split(/,/,$original_column);
my $original_column_header = "";
foreach my $col (@original_columns) {
	    $original_column_header .= "column_$col\t";
}
&main();

sub main {
	
	%tad = &readTAD($tad_file) if ($tad_file);

	&readCoordinate($gene_coor_file);	
	
	$output_dir = &formatDir($output_dir);
	
	if (!$output_prefix) {
		$output_prefix = basename($input_file, "bed");
	}
	
	runCommand("mkdir -p $output_dir");
	if ($split_col) {
		runCommand(dirname($0)."/bedSplit.pl -i $input_file -j $split_col -o $output_dir");
		my $file_pattern = $output_dir.basename($input_file,".bed")."_*.bed";
		my @files = glob($file_pattern);
		foreach my $file (@files) {
			$output_prefix = basename($file,".bed");
			&runEDEN($file, $output_dir, $output_prefix);
		}
	} else {
		&runEDEN($input_file, $output_dir, $output_prefix);
	}
	
	
}

sub runEDEN {
    my ($input_file, $output_dir, $output_prefix) = @_;
	# Output files	
	
	my $is_stitched = defined($super_loci_cutoff);	
	
	@beds = readBED($input_file);
    my @fpkm_cutoffs = ();
    push(@fpkm_cutoffs, 0);
    push(@fpkm_cutoffs, $fpkm_cutoff);
	foreach my $f (@fpkm_cutoffs) {
	    $fpkm_cutoff = $f;
	    my $multi_gene_file = $output_dir.$output_prefix."_fpkm$fpkm_cutoff"."_multi-genes.txt";
		my $nearest_gene_file = $output_dir.$output_prefix."_fpkm$fpkm_cutoff"."_nearest-genes.txt";
		my $max_gene_file = $output_dir.$output_prefix."_fpkm$fpkm_cutoff"."_max-genes.txt";
		%exp = &readExpression($exp_file, $exp_type);
        &findNearestGene($multi_gene_file, $original_column_header, @beds);
	    &findNearestMaxGeneWithDistance($nearest_gene_file, $max_gene_file, 0, $original_column_header, @beds);
	    if ($is_stitched) {
	        my $super_loci_cutoff_str = ($super_loci_cutoff/1000)."k";
	        my $super_loci_max_file = $output_dir.$output_prefix."_fpkm$fpkm_cutoff"."_$super_loci_cutoff_str.superloci.max.bed";
	        my $super_loci_nearest_file = $output_dir.$output_prefix."_fpkm$fpkm_cutoff"."_$super_loci_cutoff_str.superloci.nearest.bed";
	        my @super_loci_list = &stitchSuperLoci(@beds);
		    &findNearestMaxGeneWithDistance($super_loci_nearest_file, $super_loci_max_file, 1, @super_loci_list);
		}
	}
}

sub findNearestGene {
    my ($multi_gene_file, $original_column_header, @beds) = @_;
	open(OUT_MULTI_GENE, ">$multi_gene_file") or die "Cannot open file $multi_gene_file";	
	my $header_str = "Chromosome\tStart\tEnd\t$original_column_header"."Up_gene\tUp_gene(dist)\tDn_gene\tDn_gene(dist)\tOverlap_gene\n";
	if ($show_exp) {
	    $header_str = "Chromosome\tStart\tEnd\t$original_column_header"."Up_gene\tUp_gene(dist)\tUp_gene(FPKM)\tDn_gene\tDn_gene(dist)\tDn_gene(FPKM)\tOverlap_gene\tOverlap_gene(FPKM)\n";
		if ($exp_type eq "DIFF") {
		    $header_str = "Chromosome\tStart\tEnd\t$original_column_header".
			    "Up_gene\tUp_gene(dist)\tUp_gene(sample1)\tUp_gene(sample2)\tUp_gene(value1)\tUp_gene(value2)\tUp_gene(logFC)\tUp_gene(pvalue)\tUp_gene(qvalue)\t".
				"Dn_gene\tDn_gene(dist)\tDn_gene(sample1)\tDn_gene(sample2)\tDn_gene(value1)\tDn_gene(value2)\tDn_gene(logFC)\tDn_gene(pvalue)\tDn_gene(qvalue)\t".
				"Overlap_gene\tOverlap_gene(sample1)\tOverlap_gene(sample2)\tOverlap_gene(value1)\tOverlap_gene(value2)\tOverlap_gene(logFC)\tOverlap_gene(pvalue)\tOverlap_gene(qvalue)\n";
		}
	}
	print OUT_MULTI_GENE $header_str;
	foreach my $locus (@beds) {
	     my $chr = $locus->chr;
	     my $start = $locus->start;
	     my $end = $locus->end;
		 my $original_column = $locus->original_column;
		 my %up_gene = ("id" => ".", "start" => 0, "end" => 0, "dist" => $max_dist, "exp_data" => 0); #id, start, end, distance, fpkm
		 my %dn_gene = ("id" => ".", "start" => 0, "end" => 0, "dist" => $max_dist, "exp_data" => 0); #id, start, end, distance, fpkm
		 my %overlap_gene = ("id" => ".", "start" => 0, "end" => 0, "exp_data" => 0); #id, start, end, distance, exp_data (either FPKM or log FC)
		 # for each gene in the same chromosome
		 my $total = 0;
		 my $not_found = 0;
		 foreach my $gene_id (keys %{$exp{$chr}}) {
			 my $tss_pos = $tss{$chr}{$gene_id}{"tss"};
			 my $g_start = $tss{$chr}{$gene_id}{"start"};
			 my $g_end = $tss{$chr}{$gene_id}{"end"};			
			 if (!$tss_pos) {
				next;				
			 }
			 my $exp_data = $exp{$chr}{$gene_id}{"fpkm"};
			 my $exp_list = "";
			 if (exists $exp{$chr}{$gene_id}{"exp_list"}) {
				$exp_list = "\t".$exp{$chr}{$gene_id}{"exp_list"};
			 }
			 if ($exp_type eq "DIFF") {
			     $exp_data = abs($exp{$chr}{$gene_id}{"lfc"});
			 }
			 
			 #upstream  
			 if ($start > $g_end) {
				  my $dist = $start - $tss_pos;
				  if ($dist < $up_gene{"dist"}) {
					   %up_gene = ("id" => $gene_id, "start" => $g_start, "end" => $g_end, "dist" => $dist, "exp_data" => $exp_data);
				  }
			 }                         
			 #dnstream  
			 if ($g_start > $end) {
				 my $dist = $tss_pos - $end;
				 if ($dist < $dn_gene{"dist"}) {
					   %dn_gene = ("id" => $gene_id, "start" => $g_start, "end" => $g_end, "dist" => $dist, "exp_data" => $exp_data);
				 }
			 }
			 #overlap
			 if (($g_start < $start && $g_end > $start) || ($g_start < $end && $g_end > $end) || ($g_start > $start && $g_end < $end)) {
				  if ($exp_data > $overlap_gene{"exp_data"}) {
					  %overlap_gene = ("id" => $gene_id, "start" => $g_start, "end" => $g_end, "exp_data" => $exp_data);
				  }
			 }
		 } #end of for each		 
		 my $up_dist = ($up_gene{"id"} eq ".")? "." : $up_gene{"dist"};
		 my $dn_dist = ($dn_gene{"id"} eq ".")? "." : $dn_gene{"dist"};
		 my $out_str = "$chr\t$start\t$end\t$original_column$up_gene{'id'}\t$up_dist\t$dn_gene{'id'}\t$dn_dist\t$overlap_gene{'id'}\n";
		 if ($show_exp) {		     		     
		 	 my $up_exp_data = ($up_gene{"id"} eq ".")? "." : $up_gene{"exp_data"};
		     my $dn_exp_data = ($dn_gene{"id"} eq ".")? "." : $dn_gene{"exp_data"};
		     my $overlap_exp_data = ($overlap_gene{"id"} eq ".")? "." : $overlap_gene{"exp_data"};
			 if ($exp_type eq "DIFF") {
			     $up_exp_data = &getCuffdiffColumns($chr, $up_gene{"id"});
				 $dn_exp_data = &getCuffdiffColumns($chr, $dn_gene{"id"});
				 $overlap_exp_data = &getCuffdiffColumns($chr, $overlap_gene{"id"});
			 }
			 $out_str = "$chr\t$start\t$end\t$original_column$up_gene{'id'}\t$up_dist\t$up_exp_data\t$dn_gene{'id'}\t$dn_dist\t$dn_exp_data\t$overlap_gene{'id'}\t$overlap_exp_data\n";
		 }
		 print OUT_MULTI_GENE $out_str;
	}
}

sub stitchSuperLoci {
    my (@beds) = @_;
    my @super_loci_list = ();
	my $super_loci;
	foreach my $locus (@beds) {
	     my $chr = $locus->chr;
		 my $start = $locus->start;
	     my $end = $locus->end;
		 #stitching
		 if ($super_loci) {
			      # if in the same TAD and within distance cutoff
				  if ($chr eq $super_loci->chr && ($start - $super_loci->end) <= $super_loci_cutoff) {
				       if ($tad_file && inSameTAD($super_loci->start, $super_loci->end, $start, $end, $tad{$chr})){
					       $super_loci->set_end($end);
					       $super_loci->set_loci_num($super_loci->loci_num + 1);
					   }
					   
				  }
				  # else new another superloci
				  else {
					   $super_loci = new Locus(chr => $chr, start => $start, end => $end, value => 0, loci_num => 1, original_column => "");
					   push @super_loci_list, $super_loci;
				  }
			  #initialize a superloci
		} else {
				  $super_loci = new Locus(chr => $chr, start => $start, end => $end, value => 0, loci_num => 1, original_column => "");
				  push @super_loci_list, $super_loci;
		} 		
	} # end of foreach
	return @super_loci_list;
}    

sub findNearestMaxGeneWithDistance {
    my ($nearest_file, $max_file, $is_stitched, $original_column_header, @input_bed) = @_;
	open(OUT_NEAREST_FILE, ">$nearest_file") or die "Cannot open file $nearest_file"; 
	open(OUT_MAX_FILE, ">$max_file") or die "Cannot open file $max_file"; 
	if (!$is_stitched) {
	     my $header_str = "Chromosome\tStart\tEnd\t$original_column_header"."Gene\tLength\n";
	     if ($show_exp) {
	         $header_str = "Chromosome\tStart\tEnd\t$original_column_header"."Gene\tFPKM\tLength$matrix_header\n";
		     if ($exp_type eq "DIFF") {
		        $header_str = "Chromosome\tStart\tEnd\t$original_column_header"."Gene\tSample1\tSample2\tValue1\tValue2\tlogFC\tpValue\tqValue\tLength\n";
		     }
	     }
		 print OUT_MAX_FILE $header_str;
		 print OUT_NEAREST_FILE $header_str;
	}
	foreach my $bed (@input_bed) {
	    my $chr = $bed->chr;
	    my $start = $bed->start;
	    my $end = $bed->end;
	    my $loci_num = $bed->loci_num;
		my $original_column = $bed->original_column;
	    my $max_exp_data = 0;
		my $max_exp_list = "";
	    my $max_gene = ".";
	    my $nearest_exp_data = 0;
		my $nearest_exp_list = "";
	    my $nearest_gene = ".";
	    my $nearest_dist = $max_dist;
	    foreach my $gene_id (keys %{$exp{$chr}}) {
			 #print "$row\n";
			 my $tss_pos = $tss{$chr}{$gene_id}{"tss"};
			 my $g_start = $tss{$chr}{$gene_id}{"start"};
			 my $g_end = $tss{$chr}{$gene_id}{"end"};
			 if (!$tss_pos) {
				next;				
			 }
			 my $exp_data = $exp{$chr}{$gene_id}{"fpkm"};
			 my $exp_list = "";
			 if (exists $exp{$chr}{$gene_id}{"exp_list"}) {
				$exp_list = "\t".$exp{$chr}{$gene_id}{"exp_list"};
			 }
             if ($exp_type eq "DIFF") {
			     $exp_data = abs($exp{$chr}{$gene_id}{"lfc"});
			 }			 
			 #upstream  
			 if ($start > $g_end) {
				  my $dist = $start - $tss_pos;
				  if ($dist < $nearest_gene_loci_cutoff) {
                      if (!$tad_file || ($tad_file && inSameTAD($start, $end, $g_start, $g_end, $tad{$chr}))) {
					      if ($dist < $nearest_dist) {
						      $nearest_dist = $dist;
						      $nearest_gene = $gene_id;
						      $nearest_exp_data = $exp_data;
							  $nearest_exp_list = $exp_list;
					      }
					      if ($exp_data > $max_exp_data) {
						      $max_exp_data = $exp_data;
						      $max_gene = $gene_id;
							  $max_exp_list = $exp_list;
					      }
					  }
				  }
			 }                         
			 #dnstream  
			 if ($g_start > $end) {
				  my $dist = $tss_pos - $end;
				  if ($dist < $nearest_gene_loci_cutoff) {
                      if (!$tad_file || ($tad_file && inSameTAD($start, $end, $g_start, $g_end, $tad{$chr}))) {
					      if ($dist < $nearest_dist) {
						      $nearest_dist = $dist;
						      $nearest_gene = $gene_id;
						      $nearest_exp_data = $exp_data;
							  $nearest_exp_list = $exp_list;
						  }					   
					      if ($exp_data > $max_exp_data) {
						      $max_exp_data = $exp_data;
							  $max_exp_list = $exp_list;
						      $max_gene = $gene_id;
					      }
					  }	   
				  }
			 }

			 #overlap
			 if (($g_start < $start && $g_end > $start) || ($g_start < $end && $g_end > $end) || ($g_start > $start && $g_end < $end)) {
				 if ($nearest_dist > 0) {
					 $nearest_dist = 0;
					 $nearest_gene = $gene_id;
					 $nearest_exp_data = $exp_data;
					 $nearest_exp_list = $exp_list;
				 } elsif ($exp_data > $nearest_exp_data) {             
					 $nearest_gene = $gene_id;
					 $nearest_exp_data = $exp_data;
					 $nearest_exp_list = $exp_list;
				  }
				  if ($exp_data > $max_exp_data) {
					 $max_exp_data = $exp_data;
					 $max_exp_list = $exp_list;
					 $max_gene = $gene_id;
				 }
			 }
	  } 
	  my $length = $end-$start;	  
	  my $max_out_str = "$chr\t$start\t$end\t$original_column$max_gene\t$length\n";
	  my $nearest_out_str = "$chr\t$start\t$end\t$original_column$nearest_gene\t$length\n";
	  if ($show_exp) {		     		     
		 	 my $max_exp_str = ($max_gene eq ".")? "." : $max_exp_data;
	         my $nearest_exp_str = ($nearest_gene eq ".")? "." : $nearest_exp_data;
		     if ($exp_type eq "DIFF") {
			     $max_exp_str = &getCuffdiffColumns($chr, $max_gene);
				 $nearest_exp_str = &getCuffdiffColumns($chr, $nearest_gene);
			 }
			 $max_out_str = "$chr\t$start\t$end\t$original_column$max_gene\t$max_exp_str\t$length$max_exp_list\n";
			 $nearest_out_str = "$chr\t$start\t$end\t$original_column$nearest_gene\t$nearest_exp_str\t$length$nearest_exp_list\n";
			 if ($is_stitched) {
			     $max_out_str = "$chr\t$start\t$end\t$original_column$max_gene\t$max_exp_str\t$loci_num$max_exp_list\n";
				 $nearest_out_str = "$chr\t$start\t$end\t$original_column$nearest_gene\t$nearest_exp_str\t$length\t$loci_num$nearest_exp_list\n";
			 }
	  }
	  print OUT_MAX_FILE $max_out_str;
	  print OUT_NEAREST_FILE $nearest_out_str;
	} # end of for each 	
}

sub readTAD {
  my ($tad_file) = @_;
  my %tad = ();
  open(FILE, $tad_file) or die "Cannot open file $tad_file";
  while (<FILE>) {
       chomp;
       my ($chr, $start, $end) = split(/\t/);
       push @{$tad{$chr}}, [$start,$end];

  }
  close(FILE);
  return %tad;
}

sub readCoordinate {
  my ($coor_file) = @_;
  open(FILE, $coor_file) or die "Cannot open file $coor_file";
  while (<FILE>) {
       chomp;
       my ($chr, $start, $end, $gene, $strand) = split(/\t/);
	   my $tss_pos = ($strand eq "+")? $start:$end;
       $tss{$chr}{$gene}{"tss"} = $tss_pos;
	   $tss{$chr}{$gene}{"start"} = $start;
	   $tss{$chr}{$gene}{"end"} = $end;
	   $coord{$gene} = "$chr:$start-$end";
	   #print $gene."\t".$tss{$chr}{$gene}{"tss"}."\t".$tss{$chr}{$gene}{"start"}."\t".$tss{$chr}{$gene}{"end"}."\n";
  }
  close(FILE);  
}

sub readBED {
  my ($bed_file) = @_;
  
  my @beds = ();
  open(INPUT_FILE, $bed_file) or die "Cannot open file $bed_file";
  while (<INPUT_FILE>) {
		 $_ =~ s/\r?\n$//;		 
		 if (/^#/) {
		     my @fields = split(/\t/, substr($_,1));
			 my $is_header = 1;
		     my $header = "";
		     foreach my $col (@original_columns) {
			    if (defined($fields[$col-1])) {
				    $header .= $fields[$col-1]."\t";
				}
				else {
				    $is_header = 0;
				    last;
				}
			 }
			 if ($is_header) {
	             $original_column_header = $header;
	         }
			 next;
		 }		 
		 my @fields = split(/\t/);
		 my $chr = $fields[0];
		 my $start = $fields[1];
		 my $end = $fields[2];
		 my $peak_id = $fields[3];
		 my $original_field = "";		 
		 foreach my $col (@original_columns) {
			if (($col - 1) > $#fields) {
				die "Column index out of bound!";
			}
		    $original_field .= $fields[$col-1]."\t";
		 }
		 my $locus = new Locus(chr => $chr, start => $start, end => $end, value => 0, loci_num => 1, original_column => $original_field);
		 push @beds, $locus;
  }
  close(INPUT_FILE);
  
  return @beds;     
}

sub getCuffdiffColumns {
  my ($chr, $gene_id) = @_;
  if ($gene_id eq ".") {
      return ".\t.\t.\t.\t.\t.\t.";
  }
  return $exp{$chr}{$gene_id}{"sample1"}."\t".$exp{$chr}{$gene_id}{"sample2"}."\t".$exp{$chr}{$gene_id}{"value1"}."\t".$exp{$chr}{$gene_id}{"value2"}."\t".$exp{$chr}{$gene_id}{"lfc"}."\t".$exp{$chr}{$gene_id}{"pvalue"}."\t".$exp{$chr}{$gene_id}{"qvalue"}
}

sub readExpression {
  my ($exp_file, $exp_type) = @_;
  if ($exp_type eq "TPM") {
      return readTPM($exp_file);
  }
  if ($exp_type eq "DIFF") {
      return readCuffdiff($exp_file);
  }
  if ($exp_type eq "MATRIX") {
      return readMatrix($exp_file);
  }
  my %exp = ();
  open(FILE, $exp_file) or die "Cannot open file $exp_file";
  my $header_str = <FILE>;
  chomp $header_str;
  my @headers = split(/\t/, $header_str);
  my $id_idx;
  my $tss_idx;
  my $locus_idx;
  my $fpkm_idx;
  for (my $i=0;$i<=$#headers;$i++) {
       $id_idx = $i if ($headers[$i] eq 'tracking_id');
       $tss_idx = $i if ($headers[$i] eq 'tss_id');
       $locus_idx = $i if ($headers[$i] eq 'locus');
       $fpkm_idx = $i if ($headers[$i] eq 'FPKM');
  }
  while (<FILE>) {
       chomp;
       my @fields = split(/\t/);
       my $gene_id = $fields[$id_idx];
       my $tss_id = $fields[$tss_idx];
       my $locus = $fields[$locus_idx];
       my $fpkm = $fields[$fpkm_idx];
       next if ($fpkm < $fpkm_cutoff);
       my ($chr, $start, $end) = $locus =~ /(.*):(\d+)\-(\d+)/;
       $chr = "chr".$chr;
       $exp{$chr}{$gene_id}{"start"} = $start;
       $exp{$chr}{$gene_id}{"end"} = $end;
       $exp{$chr}{$gene_id}{"fpkm"} = $fpkm;
  }
  close(FILE);
  return %exp;     
}

sub readTPM {
  my ($exp_file) = @_;
  my %exp = ();
  open(FILE, $exp_file) or die "Cannot open file $exp_file";
  my $header_str = <FILE>;
  while (<FILE>) {
       chomp;
       my ($chr, $start, $end, $gene_id, $tpm) = split(/\t/);	   
       next if ($tpm < $fpkm_cutoff);
       $chr = "chr".$chr;
       $exp{$chr}{$gene_id}{"start"} = $start;
       $exp{$chr}{$gene_id}{"end"} = $end;
       $exp{$chr}{$gene_id}{"fpkm"} = $tpm;
  }
  close(FILE);
  return %exp;     
}

sub readCuffdiff {
  my ($exp_file) = @_;
  my %exp = ();
  open(FILE, $exp_file) or die "Cannot open file $exp_file";
  my $header_str = <FILE>;
  chomp $header_str;
  my @headers = split(/\t/, $header_str);
  my $id_idx;
  my $sample1_idx;
  my $sample2_idx;
  my $value1_idx;
  my $value2_idx;
  my $lfc_idx;
  my $pvalue_idx;
  my $qvalue_idx;
  my $locus_idx;
  
  for (my $i=0;$i<=$#headers;$i++) {
       $id_idx = $i if ($headers[$i] eq 'gene');
       $locus_idx = $i if ($headers[$i] eq 'locus');
	   $sample1_idx = $i if ($headers[$i] eq 'sample_1');
	   $sample2_idx = $i if ($headers[$i] eq 'sample_2');
	   $value1_idx = $i if ($headers[$i] eq 'value_1');
	   $value2_idx = $i if ($headers[$i] eq 'value_2');
       $lfc_idx = $i if ($headers[$i] eq 'log2(fold_change)');
	   $pvalue_idx = $i if ($headers[$i] eq 'p_value');
	   $qvalue_idx = $i if ($headers[$i] eq 'q_value');
  } 
  while (<FILE>) {
       chomp;
       my @fields = split(/\t/);
       my $gene_id = $fields[$id_idx];       
       my $locus = $fields[$locus_idx];
       my $sample1 = $fields[$sample1_idx];
	   my $sample2 = $fields[$sample2_idx];
	   my $value1 = $fields[$value1_idx];
	   my $value2 = $fields[$value2_idx];
	   my $lfc = $fields[$lfc_idx];
	   my $pvalue = $fields[$pvalue_idx];
	   my $qvalue = $fields[$qvalue_idx];
       next if ($value1 < $fpkm_cutoff && $value2 < $fpkm_cutoff);
	   next unless ($lfc eq "inf" || $lfc eq "-inf" || abs($lfc) > $lfc_cutoff);
	   if ($pvalue_cutoff) {
			if ($pvalue > $pvalue_cutoff) {
				next;
			}
	   } elsif ($qvalue > $qvalue_cutoff) {
	       next;
	   }
       my ($chr, $start, $end) = $locus =~ /(.*):(\d+)\-(\d+)/;
       $chr = "chr".$chr;
       $exp{$chr}{$gene_id}{"start"} = $start;
       $exp{$chr}{$gene_id}{"end"} = $end;
       $exp{$chr}{$gene_id}{"sample1"} = $sample1;
	   $exp{$chr}{$gene_id}{"sample2"} = $sample2;
	   $exp{$chr}{$gene_id}{"value1"} = $value1;
	   $exp{$chr}{$gene_id}{"value2"} = $value2;
	   $exp{$chr}{$gene_id}{"lfc"} = $lfc;
	   $exp{$chr}{$gene_id}{"pvalue"} = $pvalue;
	   $exp{$chr}{$gene_id}{"qvalue"} = $qvalue;
	   #print $chr."\t".$gene_id."\t".$sample1."\t".$sample2."\t".$value1."\t".$sample2."\n";
  }
  close(FILE);
  return %exp;     
}

sub readMatrix {
  my ($exp_file) = @_;
  my %exp = ();
  open(FILE, $exp_file) or die "Cannot open file $exp_file";
  my $header_str = <FILE>;  
  chomp $header_str;
  my @headers = split(/\t/, $header_str);
  $matrix_header = join ("\t", splice(@headers,1));
  $matrix_header = "\t".$matrix_header;
  my $id_idx = 0;
  while (<FILE>) {
	chomp;
	my @fields = split(/\t/);
	my $g = $fields[$id_idx];
	my $gene_id = $g;
	if ($g =~ /(.*?)_/) {
		$gene_id = $1;	
	}
	my $fpkm = $fields[$matrix_col];
	next if ($fpkm < $fpkm_cutoff);
	if (exists $coord{$gene_id}) {
		my $locus = $coord{$gene_id};
		my ($chr, $start, $end) = $locus =~ /(.*):(\d+)\-(\d+)/;
		$exp{$chr}{$gene_id}{"start"} = $start;
		$exp{$chr}{$gene_id}{"end"} = $end;
		$exp{$chr}{$gene_id}{"fpkm"} = $fpkm;
		$exp{$chr}{$gene_id}{"exp_list"} = join("\t", splice(@fields,1));
	}
  }
  close(FILE);
  return %exp; 
}

sub inSameTAD {
    my ($g_start, $g_end, $e_start, $e_end, $tad_ref) = @_;
    foreach my $tad (@{$tad_ref}) {
        my $start = @{$tad}[0];
        my $end = @{$tad}[1];
        if ((($g_start < $start && $g_end > $start && (($g_end - $start)/($g_end - $g_start) > 0.5)) || ($g_start < $end && $g_end > $end && (($end - $g_start)/($g_end - $g_start) > 0.5)) || ($g_start > $start && $g_end < $end)) && 
            (($e_start < $start && $e_end > $start) || ($e_start < $end && $e_end > $end) || ($e_start > $start && $e_end < $end))){        
#if (@{$tad}[0] <= $g_start && @{$tad}[1] >= $g_end && @{$tad}[0] <= $e_start && @{$tad}[1] >= $e_end) { 
            return 1;
        }
    } 
    return 0;  
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

package Locus;

use Class::MethodMaker
    get_set       => [qw(chr start end value loci_num original_column)],
    new_hash_init => 'new';
    
