#!/usr/bin/perl -w
use strict;
use Getopt::Long qw(GetOptions);
use Time::Piece;
use DBI;
use File::Basename;

my $tad_file;
my $exp_file;
my $enhancer_file;
my $fpkm_cutoff = 5;
my $super_loci_cutoff = 150000;
my $nearest_gene_loci_cutoff = 500000;
my $max_dist = 249250621;              #this is the max chromosome length
my $usage = <<__EOUSAGE__;

Usage:

$0 [options]

required options:

  -e    Cufflink expression file
  -t    Enhancer table file (ROSE format)
  -d    TAD file

optional options:

  -f    FPKM cutoff (default: $fpkm_cutoff)
  -s    Super loci length cutoff (default: $super_loci_cutoff)
  -n    The distance of nearest gene to the super loci (default: $nearest_gene_loci_cutoff)

  
__EOUSAGE__



GetOptions (
  'e=s' => \$exp_file,
  'd=s' => \$tad_file,
  't=s' => \$enhancer_file,
  'f=i' => \$fpkm_cutoff,
  's=i' => \$super_loci_cutoff,
);


unless ($exp_file && $enhancer_file && $tad_file) {
    die $usage;
}


my %exp = ();
my %tad = ();

%exp = &readExpression($exp_file);
%tad = &readTAD($tad_file);

my $super_loci_cutoff_str = ($super_loci_cutoff/1000)."k";
my %super = ();
my $output_dir = dirname($enhancer_file)."/";
my $cmd_log_file = $output_dir."/log/cmd_".localtime->strftime('%Y-%m-%d_%H:%M:%S').".log";
my $out_regular = $output_dir.basename($enhancer_file, ".txt")."_fpkm$fpkm_cutoff.regular.bed";
my $out_super = $output_dir.basename($enhancer_file, ".txt")."_fpkm$fpkm_cutoff.super.bed";
my $out_super_loci = $output_dir.basename($enhancer_file, ".txt")."_fpkm$fpkm_cutoff.$super_loci_cutoff_str.superloci.max.bed";
my $out_super_loci_nearest = $output_dir.basename($enhancer_file, ".txt")."_fpkm$fpkm_cutoff.$super_loci_cutoff_str.superloci.nearest.bed";
open(OUT_REGULAR, ">$out_regular") or die "Cannot open file $out_regular";
open(OUT_SUPER, ">$out_super") or die "Cannot open file $out_super";
open(ENHANCER_FILE, $enhancer_file) or die "Cannot open file $enhancer_file";
#&print_log("output reqular:".$out_regular);
#&print_log("output super:".$out_super);
my @super_loci_list = ();
my $super_loci;
my $current_chr = "";
my $resultset_ref;
my $tad_ref;
my $chip_sample = basename($enhancer_file);
# the ROSE table file must be sorted
# use command: sort -k2,2 -k3n sample_AllEnhancers.table.txt > sample_AllEnhancers.sorted.table.txt
while (<ENHANCER_FILE>) {
     chomp;
     next if (/^#/);
     my @fields = split(/\t/);
     my $peak_id = $fields[0];
     next if ($peak_id eq "REGION_ID");
     my $chr = $fields[1];
     my $start = $fields[2];
     my $end = $fields[3];
     my $sample_signal = $fields[6];
     my $input_signal = $fields[7];
     my $rank = $fields[8];
     my $isuper = $fields[9];
     
     my %up_gene = ("id" => ".", "start" => 0, "end" => 0, "dist" => $max_dist, "fpkm" => 0); #id, start, end, distance, fpkm
     my %dn_gene = ("id" => ".", "start" => 0, "end" => 0, "dist" => $max_dist, "fpkm" => 0); #id, start, end, distance, fpkm
     my %overlap_gene = ("id" => ".", "start" => 0, "end" => 0, "fpkm" => 0); #id, start, end, distance, fpkm
     #&print_log("Processing peak: $peak_id");
     foreach my $gene_id (keys %{$exp{$chr}}) {
         #print "$row\n";
         my $g_start = $exp{$chr}{$gene_id}{"start"};
         my $g_end = $exp{$chr}{$gene_id}{"end"};
         my $fpkm = $exp{$chr}{$gene_id}{"fpkm"};
         #upstream  
         if ($start > $g_end) {
              my $dist = $start - $g_end;
              if ($dist < $up_gene{"dist"}) {
                   %up_gene = ("id" => $gene_id, "start" => $g_start, "end" => $g_end, "dist" => $dist, "fpkm" => $fpkm);
              }
         }                         
         #dnstream  
         if ($g_start > $end) {
             my $dist = $g_start - $end;
             if ($dist < $dn_gene{"dist"}) {
                   %dn_gene = ("id" => $gene_id, "start" => $g_start, "end" => $g_end, "dist" => $dist, "fpkm" => $fpkm);
             }
         }
         #overlap
         if (($g_start < $start && $g_end > $start) || ($g_start < $end && $g_end > $end) || ($g_start > $start && $g_end < $end)) {
              if ($fpkm > $overlap_gene{"fpkm"}) {
                  %overlap_gene = ("id" => $gene_id, "start" => $g_start, "end" => $g_end, "fpkm" => $fpkm);
              }
         }
     } #end of for each
     my $up_fpkm_str = ($up_gene{"fpkm"} == 0)? "." : $up_gene{"fpkm"};
     my $dn_fpkm_str = ($dn_gene{"fpkm"} == 0)? "." : $dn_gene{"fpkm"};
     my $overlap_fpkm_str = ($overlap_gene{"fpkm"} == 0)? "." : $overlap_gene{"fpkm"};
     my $out_str = "$chr\t$start\t$end\t$peak_id\t$up_gene{'id'}\t$up_gene{'dist'}\t$up_fpkm_str\t$dn_gene{'id'}\t$dn_gene{'dist'}\t$dn_fpkm_str\t$overlap_gene{'id'}\t$overlap_fpkm_str\t$sample_signal\t$input_signal\n";
     if ($isuper) {
          print OUT_SUPER $out_str;
          if ($super_loci) {
              if (($start - $super_loci->end) <= $super_loci_cutoff && inSameTAD($super_loci->start, $super_loci->end, $start, $end, $tad{$chr})) {
                   $super_loci->set_end($end);
                   $super_loci->set_signal($super_loci->signal + $sample_signal);
              }
              else {
                   $super_loci = new Super_loci(chr => $chr, start => $start, end => $end, signal => $sample_signal);
                   push @super_loci_list, $super_loci;
              }
          } else {
              $super_loci = new Super_loci(chr => $chr, start => $start, end => $end, signal => $sample_signal);
              push @super_loci_list, $super_loci;
          } 
          #&findSuperLoci($up_gene{"id"}, $chr, $up_gene{"start"}, $up_gene{"end"}, $start, $end, $up_gene{"fpkm"}, $sample_signal, $tad{$chr});
          #&findSuperLoci($dn_gene{"id"}, $chr, $dn_gene{"start"}, $dn_gene{"end"}, $start, $end, $dn_gene{"fpkm"}, $sample_signal, $tad{$chr});
          #&findSuperLoci($overlap_gene{"id"}, $chr, $overlap_gene{"start"}, $overlap_gene{"end"}, $start, $end, $overlap_gene{"fpkm"}, $sample_signal, $tad{$chr});

     }
     else {
          print OUT_REGULAR $out_str;
     }
     
} # end of while

close(OUT_SUPER);
close(OUT_REGULAR);
close(ENHANCER_FILE);

#find genes associated with super loci
open(OUT_SUPER_LOCI, ">$out_super_loci") or die "Cannot open file $out_super_loci";
open(OUT_SUPER_LOCI_NEAR, ">$out_super_loci_nearest") or die "Cannot open file $out_super_loci_nearest";
foreach my $super_loci (@super_loci_list) {
  my $chr = $super_loci->chr;
  my $start = $super_loci->start;
  my $end = $super_loci->end;
  my $signal = $super_loci->signal;
  my $max_fpkm = 0;
  my $max_gene = ".";
  my $nearest_fpkm = 0;
  my $nearest_gene = ".";
  my $nearest_dist = $max_dist;
  #&print_log("$chr $start $end");
  foreach my $gene_id (keys %{$exp{$chr}}) {
         #print "$row\n";
         my $g_start = $exp{$chr}{$gene_id}{"start"};
         my $g_end = $exp{$chr}{$gene_id}{"end"};
         my $fpkm = $exp{$chr}{$gene_id}{"fpkm"};  
=old       
         #upstream  
         if ($start > $g_end) {
              my $dist = $start - $g_end;
              if ($dist < $up_min_dist && inSameTAD($start, $end, $g_start, $g_end, $tad{$chr})) {
                   if ($fpkm > $max_fpkm) {
                       $max_fpkm = $fpkm;
                       $max_gene = $gene_id;
                   }
              }
         }                         
         #dnstream  
         if ($g_start > $end) {
             my $dist = $g_start - $end;
             if ($dist < $dn_min_dist && inSameTAD($start, $end, $g_start, $g_end, $tad{$chr})) {
                   if ($fpkm > $max_fpkm) {
                       $max_fpkm = $fpkm;
                       $max_gene = $gene_id;
                   }
             }
         }
=cut
         #upstream  
         if ($start > $g_end) {
              my $dist = $start - $g_end;
              if ($dist < $nearest_gene_loci_cutoff && inSameTAD($start, $end, $g_start, $g_end, $tad{$chr})) {
                   if ($dist < $nearest_dist) {
                       $nearest_dist = $dist;
                       $nearest_gene = $gene_id;
                       $nearest_fpkm = $fpkm;
                   }
                   if ($fpkm > $max_fpkm) {
                       $max_fpkm = $fpkm;
                       $max_gene = $gene_id;
                   }
              }
         }                         
         #dnstream  
         if ($g_start > $end) {
             my $dist = $g_start - $end;
              if ($dist < $nearest_gene_loci_cutoff && inSameTAD($start, $end, $g_start, $g_end, $tad{$chr})) {
                   if ($dist < $nearest_dist) {
                       $nearest_dist = $dist;
                       $nearest_gene = $gene_id;
                       $nearest_fpkm = $fpkm;
                   }
                   if ($fpkm > $max_fpkm) {
                       $max_fpkm = $fpkm;
                       $max_gene = $gene_id;
                   }
              }
         }

         #overlap
         if (($g_start < $start && $g_end > $start) || ($g_start < $end && $g_end > $end) || ($g_start > $start && $g_end < $end)) {
             if ($nearest_dist > 0) {
                 $nearest_dist = 0;
                 $nearest_gene = $gene_id;
                 $nearest_fpkm = $fpkm;
             } elsif ($fpkm > $nearest_fpkm) {             
                 $nearest_gene = $gene_id;
                 $nearest_fpkm = $fpkm;
              }
              if ($fpkm > $max_fpkm) {
                 $max_fpkm = $fpkm;
                 $max_gene = $gene_id;
             }
         }
  }     
  #print OUT_SUPER_LOCI $super_loci{$gene_id}{"chr"}."\t".$super_loci{$gene_id}{"start"}."\t".$super_loci{$gene_id}{"end"}."\t$gene_id\t".($super_loci{$gene_id}{"end"} - $super_loci{$gene_id}{"start"})."\t".$super_loci{$gene_id}{"signal"}."\t$fpkm_str\n";
  my $max_fpkm_str = ($max_fpkm == 0)? "." : $max_fpkm;
  my $nearest_fpkm_str = ($nearest_fpkm == 0)? "." : $nearest_fpkm;
  print OUT_SUPER_LOCI "$chr\t$start\t$end\t$max_gene\t".($end-$start)."\t$max_fpkm_str\t$signal\n";
  print OUT_SUPER_LOCI_NEAR "$chr\t$start\t$end\t$nearest_gene\t".($end-$start)."\t$nearest_fpkm_str\t$signal\n";
}
close(OUT_SUPER_LOCI);
close(OUT_SUPER_LOCI_NEAR);

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

sub readExpression {
  my ($exp_file) = @_;
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

=old
sub findSuperLoci {
    my ($gene, $chr, $g_start, $g_end, $e_start, $e_end, $fpkm, $sample_signal, $tad_ref) = @_;
    return if ($gene eq ".");
    if (inSameTAD($g_start, $g_end, $e_start, $e_end, $tad_ref)){            
            if ($super_loci{$gene}) {
                if ($super_loci{$gene}{"start"} > $e_start) {
                    $super_loci{$gene}{"start"} = $e_start;
                }
                if ($super_loci{$gene}{"end"} < $e_end) {
                    $super_loci{$gene}{"end"} = $e_end;
                }
                $super_loci{$gene}{"signal"} += $sample_signal;
            }
            else {
                $super_loci{$gene}{"chr"} = $chr;
                $super_loci{$gene}{"start"} = $e_start;
                $super_loci{$gene}{"end"} = $e_end; 
                $super_loci{$gene}{"signal"} = $sample_signal;
                $super_loci{$gene}{"fpkm"} = $fpkm ;
            }
    }
}
=cut
sub inSameTAD {
    my ($g_start, $g_end, $e_start, $e_end, $tad_ref) = @_;
    foreach my $tad (@{$tad_ref}) {
        #&print_log("@{$tad}[0] @{$tad}[0] $g_start $g_end $e_start $e_end");
        my $start = @{$tad}[0];
        my $end = @{$tad}[1];
        if ((($g_start < $start && $g_end > $start) || ($g_start < $end && $g_end > $end) || ($g_start > $start && $g_end < $end)) && 
            (($e_start < $start && $e_end > $start) || ($e_start < $end && $e_end > $end) || ($e_start > $start && $e_end < $end))){        
#if (@{$tad}[0] <= $g_start && @{$tad}[1] >= $g_end && @{$tad}[0] <= $e_start && @{$tad}[1] >= $e_end) { 
            return 1;
        }
    } 
    return 0;  
}

sub print_log {
    my ($msg) = @_;
    open CMD_FILE, ">>$cmd_log_file" || print "cannot create command log file";
    print CMD_FILE "[".localtime->strftime('%Y-%m-%d %H:%M:%S')."] $msg\n";
    close(CMD_FILE);
}

sub getPathNoExtension {
    my ($path, $ext) = @_;
    return dirname($enhancer_file)."/".basename($enhancer_file, ".txt");
}

package Super_loci;

use Class::MethodMaker
    get_set       => [qw(chr start end signal)],
    new_hash_init => 'new';
    
