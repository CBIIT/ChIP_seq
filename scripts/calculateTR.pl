#!/usr/bin/perl -w
use strict;
use Getopt::Long qw(GetOptions);
use Time::Piece;
use DBI;
use File::Basename;
use Cwd;

my $bam_file;
my $out_file;
my $ref_file = "/data/khanlab/projects/ChIP_seq/data_by_file_type/test/data_by_file_type/ref/gene_coord.tsv";
my $pro_start = -800;
my $pro_end = -30;
my $tssr_start = -30;
my $tssr_end = 300;
my $gene_start = 300;
my $gene_end = 0;
my $tesr_start = 0;
my $tesr_end = 4000;
my $exp_file;
my $usage = <<__EOUSAGE__;

Usage:

$0 [options]

required options:

  -b  <string>  input BAM file
  -o  <string>  Output file
  
optional options:  
  -r  <string>  reference file: (default: $ref_file)
  -e  <string>  expression file (Khanlab TPM format): (default: None)
  -i  <integer> Promoter start position around TSS(default: $pro_start)
  -j  <integer> Promoter end position around TSS(default: $pro_end)
  -p  <integer> TSSR start position around TSS(default: $tssr_start)
  -q  <integer> TSSR end position around TSS(default: $tssr_end)  
  -u  <integer> Gene body start position around TSS(default: $gene_start)
  -v  <integer> Gene body end position around TES(default: $gene_end)  
  -x  <integer> TESR start position around TES(default: $tesr_start)
  -y  <integer> TESR end position around TES(default: $tesr_end)  
  
__EOUSAGE__



GetOptions (
  'b=s' => \$bam_file,
  'o=s' => \$out_file,
  'e=s' => \$exp_file,
  'r=s' => \$ref_file,
  'i=i' => \$pro_start,
  'j=i' => \$pro_end,
  'p=i' => \$tssr_start,
  'q=i' => \$tssr_end,
  'u=i' => \$gene_start,
  'v=i' => \$gene_end,
  'x=i' => \$tesr_start,
  'y=i' => \$tesr_end  
);

unless ($bam_file && $out_file) {
    die "Please input BAM and output TRV file\n$usage";
}

&main();
my %exp = ();

sub main {
	my $bed_file =  basename($bam_file).".bed";
	my $cov_file = basename($bam_file).".cov";
	&makeBEDFile($ref_file, $bed_file, $pro_start, $pro_end, $tssr_start, $tssr_end, $gene_start, $gene_end, $tesr_start, $tesr_end);
	if (! -e $cov_file) {
		runCommand("bedtools multicov -bams $bam_file -bed $bed_file > $cov_file");	
	}
	if ($exp_file) {
		%exp = &readExpression($exp_file);
	}
	&makeTRVFile($cov_file, $out_file);
	
}

sub makeBEDFile {
  my ($ref_file, $bed_file, $pro_start, $pro_end, $tssr_start, $tssr_end, $gene_start, $gene_end, $tesr_start, $tesr_end) = @_;
  open(REF_FILE, $ref_file) or die "Cannot open file $ref_file";
  open(BED_FILE, ">$bed_file") or die "Cannot open file $bed_file";
  <REF_FILE>;
  while (<REF_FILE>) {
       chomp;
       my ($chr, $start, $end, $strand, $gene, $overlap) = split(/\t/);
	   my $pro_start_pos = ($strand eq "+")? $start + $pro_start : $end - $pro_end;
	   my $pro_end_pos = ($strand eq "+")? $start + $pro_end : $end - $pro_start;
	   my $tssr_start_pos = ($strand eq "+")? $start + $tssr_start : $end - $tssr_end;
	   my $tssr_end_pos = ($strand eq "+")? $start + $tssr_end : $end - $tssr_start;
	   my $gene_start_pos = ($strand eq "+")? $start + $gene_start : $start - $gene_end;
	   my $gene_end_pos = ($strand eq "+")? $end + $gene_end : $end - $gene_start;
	   my $tesr_start_pos = ($strand eq "+")? $end + $tesr_start : $start - $tesr_end;
	   my $tesr_end_pos = ($strand eq "+")? $end + $tesr_end : $start - $tesr_start;
	   next if ($gene_end_pos <= $gene_start_pos);
       print BED_FILE "$chr\t$pro_start_pos\t$pro_end_pos\t$strand\t$gene\tPromoter\t$overlap\n";	   
	   print BED_FILE "$chr\t$tssr_start_pos\t$tssr_end_pos\t$strand\t$gene\tTSSR\t$overlap\n";
	   print BED_FILE "$chr\t$gene_start_pos\t$gene_end_pos\t$strand\t$gene\tGene Body\t$overlap\n";
	   print BED_FILE "$chr\t$tesr_start_pos\t$tesr_end_pos\t$strand\t$gene\tTESR\t$overlap\n";
  }
  close(REF_FILE);
  close(BED_FILE);  
}

sub makeTRVFile {
  my ($cov_file, $trv_file) = @_;
  open(COV_FILE, $cov_file) or die "Cannot open file $cov_file";
  open(TRV_FILE, ">$trv_file") or die "Cannot open file $trv_file";
  while (<COV_FILE>) {
       chomp;
	   my ($chr_pro, $start_pro, $end_pro, $strand_pro, $gene_pro, $type_pro, $overlap_pro, $pro_value) = split(/\t/);
	   my $line_tssr = <COV_FILE>;
	   chomp $line_tssr;
       my ($chr_tssr, $start_tssr, $end_tssr, $strand_tssr, $gene_tssr, $type_tssr, $overlap_tssr, $tssr_value) = split(/\t/, $line_tssr);
	   my $line_gene = <COV_FILE>;
	   chomp $line_gene;
       my ($chr_gene, $start_gene, $end_gene, $strand_gene, $gene_gene, $type_gene, $overlap_gene, $gene_value) = split(/\t/, $line_gene);
	   my $line_tesr = <COV_FILE>;
	   chomp $line_tesr;
       my ($chr_tesr, $start_tesr, $end_tesr, $strand_tesr, $gene_tesr, $type_tesr, $overlap_tesr, $tesr_value) = split(/\t/, $line_tesr);	   
	   #print "$gene1\t$gene2\t$tssr_value\t$gene_body_value\n";
	   my $ratio = ($gene_value == 0)? 0 : ($tssr_value/($end_tssr-$start_tssr+1))/($gene_value/($end_gene-$start_gene+1));
	   print TRV_FILE "$chr_tssr\t$start_tssr\t$end_tssr\t$start_gene\t$end_gene\t$strand_tssr\t$gene_tssr\t$overlap_tssr\t$ratio\t$pro_value\t$tssr_value\t$gene_value\t$tesr_value";
	   if ($exp_file) {
			my $exp_value = "NA";
			if (exists $exp{$chr_tssr}{$gene_tssr}) {
				$exp_value = $exp{$chr_tssr}{$gene_tssr};				
			}
			print TRV_FILE "\t$exp_value\n";
	   } else {
			print TRV_FILE "\n";
	   }
  }
  close(COV_FILE);
  close(TRV_FILE);  
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

sub readExpression {
  my ($exp_file) = @_;
  my %exp = ();
  open(FILE, $exp_file) or die "Cannot open file $exp_file";
  my $header_str = <FILE>;
  chomp $header_str;
  my @headers = split(/\t/, $header_str);
  while (<FILE>) {
       chomp;
       my @fields = split(/\t/);
       my $chr = $fields[0];
	   my $gene_id = $fields[3];
	   my $tpm = $fields[4];
	   $chr = "chr".$chr;
       $exp{$chr}{$gene_id} = $tpm;
  }
  close(FILE);
  return %exp;     
}

sub readCufflinks {
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
       my ($chr, $start, $end) = $locus =~ /(.*):(\d+)\-(\d+)/;
       $chr = "chr".$chr;
       $exp{$chr}{$gene_id} = $fpkm;
  }
  close(FILE);
  return %exp;     
}
