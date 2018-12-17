#!/usr/bin/perl -w
use strict;
use Getopt::Long qw(GetOptions);
use Time::Piece;
use DBI;
use File::Basename;

my $db_file;
my $exp_file;
my $bed_file;
my $super_file;
my $sample_id;
my $fpkm_cutoff = 5;

my $usage = <<__EOUSAGE__;

Usage:

$0 [options]

required options:

  -d    the SQLite database file
  -b    bedfile

optional options:

  -e    cufflink expression file
  -s    sample id
  -t    Super enhancer table file
  -f    FPKM cutoff (default: $fpkm_cutoff)

  
__EOUSAGE__



GetOptions (
  'd=s' => \$db_file,
  'b=s' => \$bed_file,
  'e=s' => \$exp_file,
  's=s' => \$sample_id,
  't=s' => \$super_file,
  'f=i' => \$fpkm_cutoff,
);


unless ($db_file && $bed_file && $out_file) {
    die $usage;
}

my $dbh = DBI->connect( "dbi:SQLite:$db_file" ) || die "Cannot connect: $DBI::errstr";
my $sth;

my %super = ();
if ($super_file) {
    $sth = $dbh->prepare("select * from Enhancer where sample_id = ?");
    $sth->execute($sample_id);
    my $found = $sth -> fetchrow_array();
    $sth->finish;
    open(SUPER_FILE, $super_file) or die "Cannot open file $super_file";
    if (!$found) {
        print "enhancer data not exists. inserting data to database...\n";
        $dbh->do("PRAGMA synchronous = OFF");
    }
    while (<SUPER_FILE>) {
         chomp;
         next if (/^#/);
         my @fields = split(/\t/);
         my $id = $fields[0];
         next if ($id eq "REGION_ID");
         my $chr = $fields[1];
         my $start = $fields[2];
         my $end = $fields[3];
         my $sample_signal = $fields[6];
         my $input_signal = $fields[7];
         my $rank = $fields[8];
         $super{$id} = [$sample_signal, $input_signal]; 
         if (!$found) {
             $sth = $dbh->prepare("insert into Enhancer values(?,?,?,?,?,?,?)");
             $sth->execute($sample_id, $chr, $start, $end, $sample_signal, $input_signal, $rank);  
         }
    }     
}

if ($exp_file) {
    if (!$sample_id) {
        $sample_id = basename($exp_file);
    }
    $sth = $dbh->prepare("select * from Expression where sample_id = ?");
    $sth->execute($sample_id);
    my $found = $sth -> fetchrow_array();
    $sth->finish;
    if (!$found) {
        print "expression data not exists. inserting data to database...\n";
        $dbh->do("PRAGMA synchronous = OFF");
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
             $sth = $dbh->prepare("insert into Expression values(?,?,?,?,?,?,?)");
             $sth->execute($sample_id, "chr$chr", $start, $end, $gene_id, $tss_id, $fpkm);  
        }
        close(FILE);   
        print "done!";
    }
}

my $out_regular = dirname($bed_file).basename($bed_file)."_fpkm$fpkm_cutoff.regular.bed";
my $out_super = dirname($bed_file).basename($bed_file)."_fpkm$fpkm_cutoff.super.bed";
open(BED_FILE, $bed_file) or die "Cannot open file $bed_file";
open(OUT_REGULAR, ">$out_regular") or die "Cannot open file $out_regular";
open(OUT_SUPER, ">$out_super") or die "Cannot open file $out_super";

my %super_loci;
my $current_chr = "";
my $resultset_ref;
while (<BED_FILE>) {
     chomp;
     next if (/^(track|#)/);
     my ($chr, $start, $end, $peak_id, @others) = split(/\t/);     
     if ($current_chr ne $chr) {
         my $sql = "select g.gene_id,g.start,g.end,g.strand,e.fpkm from Gene_Coordinates g, Expression e where g.chr=? and e.chr = ? and g.gene_id=e.gene_id and e.fpkm >= ? and sample_id = ? order by g.start";
         $sth = $dbh->prepare($sql);
         $sth->execute($chr,$chr,$fpkm_cutoff,$sample_id);
         $resultset_ref = $sth->fetchall_arrayref({});
         $sth->finish;
         $current_chr = $chr
     }   
     my $up_dist = 9999999999;
     my $up_fpkm = 0;
     my $up_gene = "";
     my $dn_dist = 9999999999;
     my $dn_fpkm = 0;
     my $dn_gene = "";
     my $overlap_fpkm = 0;
     my $overlap_gene = ".";
     foreach my $row (@{$resultset_ref}) {
         #print "$row\n";
         my $gene_id = $row->{gene_id};
         my $g_start = $row->{start};
         my $g_end = $row->{end};
         my $strand = $row->{strand};
         my $fpkm = $row->{fpkm};
         #print "$gene_id,$g_start, $g_end, $strand, $fpkm\n";
         #upstream  
         if ($start > $g_end) {
              my $dist = $start - $g_end;
              if ($dist < $up_dist) {
                   $up_dist = $dist;
                   $up_gene = $gene_id;
                   $up_fpkm = $fpkm;
              }
         }                         
         #dnstream  
         if ($g_start > $end) {
             my $dist = $g_start - $end;
             if ($dist < $dn_dist) {
                 $dn_dist = $dist;
                 $dn_gene = $gene_id;
                 $dn_fpkm = $fpkm;
             }
         }
         #overlap
         if (($g_start < $start && $g_end > $start) || ($g_start < $end && $g_end > $end) || ($g_start > $start && $g_end < $end)) {
              if ($fpkm > $overlap_fpkm) {
                  $overlap_fpkm = $fpkm;
                  $overlap_gene = $gene_id;
              }
         }
     } #end of for each
     my $overlap_fpkm_str = ($overlap_fpkm == 0)? "." : $overlap_fpkm;
     my $signal_str = ".\t.";
     if ($super{$peak_id}) {
         $signal_str = $super{$peak_id}[0]."\t".$super{$peak_id}[1];
     }
     print OUT_FILE "$chr\t$start\t$end\t$peak_id\t$up_gene\t$up_dist\t$up_fpkm\t$dn_gene\t$dn_dist\t$dn_fpkm\t$overlap_gene\t$overlap_fpkm_str\t$signal_str\n";
} # end of while

=old
while (<BED_FILE>) {
     chomp;
     next if (/^(track|#)/);
     my ($chr, $start, $end, $peak_id, @others) = split(/\t/);
     #find upstream
     my $sql = "select g.*,e.fpkm from (select *,start-? as dist from Gene_Coordinates where start > ? and chr=? and strand = ? union select *,?-end as dist from Gene_Coordinates where ? > end and chr=? and strand = ? ) g, Expression e where g.gene_id=e.gene_id and e.chr=? and e.fpkm > ? and sample_id = ? order by dist limit 1";
     $sth = $dbh->prepare($sql);
     $sth->execute($end,$end,$chr,'-',$start,$start,$chr,'+',$chr,$fpkm_cutoff,$sample_id);
     my $upstream_result = "";
     if (my @row = $sth -> fetchrow_array()) {
         $upstream_result = join("\t",splice(@row, 1, 6));
     }
     $sth->finish;
     #find downstream
     my $downstream_result = "";
     $sth->execute($end,$end,$chr,'+',$start,$start,$chr,'-',$chr,$fpkm_cutoff,$sample_id);
     if (my @row = $sth -> fetchrow_array()) {
         $downstream_result = join("\t",splice(@row, 1, 6));
     }
     $sth->finish;
     #find intronic/overlapping
     $sql = "select g.*,e.fpkm from Gene_Coordinates g, Expression e where g.chr=? and e.chr=? and g.gene_id=e.gene_id and ((g.start < ? and g.end > ?) or (g.start < ? and g.end > ?) or (g.start > ? and g.end < ?)) and e.fpkm > ? and sample_id = ? order by fpkm limit 1";
     $sth = $dbh->prepare($sql);
     my $overlap_result = ".\t.\t.\t.\t.";
     $sth->execute($chr,$chr,$start,$start,$end,$end,$start,$end,$fpkm_cutoff,$sample_id);
     if (my @row = $sth -> fetchrow_array()) {
         $overlap_result = join("\t",splice(@row, 1, 5));
     }
     $sth->finish;

     print OUT_FILE "$chr\t$start\t$end\t$peak_id\t$upstream_result\t$downstream_result\t$overlap_result\n";
 
}
=cut
close(OUT_FILE);
close(BED_FILE);
$dbh->disconnect;

sub inTAD {
    my ($chr, $start, $end) = @_;
    my $sql = "select start,end from TAD where chr=? and start >= ? and end <= ?";
    my $sth = $dbh->prepare($sql);
    $sth->execute($chr,$start, $end);
    return $sth->fetchrow_array();
}

