#!/usr/bin/perl -w
use strict;
use Getopt::Long qw(GetOptions);
use Time::Piece;
use DBI;
use File::Basename;
use Config::Simple;

my $config_file = $ARGV[0];

my $config = new Config::Simple($config_file);

my @cutoff = $config->param("MACS.pvalue_cutoff");

foreach my $cutoff (@cutoff) {
   print $cutoff."\n";
}