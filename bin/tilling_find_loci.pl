#!/usr/bin/env perl 
#===============================================================================
#
#         FILE:  tilling_find_loci.pl
#
#        USAGE:  ./tilling_find_loci.pl  
#
#  DESCRIPTION:  Script to read input vcf file and identify clusters of SNP's
#                   which may correspond to loci important for a 
#                   specific phenotype.
#
#      OPTIONS:  ---
# REQUIREMENTS:  ---
#         BUGS:  ---
#        NOTES:  ---
#       AUTHOR:  Dr. Scott A. Givan (sag), givans@missouri.edu
#      COMPANY:  University of Missouri, USA
#      VERSION:  1.0
#      CREATED:  11/30/15 16:22:21
#     REVISION:  ---
#===============================================================================

use 5.010;      # Require at least Perl version 5.10
use autodie;
use Getopt::Long; # use GetOptions function to for CL args
use Statistics::Descriptive;
use warnings;
use strict;

my ($debug,$verbose,$help,$infile,$minlocs);

my $result = GetOptions(
    "debug"     =>  \$debug,
    "verbose"   =>  \$verbose,
    "help"      =>  \$help,
    "infile:s"      =>  \$infile,
    "minlocs:i"     =>  \$minlocs,
);

$minlocs = 10 unless ($minlocs);

if ($help) {
    help();
    exit(0);
}

open (my $IN,"<",$infile);

my $stat = Statistics::Descriptive::Full->new();

say "ID\tvariance\tkurtosis\tCv";

my ($init,$lastID) = (0,'noname');
while (<$IN>) {
    next if ($_ =~ /^#/);
    my @vals = split /\t/, $_;

    if ($init) {
        $stat->clear();
        $stat->add_data($vals[1]);
    } elsif ($vals[0] eq $lastID) {
        $stat->add_data($vals[1]);
    } else {
        say "$lastID\t" . $stat->variance . "\t" . $stat->kurtosis() . "\t" . $stat->standard_deviation()/$stat->mean() if ($stat->count() >= $minlocs);
        $stat->clear();
        $stat->add_data($vals[1]);
    }


    $lastID = $vals[0];
}

say "$lastID\t" . $stat->variance . "\t" . $stat->kurtosis() . "\t" . $stat->standard_deviation()/$stat->mean() if ($stat->count() >= $minlocs);

sub help {

say <<HELP;

    "debug"     =>  \$debug,
    "verbose"   =>  \$verbose,
    "help"      =>  \$help,
    "infile:s"      =>  \$infile,
    "minlocs:i"     =>  \$minlocs,

HELP

}



