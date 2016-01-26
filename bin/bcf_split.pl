#!/usr/bin/env perl 
#===============================================================================
#
#         FILE:  bcf_split.pl
#
#        USAGE:  ./bcf_split.pl  
#
#  DESCRIPTION:  Script to create multiple output bcf files for regions of
#                   input bcf files.
#
#      OPTIONS:  ---
# REQUIREMENTS:  ---
#         BUGS:  ---
#        NOTES:  ---
#       AUTHOR:  Dr. Scott A. Givan (sag), givans@missouri.edu
#      COMPANY:  University of Missouri, USA
#      VERSION:  1.0
#      CREATED:  01/26/16 12:27:24
#     REVISION:  ---
#===============================================================================

use 5.010;      # Require at least Perl version 5.10
use autodie;
use Getopt::Long; # use GetOptions function to for CL args
use warnings;
use strict;

my ($debug,$verbose,$help,$infile,$regionsize,$refseqname);

my $result = GetOptions(
    "debug"     =>  \$debug,
    "verbose"   =>  \$verbose,
    "infile:s"      =>  \$infile,
    "regionsize:i"  =>  \$regionsize,
    "refseqname:s"  =>  \$refseqname,
    "help"      =>  \$help,
);

if ($help) {
    help();
    exit(0);
}

my $bcftools = 'bcftools';
my $tabix = 'tabix';
$regionsize = 1000000 unless ($regionsize);
$infile = 'infile' unless ($infile);
$refseqname = 'refseq' unless ($refseqname);

my $outfile = $infile;
$outfile =~ s/\.gz//;
$outfile =~ s/\.[bv]cf//;
say "outfile stub:  '$outfile'" if ($debug);

my $cnt = 0;
#my ($j,$k) = (0,0)
for (my ($j,$k,$i) = (1,$regionsize,0); $i != -1; $j += $regionsize, $k += $regionsize) {


    my $regionoutfile = "$outfile.$j-$k.bcf.gz";
    #say "outfile $cnt: '$regionoutfile'" if ($debug || $verbose);

    open(BCF,"|-","$bcftools view -O b -o $regionoutfile $infile $refseqname:$j-$k");
    close(BCF);

    open(FTEST,"-|","bcftools view -O v $regionoutfile | tail -n 1");
    chomp(my $ftest = <FTEST>);
    close(FTEST);
    say "\$ftest = '$ftest'" if ($debug);
    $i = -1 if (substr($ftest,0,1) eq '#');# exit loop when you reach the end of the reference sequence
    unlink($regionoutfile) if ($i == -1);
    say "outfile $cnt: '$regionoutfile'" if (($debug || $verbose) && $i != -1);

    unless ($i == -1) {
        open(TABIX,"|-","tabix $regionoutfile");
        close(TABIX);
    }

    $i = -1 if ((++$cnt >= 10) && $debug);
    say "\$i = '$i'" if ($debug);
}

sub help {

say <<HELP;

    "debug"     =>  \$debug,
    "verbose"   =>  \$verbose,
    "infile:s"      =>  \$infile,
    "regionsize:i"  =>  \$regionsize,
    "refseqname:s"  
    "help"      =>  \$help,

HELP

}



