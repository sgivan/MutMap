#!/usr/bin/env perl 
#===============================================================================
#
#         FILE:  calc_SNPindex.pl
#
#        USAGE:  ./calc_SNPindex.pl  
#
#
#  DESCRIPTION:  Calculate the SNP index for the alleles in a vcf file.
#                   The SNP index is defines as "the ratio between the
#                   number of reads of a mutant SNP and the total number
#                   of reads corresponding to the SNP.
#                   doi:10.1038/nbt.2095
#
#      OPTIONS:  ---
# REQUIREMENTS:  ---
#         BUGS:  ---
#        NOTES:  ---
#       AUTHOR:  Dr. Scott A. Givan (sag), givans@missouri.edu
#      COMPANY:  University of Missouri, USA
#      VERSION:  1.0
#      CREATED:  12/16/15 22:10:01
#     REVISION:  ---
#===============================================================================

use 5.010;      # Require at least Perl version 5.10
use autodie;
use Getopt::Long; # use GetOptions function to for CL args
use warnings;
use strict;
use lib '/share/apps/perl5/vcftools/lib/site_perl/5.14.2';
use Vcf;

my ($debug,$verbose,$help,$infile,$outfile);

my $result = GetOptions(
    "infile:s"  =>  \$infile,
    "outfile:s" =>  \$outfile,
    "debug"     =>  \$debug,
    "verbose"   =>  \$verbose,
    "help"      =>  \$help,
);

if ($help) {
    help();
    exit(0);
}

$infile = 'infile' unless ($infile);

my $vcf = Vcf->new(file=>$infile);
$vcf->parse_header();

my $cnt = 0;
while (my $x=$vcf->next_data_array()) {
    my $info_column = $vcf->get_column($x,'INFO');
    my $DP4 = $vcf->get_info_field($info_column,'DP4');
    #say $x->[1], "\t", $info_column, "\t", $DP4;
    my @read_cnt = split /,/, $DP4;
    my $SNP_index = ($read_cnt[2] + $read_cnt[3])/($read_cnt[0] + $read_cnt[1] + $read_cnt[2] + $read_cnt[3]);
    #say $x->[1], "\t", $DP4, "\t", $SNP_index;
    printf "%i\t%s\t%3.2f\n", $x->[1], $DP4, $SNP_index;
    last if ($debug && ++$cnt >= 10);
}

sub help {

    say <<HELP;


HELP

}



