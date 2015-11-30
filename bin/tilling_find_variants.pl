#!/usr/bin/env perl 
#===============================================================================
#
#         FILE:  tilling_find_variants.pl
#
#        USAGE:  ./tilling_find_variants.pl  
#
#  DESCRIPTION:  Script to identify variants with ~2X higher prevalence in 
#                   one of two populations.
#
#      OPTIONS:  ---
# REQUIREMENTS:  ---
#         BUGS:  ---
#        NOTES:  ---
#       AUTHOR:  Dr. Scott A. Givan (sag), givans@missouri.edu
#      COMPANY:  University of Missouri, USA
#      VERSION:  1.0
#      CREATED:  11/30/15 13:59:49
#     REVISION:  ---
#===============================================================================

use 5.010;      # Require at least Perl version 5.10
use autodie;
use Getopt::Long; # use GetOptions function to for CL args
use warnings;
use strict;

my ($debug,$verbose,$help,$heterofile,$homofile,$threshold);

my $result = GetOptions(
    "debug"     =>  \$debug,
    "verbose"   =>  \$verbose,
    "help"      =>  \$help,
    "heterofile:s"      =>  \$heterofile,
    "homofile:s"        =>  \$homofile,
    "threshold:f"       =>  \$threshold,
);

if ($help) {
    help();
    exit(0);
}

unless ($heterofile && $homofile) {
    say "you must provide both --heterofile and --homofile";
    help();
    exit();
}

open(my $HET, "<", $heterofile);
open(my $HOMO, "<", $homofile);

my (%HET,%HOMO) = ();


while (<$HET>) {
    next unless ($_ !~ /^#/);
    my @vals = split /\t/, $_;

    my $hashkey = $vals[0] . ":" . $vals[1];

    if (!$HET{$hashkey}) {
       $HET{$hashkey} = [$vals[3], $vals[4]]; 
    } else {
        die "'$hashkey' present more than 1X";
    }
}

while (<$HOMO>) {
    next unless ($_ !~ /^#/);
    my @vals = split /\t/, $_;

    my $hashkey = $vals[0] . ":" . $vals[1];

    if (!$HOMO{$hashkey}) {
       $HOMO{$hashkey} = [$vals[3], $vals[4]]; 
    } else {
        die "'$hashkey' present more than 1X";
    }
}

if ($verbose) {
    say "het loci: " . keys(%HET);
    say "homo loci: " . keys(%HOMO);
}

for my $homokey (keys %HOMO) {
    if ($HET{$homokey}) {
        say "$homokey\tHET variant: " . $HET{$homokey}[0] . " -> " . $HET{$homokey}[1] . "\tHOMO variant: " . $HOMO{$homokey}[0] . " -> " . $HOMO{$homokey}[1];
    }
}

sub help {

say <<HELP;

    "debug"     =>  \$debug,
    "verbose"   =>  \$verbose,
    "help"      =>  \$help,
    "heterofile:s"      =>  \$heterofile,
    "homofile:s"        =>  \$homofile,
    "threshold:f"       =>  \$threshold,

HELP

}



