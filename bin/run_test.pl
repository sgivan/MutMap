#!/usr/bin/env perl 
#===============================================================================
#
#         FILE:  run_test.pl
#
#        USAGE:  ./run_test.pl  
#
#  DESCRIPTION:  Run all test scripts
#
#      OPTIONS:  ---
# REQUIREMENTS:  ---
#         BUGS:  ---
#        NOTES:  ---
#       AUTHOR:  Dr. Scott A. Givan (sag), givans@missouri.edu
#      COMPANY:  University of Missouri, USA
#      VERSION:  1.0
#      CREATED:  02/17/16 15:31:18
#     REVISION:  ---
#===============================================================================

use 5.010;      # Require at least Perl version 5.10
use autodie;
use Getopt::Long; # use GetOptions function to for CL args
use warnings;
use strict;
use TAP::Harness;
use File::chdir;

$CWD = "../test";

my ($debug,$verbose,$help);

my $result = GetOptions(
    "debug"     =>  \$debug,
    "verbose"   =>  \$verbose,
    "help"      =>  \$help,
);

my %args = (
    verbosity =>    1,
);

my $harness = TAP::Harness->new(\%args);

$harness->rules({
        seq     =>  [
            { par   =>  ['*.t'], },
        ]
    }
);
my @test_scripts = qw( calc_SNPindex.t calc_SNPindex+.t );
my $aggregator = $harness->runtests(@test_scripts);

if ($help) {
    help();
    exit(0);
}

sub help {

    say <<HELP;


HELP

}



