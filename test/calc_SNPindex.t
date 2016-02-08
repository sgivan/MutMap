#
#===============================================================================
#
#         FILE:  calc_SNP.t
#
#  DESCRIPTION:  Test script for calc_SNP.pl
#
#        FILES:  ---
#         BUGS:  ---
#        NOTES:  ---
#       AUTHOR:  Dr. Scott A. Givan (sag), givans@missouri.edu
#      COMPANY:  University of Missouri, USA
#      VERSION:  1.0
#      CREATED:  02/08/16 06:15:41
#     REVISION:  ---
#===============================================================================

use 5.010;      # Require at least Perl version 5.10
use autodie;
use strict;
use warnings;

# declare number of tests to run
#use Test::More tests => 1;
use Test::More;
# 
use lib '/share/apps/perl5/vcftools/lib/site_perl/5.14.2';

use_ok('Vcf');

my $script = '../bin/calc_SNPindex.pl';

pass('calc_SNPindex.pl exists') if (-e -x '../bin/calc_SNPindex.pl');

open(CMD01,"-|","$script --help");
my $cmd_out = <CMD01>;
close(CMD01);

like($cmd_out, qr/usage/i, 'help menu');

open(CMD02,"-|", "$script --infile homo.vcf --window 20000");
my @cmd_out = <CMD02>;
close(CMD02);

is(scalar(@cmd_out),1,'one output line');

chomp($cmd_out[0]);
is($cmd_out[0],"53607763\t0.528430512072431\tconcatseq",'data valid');

open(CMD03,"-|","$script --infile homo.vcf --window 20000 --median");
@cmd_out = <CMD03>;
close(CMD03);

chomp($cmd_out[0]);
is($cmd_out[0],"47036490\t0.528430512072431\tconcatseq",'median coordinate');


done_testing();
