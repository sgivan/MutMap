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
#use lib '/share/apps/perl5/vcftools/lib/site_perl/5.14.2';
use lib '../vcftools/src/perl/';

use_ok('Vcf');

my $script = '../bin/calc_SNPindex+.pl';
my (@cmd_out,$cmd_out);

pass('calc_SNPindex+.pl exists') if (-e -x $script);

ok(
    eval {
        open(CMD01,"-|","$script --help");
        $cmd_out = <CMD01>;
        close(CMD01);
    },'CMD01'
);

like($cmd_out, qr/usage/i, 'help menu');

diag("the following step may take several minutes to finish ...");

ok(
    eval {
        open(CMD02,"-|","$script --ignore ignore.txt --maskedseq masked.fa --homoinfile homo.80000001-88000000.bcf.gz --hetinfile het.bcf.gz --refseqfile refseq.fa --tab --gff filter.gff --mincoord --increment 200000 --window 2000000 --refseqname concatseq --start 80000001 --stop 88000000 > homo.80000001-88000000.bcf.gz.SNPindex");
        chomp(@cmd_out = <CMD02>);
        close(CMD02);
    }, 'CMD05'
);

open(CMPR,"-|","diff -q homo.80000001-88000000.bcf.gz.SNPindex homo.80000001-88000000.bcf.gz.SNPindex.ref");
chomp(@cmd_out = <CMPR>);
close(CMPR);

is(scalar(@cmd_out),0,'reference output confirmed');

unlink('homo.80000001-88000000.bcf.gz.SNPindex');

#(@col1,@col2,@col3) = ();
#
#for (my $i = 0; $i < scalar(@cmd_out); ++$i) {
#    ($col1[$i], $col2[$i], $col3[$i]) = split "\t", $cmd_out[$i];
#}
#
#is($col1[0],38412004,'col1 val1');
#is($col1[1],87766671,'col1 val2');
##is($col1[2],91529414,'col1 val3');
#ok(! exists($col1[2]), 'no value for row 3');
#
#is($col2[0],0.536230545468286,'col2 val1');
#is($col2[1],0.506347647280669,'col2 val2');
##is($col2[2],0.500364552386319,'col2 val2');
#ok(! exists($col2[2]), 'no value for row 3');
#
#is($col3[0],'concatseq','col3 val1');
#is($col3[1],'concatseq','col3 val2');
##is($col3[2],'concatseq','col3 val2');
#ok(! exists($col3[2]), 'no value for row 3');


done_testing();
