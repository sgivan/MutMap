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

my $script = '../bin/calc_SNPindex.pl';
my (@cmd_out,$cmd_out);

pass('calc_SNPindex.pl exists') if (-e -x '../bin/calc_SNPindex.pl');

ok(
    eval {
        open(CMD01,"-|","$script --help");
        $cmd_out = <CMD01>;
        close(CMD01);
    },'CMD01'
);

like($cmd_out, qr/usage/i, 'help menu');

ok(
    eval {
        open(CMD02,"-|", "$script --infile homo.vcf --window 20000");
        @cmd_out = <CMD02>;
        close(CMD02);
    }, 'CMD02'
);


is(scalar(@cmd_out),1,'one output line');

chomp($cmd_out[0]);
is($cmd_out[0],"53607763\t0.528430512072431\tconcatseq",'data valid');

ok(
    eval{
        open(CMD03,"-|","$script --infile homo.vcf --window 20000 --median");
        @cmd_out = <CMD03>;
        close(CMD03);
    }, 'CMD03'
);

chomp($cmd_out[0]);
is($cmd_out[0],"47036490\t0.528430512072431\tconcatseq",'median coordinate');

ok(
    eval {
        open(CMD03,"-|","$script --infile homo.vcf --window 2000 --median");
        chomp(@cmd_out = <CMD03>);
        close(CMD03);
    }, 'CMD03'
);

is(scalar(@cmd_out),3,'2K window data');

#say @cmd_out;

#26534690        0.542205934956505       concatseq
#53815659        0.521581173856048       concatseq
#89416463        0.517232315358291       concatseq

my (@col1,@col2,@col3);

for (my $i = 0; $i < scalar(@cmd_out); ++$i) {
    ($col1[$i], $col2[$i], $col3[$i]) = split "\t", $cmd_out[$i];
}

is($col1[0],26534690,'col1 val1');
is($col1[1],53815659,'col1 val2');
is($col1[2],89416463,'col1 val3');

is($col2[0],0.542205934956505,'col2 val1');
is($col2[1],0.521581173856048,'col2 val2');
is($col2[2],0.517232315358291,'col2 val2');

is($col3[0],'concatseq','col3 val1');
is($col3[1],'concatseq','col3 val2');
is($col3[2],'concatseq','col3 val2');

ok(
    eval {
        open(CMD04,"-|","$script --infile homo.vcf --window 2000 --median --ignore ignore.txt");
        chomp(@cmd_out = <CMD04>);
        close(CMD04);
    }, 'CMD04'
);

#27701888        0.541244060869365       concatseq
#66071011        0.525455122605403       concatseq
#91529414        0.500364552386319       concatseq

(@col1,@col2,@col3) = ();

for (my $i = 0; $i < scalar(@cmd_out); ++$i) {
    ($col1[$i], $col2[$i], $col3[$i]) = split "\t", $cmd_out[$i];
}

is($col1[0],27701888,'col1 val1');
is($col1[1],66071011,'col1 val2');
is($col1[2],91529414,'col1 val3');

is($col2[0],0.541244060869365,'col2 val1');
is($col2[1],0.525455122605403,'col2 val2');
is($col2[2],0.500364552386319,'col2 val2');

is($col3[0],'concatseq','col3 val1');
is($col3[1],'concatseq','col3 val2');
is($col3[2],'concatseq','col3 val2');


#[02/08/16 15:03:43] stahl test/$ ../bin/calc_SNPindex.pl --infile homo.vcf --window 2000 --maskedseq masked.fa
#38412004        0.536230545468286       concatseq
#87766671        0.506347647280669       concatseq

ok(
    eval {
        open(CMD05,"-|","$script --infile homo.vcf --window 2000 --maskedseq masked.fa");
        chomp(@cmd_out = <CMD05>);
        close(CMD05);
    }, 'CMD05'
);

(@col1,@col2,@col3) = ();

for (my $i = 0; $i < scalar(@cmd_out); ++$i) {
    ($col1[$i], $col2[$i], $col3[$i]) = split "\t", $cmd_out[$i];
}

is($col1[0],38412004,'col1 val1');
is($col1[1],87766671,'col1 val2');
#is($col1[2],91529414,'col1 val3');
ok(! exists($col1[2]), 'no value for row 3');

is($col2[0],0.536230545468286,'col2 val1');
is($col2[1],0.506347647280669,'col2 val2');
#is($col2[2],0.500364552386319,'col2 val2');
ok(! exists($col2[2]), 'no value for row 3');

is($col3[0],'concatseq','col3 val1');
is($col3[1],'concatseq','col3 val2');
#is($col3[2],'concatseq','col3 val2');
ok(! exists($col3[2]), 'no value for row 3');


done_testing();
