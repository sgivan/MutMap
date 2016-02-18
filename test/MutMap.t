#
#===============================================================================
#
#         FILE:  MutMap.t
#
#  DESCRIPTION:  Test script for MutMap class.
#
#        FILES:  ---
#         BUGS:  ---
#        NOTES:  ---
#       AUTHOR:  Dr. Scott A. Givan (sag), givans@missouri.edu
#      COMPANY:  University of Missouri, USA
#      VERSION:  1.0
#      CREATED:  02/18/16 09:30:49
#     REVISION:  ---
#===============================================================================

use 5.010;      # Require at least Perl version 5.10
use autodie;
use strict;
use warnings;
use lib '../lib';

# declare number of tests to run
use Test::More;

use_ok('MutMap');

# passing a fastafile to the constructor should build the index
my $mutmap = MutMap->new(fastafile=>'masked.fa');

isa_ok($mutmap,'MutMap');

is($mutmap->fastadb->filename(),'masked.fa',"filename passed in constructor");

my $fastadb = $mutmap->fastadb();

isa_ok($fastadb,'MutMap::FastaDB');

if (-e "masked.fa.index") {
    pass('fasta index created');
} else {
    fail('no fasta index created');
}

unlink("masked.fa.index");

my $mutmap2 = MutMap->new();

$mutmap2->fastadb->filename('masked.fa');

is($mutmap2->fastadb->filename(),'masked.fa','filename set with attribute method');

$mutmap2->makeFastaDB();

if (-e "masked.fa.index") {
    pass('fasta index created');
} else {
    fail('no fasta index created');
}

unlink("masked.fa.index");


done_testing();
