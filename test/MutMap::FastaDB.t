#
#===============================================================================
#
#         FILE:  MutMap::FastaDB.t
#
#  DESCRIPTION:  Script to test MutMap::FastaDB module.
#
#        FILES:  ---
#         BUGS:  ---
#        NOTES:  ---
#       AUTHOR:  Dr. Scott A. Givan (sag), givans@missouri.edu
#      COMPANY:  University of Missouri, USA
#      VERSION:  1.0
#      CREATED:  02/17/16 17:20:47
#     REVISION:  ---
#===============================================================================

use 5.010;      # Require at least Perl version 5.10
use autodie;
use strict;
use warnings;
use Test::More;
use lib '../lib';

# declare number of tests to run
#use Test::More tests => 1;

use_ok( 'MutMap::FastaDB' );

my $fastafile = 'masked.fa';

my $fastaDB = MutMap::FastaDB->new( filename => $fastafile );

isa_ok($fastaDB, 'MutMap::FastaDB');

my $db = $fastaDB->makeFastaDB();

isa_ok($db,'Bio::DB::Fasta');

done_testing();
