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

my $mutmap = MutMap->new(fastafile=>'masked.fa');

isa_ok($mutmap,'MutMap');

$mutmap->makeFastaDB();

is($mutmap->fastadb->filename(),'masked.fa',"filename OK");

done_testing();
