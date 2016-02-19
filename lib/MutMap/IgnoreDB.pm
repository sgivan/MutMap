package MutMap::IgnoreDB;
#
#===============================================================================
#
#         FILE:  IgnoreDB.pm
#
#  DESCRIPTION:  Class to manage/query ignored positions in genome.
#
#        FILES:  ---
#         BUGS:  ---
#        NOTES:  ---
#       AUTHOR:  Dr. Scott A. Givan (sag), givans@missouri.edu
#      COMPANY:  University of Missouri, USA
#      VERSION:  1.0
#      CREATED:  02/19/16 14:25:12
#     REVISION:  ---
#===============================================================================

use 5.010;
use strict;
use warnings;
use autodie;
use Moose;

1;

has ignorefile =>   (
    isa         =>  'Str',
    is          =>  'rw',
    default     =>  'ignorefile.txt',
);


sub makeIgnoreDB {
    my $self = shift;
    my $ignoreDB = shift;


    open(my $IG, "<", $self->ignorefile);
    for my $line (<$IG>) {
        chomp($line);
        my @linevals = split "\t", $line;
        # keys are "refID:coordinate"
        ++$ignoreDB->{$linevals[2] . ":" . $linevals[0]};
    }
    close($IG);
    return $ignoreDB;
}

no Moose;
__PACKAGE__->meta->make_immutable;
