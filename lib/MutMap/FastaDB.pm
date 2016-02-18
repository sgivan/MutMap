package MutMap::FastaDB;
#
#===============================================================================
#
#         FILE:  FastaDB.pm
#
#  DESCRIPTION:  Module to handle reading/creating a database
#                   for a fasta file
#
#        FILES:  ---
#         BUGS:  ---
#        NOTES:  ---
#       AUTHOR:  Dr. Scott A. Givan (sag), givans@missouri.edu
#      COMPANY:  University of Missouri, USA
#      VERSION:  1.0
#      CREATED:  02/17/16 16:39:29
#     REVISION:  ---
#===============================================================================

use 5.010;
use strict;
use warnings;
use autodie;
use Bio::DB::Fasta;
use Moose;

1;

has 'filename' => (
    is      =>  'rw',
    isa     =>  'Str',
    default =>  'fastafile',
);

has 'debug' =>  (
    is      =>  'rw',
    isa     =>  'Int',
    default =>  0,
);

sub makeFastaDB {
    my $self = shift;

    my $db = Bio::DB::Fasta->new($self->filename());
    if ($self->debug) {
        say "reading masked sequences in '" . $self->filname . "'";
    }
    return $db;
}

no Moose;
__PACKAGE__->meta->make_immutable;
