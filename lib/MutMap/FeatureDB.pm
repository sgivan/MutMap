package MutMap::FeatureDB;
#
#===============================================================================
#
#         FILE:  FeatureDB.pm
#
#  DESCRIPTION:  Class to create/query coordinates of genome features
#
#        FILES:  ---
#         BUGS:  ---
#        NOTES:  ---
#       AUTHOR:  Dr. Scott A. Givan (sag), givans@missouri.edu
#      COMPANY:  University of Missouri, USA
#      VERSION:  1.0
#      CREATED:  02/19/16 14:55:46
#     REVISION:  ---
#===============================================================================

use 5.010;
use strict;
use warnings;
use autodie;
use Moose;
use Bio::DB::SeqFeature::Store;

has 'featurefile' => (
    is      =>  'rw',
    isa     =>  'Str',
    default =>  'featurefile.gff',
);

sub makeFeatureDB {
    my $self = shift;


    my $gffdb = Bio::DB::SeqFeature::Store->new( -adaptor => 'memory', -dsn => $self->featurefile);

    return $gffdb;
}

no Moose;
__PACKAGE__->meta->make_immutable;
