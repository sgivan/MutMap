package MutMap;
#
#===============================================================================
#
#         FILE:  MutMap.pm
#
#  DESCRIPTION:  
#
#        FILES:  ---
#         BUGS:  ---
#        NOTES:  ---
#       AUTHOR:  Dr. Scott A. Givan (sag), givans@missouri.edu
#      COMPANY:  University of Missouri, USA
#      VERSION:  1.0
#      CREATED:  02/17/16 17:29:46
#     REVISION:  ---
#===============================================================================

use 5.010;
use strict;
use warnings;
use autodie;
use Moose;
use Data::Dump qw/dump/;

1;

has 'fastadb' => (
    is      => 'ro',
    isa     =>  'MutMap::FastaDB',
    lazy_build  =>  1,
    handles =>  [qw( makeFastaDB )],
);

sub _build_fastadb {
    my $self = shift;
    my $fastafile = $self->{_fastafile};
    require MutMap::FastaDB;

    my $fastadb = MutMap::FastaDB->new(filename => $fastafile);
    $fastadb->makeFastaDB();
    return $fastadb;
}

sub BUILD {
    my $self = shift;
    my $args = shift;

#    say dump($args);

    $self->{_fastafile} = $args->{fastafile};
}

no Moose;
__PACKAGE__->meta->make_immutable;
