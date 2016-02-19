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

    my $fastadb = MutMap::FastaDB->new();;
    if ($fastafile) {
        $fastadb->filename($fastafile);
        $fastadb->makeFastaDB();
    }# else {
#        my $fastadb = MutMap::FastaDB->new();
#    }
    return $fastadb;
}

has 'ignoredb' => (
    is          =>  'ro',
    isa         =>  'MutMap::IgnoreDB',
    lazy_build  =>  1,
    handles     =>  [qw( makeIgnoreDB )],
);

sub _build_ignoredb {
    my $self = shift;
    my $ignorefile = $self->{_ignorefile};
    require MutMap::IgnoreDB;

    my $ignoredb = MutMap::IgnoreDB->new();
    if ($ignorefile) {
        $ignoredb->filename($ignorefile);
        $ignoredb->makeIgnoreDB();
    }
    return $ignoredb;
}

has 'featuredb' => (
    is          =>  'ro',
    isa         =>  'MutMap::FeatureDB',
    lazy_build  =>  1,
    handles     =>  [qw( makeFeatureDB )],
);

sub _build_featuredb {
    my $self = shift;
    require MutMap::FeatureDB;

    my $featuredb = MutMap::FeatureDB->new();

    return $featuredb;
}

sub BUILD {
    my $self = shift;
    my $args = shift;

#    say dump($args);

    $self->{_featurefile} = $args->{featurefile};
    $self->{_fastafile} = $args->{fastafile};
    $self->{_ignorefile} = $args->{ignorefile};
}

no Moose;
__PACKAGE__->meta->make_immutable;
