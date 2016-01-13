#!/usr/bin/env perl 
#===============================================================================
#
#         FILE:  calc_SNPindex.pl
#
#        USAGE:  ./calc_SNPindex.pl  
#
#
#  DESCRIPTION:  Calculate the SNP index for the alleles in a vcf file.
#                   The SNP index is defines as "the ratio between the
#                   number of reads of a mutant SNP and the total number
#                   of reads corresponding to the SNP.
#                   doi:10.1038/nbt.2095
#
#      OPTIONS:  ---
# REQUIREMENTS:  ---
#         BUGS:  ---
#        NOTES:  ---
#       AUTHOR:  Dr. Scott A. Givan (sag), givans@missouri.edu
#      COMPANY:  University of Missouri, USA
#      VERSION:  1.0
#      CREATED:  12/16/15 22:10:01
#     REVISION:  ---
#===============================================================================

use 5.010;      # Require at least Perl version 5.10
use autodie;
use Getopt::Long; # use GetOptions function to for CL args
use warnings;
use strict;
use Data::Dumper;
use Statistics::Descriptive;
use Bio::SeqIO;
use Bio::DB::Fasta;
#use Bio::DB::GFF;
use Bio::DB::SeqFeature::Store;
use Bio::DB::SeqFeature::Store::GFF3Loader;
use lib '/share/apps/perl5/vcftools/lib/site_perl/5.14.2';
use Vcf;
$Data::Dumper::Deepcopy = 1;

my ($debug,$verbose,$help,$infile,$outfile,$window,$blockcnt,$tabstdout,$usemax,$maskedseqs,$gff,$overlap,$usemedian);

my $result = GetOptions(
    "infile:s"  =>  \$infile,
    "outfile:s" =>  \$outfile,
    "window:i"  =>  \$window,
    "max"       =>  \$usemax,
    "maskedseqs:s"   =>  \$maskedseqs,
    "gff:s"         =>  \$gff,
    "tab"       =>  \$tabstdout,
    "overlap"       =>  \$overlap,
    "median"        =>  \$usemedian,
    "debug"     =>  \$debug,
    "verbose"   =>  \$verbose,
    "help"      =>  \$help,
);

if ($help) {
    help();
    exit(0);
}

$infile = 'infile' unless ($infile);
$window = 5 unless ($window);

my $vcf = Vcf->new(file=>$infile);
$vcf->parse_header();

my $db;
if ($maskedseqs) {
    # Bio::DB::SeqFeature::Store
    $db      = Bio::DB::Fasta->new($maskedseqs);
    if ($debug) {
        say "reading masked sequence files in '$maskedseqs'";
        say "\$db isa " . ref($db);
    }
}

my $gffdb;
if ($gff) {
    # Bio::DB::SeqFeature::Store
    say "loading GFF file '$gff'" if ($debug);
    $gffdb = Bio::DB::SeqFeature::Store->new( -adaptor => 'memory', -dsn => $gff);
#    my $loader = Bio::DB::SeqFeature::Store::GFFLoader->new( -store => $gffdb );
#    $loader->load($gff);

    if ($debug) {
        say "loaded GFF file";
        say "\$gffdb isa '" . ref($gffdb) . "'";
    }
}

my ($cnt,$i,$windowcnt) = (0,0,1);
my (@SNPblock,@SNPlocalblock) = ();
my $seq;
while (my $x=$vcf->next_data_array()) {
    my $info_column = $vcf->get_column($x,'INFO');
    my $DP4 = $vcf->get_info_field($info_column,'DP4');
    my $refID = $vcf->get_column($x,'CHROM');
    #say "refID = '$refID'" if ($debug);
    my @read_cnt = split /,/, $DP4;
   
    my $SNP_index = ($read_cnt[2] + $read_cnt[3])/($read_cnt[0] + $read_cnt[1] + $read_cnt[2] + $read_cnt[3]);

    if ($maskedseqs) {
        $seq = $db->get_Seq_by_id($refID);
        my $nt = $seq->subseq($x->[1] => $x->[1]);
        say "nt from masked file: '$nt'" if ($debug);
        if ($nt eq 'N') {
            say "skipping" if ($debug);
            next;
        }
    }

    if ($gff) {
        say "checking for features that overlap location " . $x->[1] if ($debug);
        my @tfeatures = $gffdb->get_features_by_location(
            -seq_id     =>  $refID,
            -start      =>  $x->[1],
            -end        =>  $x->[1],
        );

        if (scalar(@tfeatures)) {
            say "overlapping features at location " . $x->[1] . ": " . scalar(@tfeatures) if ($debug);
        } else {
            say "no overlapping features at location " . $x->[1] if ($debug);
            next;
        }

    }
    
    if ($debug) {
        say "\nwindow cnt = $windowcnt";
        say "i = $i";
        printf "%i\t%s\t%3.2f\n", $x->[1], $DP4, $SNP_index;
    }

    # build array of windows (array of arrays)
    last if ($debug && $windowcnt > 5);
    if ($i >= $window) {
        if (!$overlap) {
            say "pushing " . Dumper(@SNPlocalblock) . " onto SNPblock" if ($debug);
            push(@SNPblock,[@SNPlocalblock]);
            @SNPlocalblock = ();
            $i = 0;
            redo;
        } else {
            # treat @SNPlocalblock as a filo 5 element buffer
            push(@SNPblock,[@SNPlocalblock]);
            say "before shift: \@SNPlocalblock: " . Dumper(@SNPlocalblock) if ($debug);
            shift(@SNPlocalblock);
            say "after shift: \@SNPlocalblock: " . Dumper(@SNPlocalblock) if ($debug);
            push(@SNPlocalblock,[$x->[1], $DP4, $SNP_index]);
            say "after push: \@SNPlocalblock: " . Dumper(@SNPlocalblock) if ($debug);
        }
        ++$windowcnt;
    } elsif ($i < $window) {
        say "adding to SNPlocalblock" if ($debug);
        push(@SNPlocalblock, [$x->[1], $DP4, $SNP_index]);
    }

    ++$i;
}

# if there is data left over ...
if (@SNPlocalblock) {
    push(@SNPblock,[@SNPlocalblock]);
}

if ($debug) {
    say "Dump of \@SNPblock: ";
    print Dumper(@SNPblock);
    say "window count: $windowcnt";

}

# now build data structure to facilitate plotting the SNPindex
#
# follow general algorithm in paper specified by doi, above
# take average SNPindex within a window and plot at middle of coordinate range
#

my $stat = Statistics::Descriptive::Full->new();
my @plotdata = ();
for my $windowdata (@SNPblock) {

    my ($x,$y) = ();
    my $lastidx = @$windowdata - 1;
    #say "lastidx = $lastidx";
    for my $datapt (@$windowdata) {
        $stat->add_data($datapt->[0]);
    }
    #$stat->add_data($windowdata->[0]->[0], $windowdata->[$lastidx]->[0]);
    if ($usemedian) {
        say "median coordinate: " . int($stat->median()) if ($debug);
        $x = int($stat->median());
    } else {
        say "mean coordinate: " . int($stat->mean()) if ($debug);
        $x = int($stat->mean());
    }
    $stat->clear();

    for my $datapt (@$windowdata) {
        $stat->add_data($datapt->[2]);
    }

    if ($usemax) {
        say "max SNPindex: " . $stat->max() if ($debug);
        $y = $stat->max();
    } else {
        say "mean SNPindex: " . $stat->mean() if ($debug);
        $y = $stat->mean();
    }


    $stat->clear();
    push(@plotdata, [$x, $y]);
}

if ($debug) {
    say "data to plot:";
    print Dumper(@plotdata);
}

if ($tabstdout) {
    for my $xy (@plotdata) {
        say $xy->[0] . "\t" . $xy->[1];
    }
}


sub help {

say <<HELP;

    "infile:s"      =>  name of input file
    "outfile:s"     =>  name of output file
    "window:i"      =>  window size; ie, the number of SNP's to include in the SNP index calculation
    "max"           =>  instead of average in a window, use the maximum SNP index
    "maskedseq:s"   =>  name of FASTA file containing masked sequences
    "gff:s"         =>  name of GFF3 file containing sequences to include -- happens after maskedseq.
    "tab"           =>  \$tabstdout,
    "overlap"       =>  overlap windows - default behavior does not overlap windows
    "median"        =>  use the median as the window genome coordinate instead of the mean
    "debug"         =>  \$debug,
    "verbose"       =>  \$verbose,
    "help"          =>  \$help,

    --maskedseq is a negative filter. It causes the script to ignore potential SNPs in the masked regions.
    --gff is a positive filter. It causes the script to accept potential SNPs only in the specified regions.

HELP

}



