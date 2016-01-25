#!/usr/bin/env perl 
#===============================================================================
#
#         FILE:  calc_SNPindex+.pl
#
#        USAGE:  ./calc_SNPindex+.pl  
#
#
#  DESCRIPTION:  Calculate the SNP index for the alleles in a vcf file.
#                   The SNP index is defines as "the ratio between the
#                   number of reads of a mutant SNP and the total number
#                   of reads corresponding to the SNP. The specific
#                    algorithm used roughly corresponds to the MutMap+
#                    protocol, as described here;
#                    DOI: 10.1371/journal.pone.0068529
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

my
($debug,$verbose,$help,$homoinfile,$hetinfile,$refseqfile,$outfile,$window,$increment,$blockcnt,$tabstdout,$usemax,$maskedseqs,$gff,$overlap,$usemedian,$ignorefile,$onlyones,$usemaxcoord);

my $result = GetOptions(
    "homoinfile:s"  =>  \$homoinfile,
    "hetinfile:s"  =>  \$hetinfile,
    "refseq:s"      =>  \$refseqfile,
    "outfile:s" =>  \$outfile,
    "window:i"  =>  \$window,
    "increment:i"   =>  \$increment,
    "max"       =>  \$usemax,
    "maskedseqs:s"   =>  \$maskedseqs,
    "ignore:s"      => \$ignorefile,
    "gff:s"         =>  \$gff,
    "tab"       =>  \$tabstdout,
    "overlap"       =>  \$overlap,
    "median"        =>  \$usemedian,
    "maxcoord"      =>  \$usemaxcoord,
    "onlyones"      =>  \$onlyones,
    "debug"     =>  \$debug,
    "verbose"   =>  \$verbose,
    "help"      =>  \$help,
);

if ($help) {
    help();
    exit(0);
}

$window = 4000000 unless ($window);
$increment = 10000 unless ($increment);

#my $homovcf = Vcf->new(file=>$homoinfile);
#$homovcf->parse_header();
#
#my %homoDB = ();
#say "reading homo vcf file" if ($debug);
#while (my $x=$homovcf->next_data_array()) {
#    my $info_column = $homovcf->get_column($x,'INFO');
#    my $DP4 = $homovcf->get_info_field($info_column,'DP4');
#    my $refID = $homovcf->get_column($x,'CHROM');
#    my $coord = $homovcf->get_column($x,'POS');
#    my $nt = $homovcf->get_column($x,'ALT');
#    my $hashkey = $refID . ":" . $coord . ":" . $nt;
#    my @read_cnt = split /,/, $DP4;
#   
#    my $SNP_index = ($read_cnt[2] + $read_cnt[3])/($read_cnt[0] + $read_cnt[1] + $read_cnt[2] + $read_cnt[3]);
#
#    $homoDB{$hashkey} = $SNP_index;
#}
#
#say "\%homoDB created with ", keys(%homoDB), " keys" if ($debug);
#
#my $hetvcf = Vcf->new(file => $hetinfile);
#$hetvcf->parse_header();
#
#my %hetDB = ();
#say "reading het vcf file" if ($debug);
#while (my $y = $hetvcf->next_data_array()) {
#    my $refID = $hetvcf->get_column($y,'CHROM');
#    my $coord = $hetvcf->get_column($y,'POS');
#    my $nt = $hetvcf->get_column($y,'ALT');
#    my $hashkey = $refID . ":" . $coord . ":" . $nt;
#
#    my $info_column = $hetvcf->get_column($y,'INFO');
#    my $DP4 = $hetvcf->get_info_field($info_column,'DP4');
#    my @read_cnt = split /,/, $DP4;
#    my $SNP_index = ($read_cnt[2] + $read_cnt[3])/($read_cnt[0] + $read_cnt[1] + $read_cnt[2] + $read_cnt[3]);
#
#    $hetDB{$hashkey} = $SNP_index;
#
#}
#say "\%hetDB created with ", keys(%hetDB), " keys" if ($debug);

my $refseqIO = Bio::SeqIO->new(
    -file   =>  $refseqfile,
    -format =>  'fasta',
);
say "\$refseqIO is a '", ref($refseqIO), "'" if ($debug);

my $db;
if ($maskedseqs) {
    # Bio::DB::SeqFeature::Store
    $db      = Bio::DB::Fasta->new($maskedseqs);
    if ($debug) {
        say "reading masked sequence files in '$maskedseqs'";
        say "\$db isa " . ref($db);
    }
}

my %ignoreDB = (); 
# expect file with coordinate in column 1 and sequence ID in column 3
if ($ignorefile) {
    open(my $IG, "<", $ignorefile);
    for my $line (<$IG>) {
        chomp($line);
        my @linevals = split "\t", $line;
        # keys are "refID:coordinate"
        ++$ignoreDB{$linevals[2] . ":" . $linevals[0]};
    }
    close($IG);
}

if ($debug) {
    my $cnt = 0;
    for my $key (keys(%ignoreDB)) {
        say "ignoreDB key: '$key'";
        last if (++$cnt > 10);
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
while (my $refseq = $refseqIO->next_seq()) {
    my $refseqlength = $refseq->length();
    my $start = 1;
    my $refID = $refseq->id();
    say "refseq: '$refID'. length: '$refseqlength'" if ($debug);
    # $window is defined above
    # $increment is defined above
    my $window_SNP_index = 0;
    for (my ($j,$k) = ($start,$start + $window - 1); $j <= $refseqlength; $j += $increment, $k += $increment) {
        print "\n" if ($debug);
        last if (++$cnt >= 100 && $debug);
        say "cnt = $cnt" if ($debug);
        
        say "\$j = '$j', \$k = '$k'" if ($debug);
        my $homo_region = "$refID:$j" . "-" . $k;
        say "homoregion = '$homo_region'" if ($debug);
        
        my $homovcf = Vcf->new( file => $homoinfile, region => $homo_region);
        $homovcf->parse_header();
        say "\$homovcf isa '", ref($homovcf), "'" if ($debug);
        
        my $homo_data_array = 0;
        while (my $x=$homovcf->next_data_array()) {
            ++$homo_data_array;
            say "homo data array count: '$homo_data_array'" if ($debug);
            my $homo_info_column = $homovcf->get_column($x,'INFO');
            my $homoDP4 = $homovcf->get_info_field($homo_info_column,'DP4');
            my $refID = $homovcf->get_column($x,'CHROM');
            my $homo_coord = $homovcf->get_column($x,'POS');
            my $homo_alt = $homovcf->get_column($x,'ALT');
            #say "coord = '$coord'" if ($debug);
            #say "refID = '$refID'" if ($debug);
            my @homo_read_cnt = split /,/, $homoDP4;

            my $homo_SNP_index = ($homo_read_cnt[2] + $homo_read_cnt[3])/($homo_read_cnt[0] + $homo_read_cnt[1] + $homo_read_cnt[2] + $homo_read_cnt[3]);

            say "homo SNP identified at coordinate $homo_coord" if ($debug);

            if ($maskedseqs) {
                $seq = $db->get_Seq_by_id($refID);
                my $nt = $seq->subseq($homo_coord => $homo_coord);
                if ($nt) {
                    say "nt from masked file: '$nt'" if ($debug);
                    if ($nt eq 'N') {
                        say "skipping" if ($debug);
                        next;
                    }
                }
            }

            if ($ignorefile) {
                my $hashkey = $refID . ":" . $homo_coord;
                say "checking if I should ignore '$hashkey'" if ($debug);
                if (exists($ignoreDB{$hashkey})) {
                    say "ignoring potential SNP at homo_coordinate " . $homo_coord . " because it is in '$ignorefile'" if ($debug);
                    next;
                }
            }

            if ($gff) {
                say "checking for features that overlap location " . $homo_coord if ($debug);
                my @tfeatures = $gffdb->get_features_by_location(
                    -seq_id     =>  $refID,
                    -start      =>  $homo_coord,
                    -end        =>  $homo_coord,
                );

                if (scalar(@tfeatures)) {
                    say "overlapping features at location " . $homo_coord . ": " . scalar(@tfeatures) if ($debug);
                } else {
                    say "no overlapping features at location " . $homo_coord if ($debug);
                    next;
                }

            }
            
#            if ($debug) {
#                say "\nwindow cnt = $windowcnt";
#                say "i = $i";
#                printf "%i\t%s\t%3.2f\n", $coord, $DP4, $homo_SNP_index;
#            }

            # only instantiate new Vcf object once above filters are passed
            my $hetvcf = Vcf->new( file => $hetinfile, region => "$refID:$homo_coord" . "-" . $homo_coord);
            $hetvcf->parse_header();
            if ($hetvcf) {
                
                my $het_data_array = 0;
                while (my $y = $hetvcf->next_data_array()) {
                    ++$het_data_array;
                    say "het data arrays: '$het_data_array'" if ($debug);
                    my $het_alt = $hetvcf->get_column($y,'ALT');
                    if ($het_alt eq $homo_alt) {
                        say "het SNP allele matches: '$homo_alt' == '$het_alt'" if ($debug);
                        my $het_info_column = $hetvcf->get_column($y,'INFO');
                        my $hetDP4 = $hetvcf->get_info_field($het_info_column,'DP4');
                        my @het_read_cnt = split /,/, $hetDP4;
                        my $het_SNP_index = ($het_read_cnt[2] + $het_read_cnt[3])/($het_read_cnt[0] + $het_read_cnt[1] + $het_read_cnt[2] + $het_read_cnt[3]);
                        say "modifying value of homo_SNP_index: $homo_SNP_index " . "-" . " $het_SNP_index = ", ($homo_SNP_index - $het_SNP_index) if ($debug);
                        $homo_SNP_index -= $het_SNP_index;
                        last;
                    } else {
                        say "het SNP allele does not match: '$homo_alt' == '$het_alt'" if ($debug);
                    }
                }
            }
            say "pushing onto SNPlocalblock [ $homo_coord $homoDP4 $homo_SNP_index $refID ]" if ($debug);
            push(@SNPlocalblock, [$homo_coord, $homoDP4, $homo_SNP_index, $refID, $k]);

#            if (0) {
#
#                 build array of windows (array of arrays)
#                last if ($debug && $windowcnt > 1000);
#                if ($i >= $window) {
#                    if (!$overlap) {
#                        say "pushing " . dumper(@snplocalblock) . " onto snpblock" if ($debug);
#                        push(@snpblock,[@snplocalblock]);
#                        @snplocalblock = ();
#                        $i = 0;
#                        redo;
#                    } else {
##                         treat @snplocalblock as a filo 5 element buffer
#                        push(@snpblock,[@snplocalblock]);
#                        say "before shift: \@snplocalblock: " . dumper(@snplocalblock) if ($debug);
#                        shift(@snplocalblock);
#                        say "after shift: \@snplocalblock: " . dumper(@snplocalblock) if ($debug);
#                        push(@snplocalblock,[$coord, $dp4, $homo_snp_index, $refid]);
#                        say "after push: \@snplocalblock: " . dumper(@snplocalblock) if ($debug);
#                    }
#                    ++$windowcnt;
#                } elsif ($i < $window) {
#                    say "adding to snplocalblock" if ($debug);
#                    push(@snplocalblock, [$coord, $dp4, $homo_snp_index, $refid]);
#                }
#
#                ++$i;
#            }
        } # end of loop through homo Vcf objects.

        say "pushing on to SNPblock: " . scalar(@SNPlocalblock) . " element array" if ($debug);
        push(@SNPblock,[@SNPlocalblock]);
        @SNPlocalblock = ();
        print "\n\n\n" if ($debug);
    }
}

# if there is data left over ...
#if (@SNPlocalblock) {
#    push(@SNPblock,[@SNPlocalblock]);
#}

if ($debug) {
    say "Dump of \@SNPblock: ";
    print Dumper(sort {$a->[0] <=> $b->[0]} @SNPblock);
    say "window count: $windowcnt\n\n";

}

# now build data structure to facilitate plotting the SNPindex
#
# follow general algorithm in paper specified by doi, above
# take average SNPindex within a window and plot at middle of coordinate range
#

my $stat = Statistics::Descriptive::Full->new();
my @plotdata = ();
for my $windowdata (@SNPblock) {

    my ($x,$y,$refmol) = ();
    my $lastidx = @$windowdata - 1;
    #say "lastidx = $lastidx";
    for my $datapt (@$windowdata) {
        if ($usemaxcoord) {
            $stat->add_data($datapt->[4]);
        } else {
            $stat->add_data($datapt->[0]);
        }
        $refmol = $datapt->[3] unless ($refmol);
    }
    #$stat->add_data($windowdata->[0]->[0], $windowdata->[$lastidx]->[0]);
    say "SNPs in window: " . $stat->count() if ($debug);
    if ($usemedian) {
        say "median coordinate: " . int($stat->median()) if ($debug);
        $x = int($stat->median());
    } elsif ($usemaxcoord) {
        say "max coordinate: " . int($stat->max()) if ($debug);
        $x = int($stat->max());
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
    push(@plotdata, [$x, $y, $refmol]);
}

if ($debug) {
    say "data to plot:";
    print Dumper(@plotdata);
}

if ($tabstdout) {
    for my $xy (@plotdata) {
        if ($onlyones) {
            next if ($xy->[1] < 1);
        }
        say $xy->[0] . "\t" . $xy->[1] . "\t" . $xy->[2]
    }
}


sub help {

say <<HELP;

    "homoinfile:s"  =>  VCF file of homozygous pool
    "hetinfile:s"  =>  VCF file of het/WT pool
    "refseq:s"      =>  fasta file of reference sequence(s) used for VCF files
    "outfile:s" =>  \$outfile,
    "window:i"  =>      size of sliding window (default = 4M nt)
    "increment:i"   =>  increment amount (default = 10K nt)
    "max"           =>  instead of average in a window, use the maximum SNP index
    "maskedseq:s"   =>  name of FASTA file containing masked sequences to ignore
    "ignore:s"      =>  name of tab-delimited file containing coordinates to ignore in first column and reference sequence ID in column 3
    "gff:s"         =>  name of GFF3 file containing sequences to include -- happens after maskedseq.
    "tab"           =>  \$tabstdout,
    "overlap"       =>  overlap windows - default behavior does not overlap windows
    "median"        =>  use the median as the window genome coordinate instead of the mean
    "maxcoord"      =>  for plotting, use last coordinate of genome window 
    "onlyones"      =>  only output sites with SNPindex = 1
    "debug"         =>  \$debug,
    "verbose"       =>  \$verbose,
    "help"          =>  \$help,

    --maskedseq & --ignore are negative filters. They cause the script to ignore potential SNPs in the masked regions or specified coordinates.
    --gff is a positive filter. It causes the script to accept potential SNPs only in the specified regions.

HELP

}

