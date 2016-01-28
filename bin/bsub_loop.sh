#!/bin/bash

window=$1
increment=$2

if [[ ! $window ]] ; then
    window=500000
else
    shift
fi

if [[ ! $increment ]] ; then
    increment=50000
else
    shift
fi

echo "window: $window"
echo "increment: $increment"

for dir in *.bcf.gz
do
    echo $dir
    start=$(echo $dir | cut -f 2 -d '.' | cut -f 1 -d '-')
    stop=$(echo $dir | cut -f 2 -d '.' | cut -f 2 -d '-')
    #echo "start: $start"
    #echo "stop: $stop"
    # bsm -J ${start}.${stop} '~/projects/mutmap/bin/calc_SNPindex+.pl --ignore ../B2SNPindex_w1_gff_masked_overlap_AD.1.txt --maskedseq \
    # ../../masked_seqs/ForrestConcat.fa.masked --homoinfile $dir --hetinfile ../het.bcf.gz --refseqfile \
    # ../../index/ForrestConcat.fa --tab --gff ../Gmax_cds.blastn.gff --mincoord --increment 500000 --window 8000000 --refseqname concatseq \
    # --start $start --stop $stop  > ${dir}.SNPindex'

    #bsub -u givans@missouri.edu -B -N -R 'span[hosts=1]' -o %J.o -e %J.e -J ${dir}.${start}-${stop} "~/projects/mutmap/bin/calc_SNPindex+.pl \
    #--ignore ../B2SNPindex_w1_gff_masked_overlap_AD.1.txt --maskedseq ../../masked_seqs/ForrestConcat.fa.masked --homoinfile $dir \
    #--hetinfile ../het.bcf.gz --refseqfile ../../index/ForrestConcat.fa --tab --gff ../Gmax_cds.blastn.gff \
    #--mincoord --increment 50000 --window 1000000 --refseqname concatseq --start $start --stop $stop  > ${dir}.SNPindex"

    #bsub -u givans@missouri.edu -B -N -R 'span[hosts=1]' -o %J.o -e %J.e -J ${dir}.${start}-${stop} "~/projects/mutmap/bin/calc_SNPindex+.pl \
    #--ignore B2SNPindex_w1_gff_masked_overlap_AD.1.txt --maskedseq ForrestConcat.fa.masked --homoinfile $dir \
    #--hetinfile het.bcf.gz --refseqfile ForrestConcat.fa --tab --gff Gmax_cds.blastn.gff \
    #--mincoord --increment 50000 --window 500000 --refseqname concatseq --start $start --stop $stop  > ${dir}.SNPindex"
    
    #bsub -u givans@missouri.edu -B -N -R 'span[hosts=1]' -o %J.o -e %J.e -J ${dir}.${start}-${stop} "~/projects/mutmap/bin/calc_SNPindex+.pl \
    #--ignore B2SNPindex_w1_gff_masked_overlap_AD.1.txt --maskedseq ForrestConcat.fa.masked --homoinfile $dir \
    #--hetinfile het.bcf.gz --refseqfile ForrestConcat.fa --tab --gff Gmax_cds.blastn.gff \
    #--mincoord --increment $increment --window $window --refseqname concatseq --start $start --stop $stop  > ${dir}.SNPindex"


    exit # use when testing
done


