#!/bin/bash

window=
increment=
while getopts w:i: opt; do
    case $opt in
        w)
            window=$OPTARG
            ;;
        i)
            increment=$OPTARG
            ;;
    esac
done
shift $((OPTIND - 1))

#window=$1
#increment=$2
#
#if [[ ! $window ]] ; then
#    window=500000
#else
#    shift
#fi
#
#if [[ ! $increment ]] ; then
#    increment=50000
#else
#    shift
#fi

echo "window: $window"
echo "increment: $increment"

#exit

for dir
do
    echo $dir
    start=$(echo $dir | cut -f 2 -d '.' | cut -f 1 -d '-')
    stop=$(echo $dir | cut -f 2 -d '.' | cut -f 2 -d '-')
    #echo "start: $start"
    #echo "stop: $stop"
    
    bsub -u givans@missouri.edu -B -N -R 'span[hosts=1] rusage[mem=10000]' -o %J.o -e %J.e -J ${dir}.${start}-${stop} "~/projects/mutmap/bin/calc_SNPindex+.pl \
    --ignore B2SNPindex_w1_gff_masked_overlap_AD.1.txt --maskedseq ForrestConcat.fa.masked --homoinfile $dir \
    --hetinfile het.bcf.gz --refseqfile ForrestConcat.fa --tab --gff Gmax_cds.blastn.gff \
    --mincoord --increment $increment --window $window --refseqname concatseq --start $start --stop $stop  > ${dir}.SNPindex.stdout"


#    exit # use when testing
done


