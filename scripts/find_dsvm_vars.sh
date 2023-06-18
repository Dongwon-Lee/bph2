#!/bin/bash

set -o errexit
set -o nounset

dsvmf=$1 # ex) ../models/dsvm.thyroid_atac_r1.1kgp3.hg19.eur.maf1.bed 
cref=$2 # ex) ../all_heart_dhs_summits.e500.extend.500.bed
percentile=$3 # 15 (recommended)

annotdsvmf=`basename $dsvmf`
annotdsvmf=${annotdsvmf%.bed}.${expid}.pos.bed
thfn=${dsvmf%.bed}.${percentile}.cutoff.txt

if [ -e $thfn ]; then
    echo "skip calculating dsvm threshold..."
    threshold=`cat $thfn`
else
    n=`cat $dsvmf| wc -l`
    l=`echo "$n*$percentile/100" | bc`
    threshold=`cut -f 8 $dsvmf |sed 's/^-//g' |sort -gr -S 3G |head -n $l |tail -n 1`
    echo "$threshold" >$thfn
fi

# select SNPs overlapping CREs AND above threshold
intersectBed -a $dsvmf -b $cref |awk -v TH=$threshold '$8<(-1)*TH || $8>TH' >$annotdsvmf
