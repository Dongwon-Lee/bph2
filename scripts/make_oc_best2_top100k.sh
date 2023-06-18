#!/bin/bash

set -o errexit
set -o nounset

expid1=$1
expid2=$2
outfile=$3

for expid in ${expid1} ${expid2}; do
    echo ${expid}
    fn=${expid}.hg19.qc.fe1000.bed
    cat $fn |sortBed | mergeBed -c 5 -o max |\
    sort -grk 4 |awk '{print $0"\t"NR}' |sortBed >${fn%.bed}.merged2.bed
done

cat ${expid1}.hg19.qc.fe1000.merged2.bed \
    ${expid2}.hg19.qc.fe1000.merged2.bed \
    |sortBed |mergeBed -c 5 -o mean \
    |sort -gk 4 |head -n 100000 \
    |sortBed >${outfile}
