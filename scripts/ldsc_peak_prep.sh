#!/bin/bash

set -o errexit
set -o nounset

module load bedtools/2.27.1

fn0=$1 # narrowpeak file
expid=$2 # ENCFFxxxx_(ATAC)

fn1=${expid}.hg19.bed
fn2=${expid}.hg19.qc.bed
fn3=${expid}.hg19.qc.fe1000.bed
fn4=${expid}.hg19.qc.fe1000.merged.bed

liftOver ${fn0} -bedPlus=5 \
    ~/projects/shared_data/liftOver/hg38ToHg19.over.chain.gz \
    ${fn1} \
    unmapped

awk '$1 ~ /chr[0-9XY]+$/' ${fn1} >${fn2}

EXT=1000
awk -v OFS="\t" -v shft=$EXT \
    '$2>shft{summit=$2+$10; print $1,summit-shft,summit+shft,$4,$8}' \
    ${fn2} >${fn3}

cat ${fn3} |sortBed | mergeBed >${fn4}

