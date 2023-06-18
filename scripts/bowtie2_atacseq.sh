#!/bin/bash

set -o errexit
set -o nounset

SEQF1=$1
SEQF2=$2
IDXPREFIX=$3 #~/projects/shared_data/bowtie2_index/mm10/mm10
OUTBAMF=$4
NCPUS=8
MAXIS=2000

#2. mapping
bowtie2 -p $NCPUS -X $MAXIS --dovetail --no-mixed --no-discordant -t -x $IDXPREFIX \
    -1 $SEQF1 -2 $SEQF2 | samtools view -b - >$OUTBAMF
