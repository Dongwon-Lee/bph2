#!/bin/bash

set -o errexit
set -o nounset

PAIRED_IN1=$1
PAIRED_IN2=$2
OUTDIR=$3

# parameters optimized for ATAC-seq
STRINGENCY=1 # Require MINLENGTH overlap between read and adapter for an adapter to be found
QUAL=20 # Trim low-quality bases from 5' and/or 3' ends of each read before adapter removal
LENGTH=35 # Discard reads shorter than this
ER=0.1 # Maximum allowed error rate as value between 0 and 1 (no. of errors divided by length of matching region).
ATACSEQ_ADAPTER=CTGTCTCTTATACACATCT

PAIRED_OUT1=`basename ${PAIRED_IN1%.fastq.gz}`
PAIRED_OUT2=`basename ${PAIRED_IN2%.fastq.gz}`
PAIRED_OUT1="$OUTDIR/${PAIRED_OUT1}.trimmed.fastq.gz"
PAIRED_OUT2="$OUTDIR/${PAIRED_OUT2}.trimmed.fastq.gz"

cutadapt -e $ER -O $STRINGENCY -q $QUAL,$QUAL -m $LENGTH \
    -a $ATACSEQ_ADAPTER -A $ATACSEQ_ADAPTER \
    -o ${PAIRED_OUT1} -p ${PAIRED_OUT2} ${PAIRED_IN1} ${PAIRED_IN2}
