#!/bin/bash

set -o errexit
set -o nounset

ANNOT=$1 #1000G_EUR_Phase3_AllHeartCRE/AllHeartCRE.${CHROM}.annot.gz
BFILE=$2 #1000G_EUR_Phase3_plink/1000G.EUR.QC.${CHROM}
SNPLIST=$3 #1000G_EUR_Phase3_baseline_snps/baseline_snp.${CHROM}.snp
OUTPREFIX=$4

python ~/ldsc/ldsc.py \
    --l2 \
    --bfile ${BFILE} \
    --ld-wind-cm 1 \
    --annot ${ANNOT} \
    --out ${OUTPREFIX} \
    --print-snps ${SNPLIST}
