#!/bin/bash

set -o errexit
set -o nounset

SUMSTATFN=$1
EXPIDS=$2
WLD=$3 # weights_hm3_no_hla/weights.
BASELINE=$4 #1000G_Phase3_baselineLD_v2.2/baselineLD.
FREQ=$5 # 1000G_Phase3_frq/1000G.EUR.QC.
OUTID=$6


TESTS=${BASELINE}
EXPIDSTR=""
for expid in ${EXPIDS//,/ } 
do
    echo $expid
    TESTS="${TESTS},1000G_EUR_Phase3_${expid}/${expid}."
    EXPIDSTR="${EXPIDSTR}.${expid}"
done

echo $TESTS
echo $EXPIDSTR

python ~/ldsc/ldsc.py \
    --h2 ${SUMSTATFN}\
    --w-ld-chr ${WLD}\
    --ref-ld-chr ${TESTS}\
    --overlap-annot\
    --frqfile-chr ${FREQ}\
    --print-coefficients\
    --out ldsc_baseLD22_${OUTID}
