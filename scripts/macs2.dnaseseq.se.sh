#!/bin/bash

#SBATCH --job-name=macs2.dhs
#SBATCH --time=04:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --partition=cpu_short
#SBATCH --mem=8G

set -o errexit
set -o nounset

GENOME=$1 #mm or hs
OUTPREFIX=$2
FN=$3 # bamfile

macs2 callpeak -n $OUTPREFIX -g $GENOME --nomodel --shift -50 --extsize 100 -f BAM -t $FN -B --SPMR 2>${OUTPREFIX}.log

gzip ${OUTPREFIX}_treat_pileup.bdg
gzip ${OUTPREFIX}_control_lambda.bdg
