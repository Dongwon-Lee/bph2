#!/bin/bash

#SBATCH --job-name=macs2.cp
#SBATCH --time=04:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --partition=cpu_short
#SBATCH --mem=8G

set -o errexit
set -o nounset

GENOME=$1 #mm or hs
OUTPREFIX=$2

shift
shift
FNs="$@"

macs2 callpeak -n $OUTPREFIX -g $GENOME --nomodel --shift -50 --extsize 100 --keep-dup all -f BED -t $FNs 2>${OUTPREFIX}.log
