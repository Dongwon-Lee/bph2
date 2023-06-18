#!/bin/bash

set -o errexit
set -o nounset

GENOME=$1 #mm9, mm10, hg19, hg38
OUTPREFIX=$2 #roadmap_uw_hrt

shift
shift
FNs="$@"

module unload python
module load macs2

gsize=${GENOME:0:2}
if [[ "$gsize" == "hg" ]]; then
    gsize="hs"
fi
if [[ "$gsize" == "rn" ]]; then
    gsize="mm"
fi

macs2 callpeak -n $OUTPREFIX -g $gsize --nomodel --shift -50 --extsize 100 --keep-dup all -f BED -t $FNs
