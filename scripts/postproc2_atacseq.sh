#!/bin/bash

set -o errexit
set -o nounset

BAMF=$1 #atacseq_neuro2a_r1.srt.rmdup.bam

cutpos_bedgz=${BAMF%.bam}.bed.gz

#remove chrM and any other random chromosomes containing '_'
#also shift +4/-5bp to find the cut positions of tagmentation reactions
samtools sort -n -m 2G -O BAM $BAMF \
| bamToBed -bedpe \
| awk -v OFS="\t" \
    '$1 !~ "chrM" && $3 !~ "_" {
        r1c=$1; r1s=$2; r2e=$6; sn=$7;
        h1=r1s+4;
        h2=r2e-5;
        print r1c,h1,h1+1,sn;
        print r1c,h2-1,h2,sn;
    }' \
| gzip -c >${cutpos_bedgz}
