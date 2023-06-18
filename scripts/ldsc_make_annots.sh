#!/bin/bash

#SBATCH --job-name=annot
#SBATCH --time=0:30:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --partition=cpu_short
#SBATCH --mem=4G

set -o errexit
set -o nounset

cref=$1 #../all_heart_dhs_summits.e500.extend.500.cardiovascular.bed
expid=$2 #ObsHrtCRE2kbCV

mkdir -p 1000G_EUR_Phase3_${expid}

#i=22
for i in `seq 22`; do
    echo $i
    annotf=~/projects/shared_data/ldsc/1000G_EUR_Phase3_baseline/baseline.${i}.annot.gz

    zcat $annotf |awk -v OFS="\t" 'NR>1{print "chr"$1,$2-1,$2,$3}' \
        |intersectBed -c -a - -b $cref \
        |awk -v EXPID=$expid 'BEGIN{print EXPID}{print $5}' \
        >annot.${i}.${expid}

    zcat $annotf |cut -f 1-4 |paste - annot.${i}.${expid} \
        |gzip -c >1000G_EUR_Phase3_${expid}/${expid}.${i}.annot.gz

    rm -f annot.${i}.${expid}
done
