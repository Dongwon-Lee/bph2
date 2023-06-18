#!/bin/bash

#SBATCH --job-name=preproc_atacseq_bam
#SBATCH --time=11:59:59
#SBATCH --ntasks=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --partition=cpu_short
#SBATCH --mem=5G

set -o errexit
set -o nounset

BAMF=$1
OUTPREFIX=$2

DBAMF=${OUTPREFIX}.bam
cutpos_bedgz=${OUTPREFIX}.bed.gz

#MAPQ >= X   #MM Q40   #MM Q20      #MM Q0    Description
#0           5         7            15        All mappable reads
#1           3         5            10        True multi w/ "good" AS, maxi of MAPQ >= 1
#2           3         5            10        No true multi, maxi of MAPQ >= 2
#3           3         5            10        No true multi,  maxi of MAPQ >= 3
#8           2         4            8         No true multi, maxi of MAPQ >= 8
#23          2         3            7         No true multi, maxi of MAPQ >= 23
#30          1         2            4         No true multi, maxi of MAPQ >= 30
#39          1         2            4         No true multi, maxi of MAPQ == 39*
#40          1         2            4         No true multi, only true uni-reads
#42          0         1            2         Only "perfect" true unireads
MAPQ=20

# It has already been sorted and duplicated reads were removed... so skip these parts
# We just need to filter some low quality reads...
samtools view -b -f 0x2 -q $MAPQ $BAMF >$DBAMF
samtools stats ${DBAMF} >${OUTPREFIX}.samstats.txt

#remove chrM and any other random chromosomes containing '_'
#also shift +4/-5bp to find the cut positions of tagmentation reactions
samtools sort -n -m 2G -O BAM $DBAMF \
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
