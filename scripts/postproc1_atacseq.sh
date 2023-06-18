#!/bin/bash

set -o errexit
set -o nounset

module load samtools
module load picard

BAMF=$1
EXPNAME=${BAMF%.bam}
TBAMF=${EXPNAME}.tmp.bam
SBAMF=${EXPNAME}.flt.bam
DBAMF=${EXPNAME}.flt.rd.bam
METRICF=${EXPNAME}.flt.rd.metric.txt

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

# filtering & sorting
samtools view -b -f 0x2 -q $MAPQ $BAMF >$TBAMF

samtools sort -m 2000M -O BAM -o $SBAMF $TBAMF

samtools index $SBAMF

java -Xmx4G -jar $PICARD_DIR/libs/picard.jar MarkDuplicates I=$SBAMF O=$DBAMF M=$METRICF REMOVE_DUPLICATES=true

samtools index $DBAMF

# Collect bam stats
samtools stats ${BAMF} >${EXPNAME}.samstats.txt
samtools stats ${SBAMF} >${SBAMF%.bam}.samstats.txt
samtools stats ${DBAMF} >${DBAMF%.bam}.samstats.txt

rm -f $TBAMF $METRICF
