#!/bin/bash

#SBATCH --job-name=preproc_dnase_bam
#SBATCH --time=11:59:59
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --partition=cpu_short
#SBATCH --mem=6G

set -o errexit
set -o nounset

if [ $# -ne 2 ]; then
    echo "Usage: preproc_dnase_bam.sh <filtered.bam> <output_root>"
    echo "1) Input: Filtered paired-end aligned reads of DNase-seq data from the ENCODE project"
    echo "2) Tasks:"
    echo " - filter out duplicated reads"
    echo " - extract the cut sites from each side"
    echo " - store them as <output_root>.bed.gz"
    exit -1; 
fi

filtered_bam=$1  # Stam's filtered DNase-seq BAM file from the ENCODE project
output_prefix=$2 # prefix for outputs. i.e. "out" will create "out.bed.gz"
output_bam=${output_prefix}.bam
cutpos_bedgz=${output_prefix}.bed.gz

echo "-- rmdup alignments file will be: '${output_bam}'"

echo "-- remove duplicated reads"
#    1 read paired
#    2 read mapped in proper pair
#    4 read unmapped
#    8 mate unmapped
#   16 read reverse strand
#   32 mate reverse strand
#   64 first in pair
#  128 second in pair
#  256 not primary alignment
#  512 read fails platform/vendor quality checks
# 1024 read is PCR or optical duplicate
# 2048 supplementary alignment

filter_flags=1024
samtools view -F ${filter_flags} -b ${filtered_bam} >${output_bam}

echo "-- Collect bam stats..."
samtools flagstat ${filtered_bam} >${filtered_bam%.bam}_flagstat.txt
samtools flagstat ${output_bam} >${output_prefix}_flagstat.txt
samtools stats ${output_bam} >${output_prefix}_samstats.txt
grep ^SN ${output_prefix}_samstats.txt | cut -f 2- > ${output_prefix}_samstats_summary.txt

echo "-- Collect bam stats..."

echo "-- extract cut positions..." 
# 1. sort by name
# 2. convert them to bampe format using bedtools
# 3. find the cut positions and convert to bed 
# 4. compress and save them

samtools sort -n -m 1024M -O BAM ${output_bam} \
| bamToBed -bedpe \
| awk -v OFS="\t" \
    '{
        r1c=$1; r1s=$2; r2e=$6; sn=$7;
        print r1c, r1s, r1s+1, sn;
        print r1c, r2e-1, r2e, sn;
    }' \
| gzip -c >${cutpos_bedgz}

if [ "$?" != "0" ]; then
    echo -e "An error occurred while processing the BAM file ${output_bam}; exiting."
fi
