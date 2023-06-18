#!/bin/bash

set -o nounset
set -o errexit

indir=$1
expid=$2

module load bedtools/2.27.1

cat ${indir}/${expid}.e300.qc.top1.bed \
    ${indir}/${expid}.e300.qc.top2.bed \
    ${indir}/${expid}.e300.qc.top3.bed \
    ${indir}/${expid}.e300.qc.top4.bed \
    ${indir}/${expid}.e300.qc.top5.bed \
    ${indir}/${expid}.e300.qc.top6.bed \
    ${indir}/${expid}.e300.qc.top7.bed \
    ${indir}/${expid}.e300.qc.top8.bed \
    ${indir}/${expid}.e300.qc.top9.bed \
    ${indir}/${expid}.e300.qc.top10.bed \
    | sortBed >${expid}.e300.qc.top100k.bed
