#!/usr/bin/env bash

id=$1; shift
bam=$1; shift

for peak in $(ls peaks/); do
    echo -e $(basename $peak ),$( bedtools intersect -a $bam -b peaks/$peak -u | samtools view -F 256 -f 64 -c - ) 
done >${id}.reads_in_peak.csv