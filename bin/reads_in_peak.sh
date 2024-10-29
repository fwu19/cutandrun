#!/usr/bin/env bash

bam=$1; shift
peak=$1; shift
echo -e $(basename $peak ),$( bedtools intersect -a $bam -b $peak -u | samtools view -F 256 -f 64 -c - ) >$(basename $peak ).reads_in_peak.csv # count only R1 reads
