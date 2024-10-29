#!/usr/bin/env bash

# https://www.protocols.io/view/cut-amp-tag-data-processing-and-analysis-tutorial-e6nvw93x7gmk/v1

## Load required modules
# module load SAMtools/1.11-GCC-10.2.0

sampleId=$1; shift
inBam=$1; shift
if [[ -f $inBam ]]; then
    samtools view -F 0x04 $inBam | awk -F'\t' 'function abs(x){return ((x < 0.0) ? -x : x)} {print abs($9)}' | sort | uniq -c | awk -v OFS="\t" '{print $2, $1/2}' >$sampleId.fragment_length.txt
fi
