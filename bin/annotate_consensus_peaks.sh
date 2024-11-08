#!/usr/bin/env bash

genome=$1; shift
bed=$1; shift 
gtf=$1; shift

prefix="${bed%.*}"
cat $bed | awk 'BEGIN {FS=OFS="\t"} {print $1,$2,$3,$4,".","+"}' >$prefix.homer.bed

if [[ -f $gtf ]]; then
    annotatePeaks.pl $prefix.homer.bed $genome -gtf $gtf >$prefix.annotation.txt 2>$prefix.annotation.log
else
    annotatePeaks.pl $prefix.homer.bed $genome >$prefix.annotation.txt 2>$prefix.annotation.log
fi
rm $prefix.homer.bed
    
