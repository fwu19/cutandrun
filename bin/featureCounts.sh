#!/usr/bin/env bash

id=$1; shift
bam=$1; shift

for bed in $(ls conp); do
    prefix="$(echo $bed | sed 's/\.bed$//')"

    cat conp/$bed | awk 'BEGIN {OFS=FS="\t"} { print $1":"$2+1"-"$3,$1,$2+1,$3,"." }' >${prefix}.saf

    featureCounts -T 6 -p -C -M --minOverlap 1 -a ${prefix}.saf -F SAF -s 0 -o ${prefix}.${id}.fragmentCounts.txt $bam
done