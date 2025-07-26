## Introduction

[This pipeline](https://github.com/fwu19/cutandrun) is to analyze sequence data from Cut and Run experiments, Cut and Tag experiments or CUTAC experiments. It was originally forked from [nf-core Cut & Run pipeline](https://nf-co.re/cutandrun/3.2.2/) and has been extensively modified in order to support the analysis of a large number of samples derived from multiple conditions and/or multiple antibodies.




### Pipeline summary

The current pipeline supports the analysis of experimental samples and the analysis of process controls in each auto-Cut & Run (CnR) or auto-Cut & Tag (CnT) experiment by running the following steps:

- Adapter trimming.

- Bowtie2 alignment to the target genome and the spike-in genome.

- Mark duplicated reads and remove duplicate reads.

- Compute genome coverage of aligned reads to the target genome and report unnormalized bedgraph and count-per-million bigwig.

- Call peaks using MACS2 with both narrow and broad modes and SEACR with stringent mode. If matching IgG is available, IgG-controlled peaks and target-only peaks will be intersected and common peaks are used for the further analysis.

For experimental samples only,

- Generate replicated peaks by group and consensus peaks by target.

- Call differential peaks between two sample groups.

- Generate MultiQC report on Trim Galore, Bowtie2, Samtools and GATK results.

- Generate an analysis report featuring a summary of read metrics, original peaks, replicated peaks, consensus peaks and differential peaks, and description of deliverable and analytical methods.

For process controls only,

- Compute precision and recall relative to multiple reference peak sets.

- Generate an updated analysis report by combining with previous data.

