#!/usr/bin/env Rscript

options(stringsAsFactors = F)
library(dplyr)

args <- commandArgs(T)

in_files <- file.path('multiqc_data', c('multiqc_bowtie2.txt', 'multiqc_bowtie2_1.txt', 'multiqc_picard_dups.txt'))

## metadata ####
bt2_target <- read.delim('multiqc_data/multiqc_bowtie2.txt') %>% 
  mutate(
    id = Sample,
    total_aligned = paired_aligned_one + paired_aligned_multi,
    overall_alignment_rate = overall_alignment_rate/100 
    )
bt2_spikein <- read.delim('multiqc_data/multiqc_bowtie2_1.txt') %>% 
  mutate(
    id = Sample,
    total_aligned = paired_aligned_one + paired_aligned_multi,
    overall_alignment_rate = overall_alignment_rate/100 
  )

meta <- bt2_target %>% 
  dplyr::select(-Sample) %>% 
  left_join(
      bt2_spikein %>% 
      dplyr::select(-Sample),
    by = 'id', suffix = c('_target', '_spikein')
  ) %>% 
  relocate(id)
colnames(meta)[2:ncol(meta)] <- paste('bt2', colnames(meta), sep = '_')[2:ncol(meta)]

if (file.exists('multiqc_data/multiqc_picard_dups.txt')){
  dedup <- read.delim('multiqc_data/multiqc_picard_dups.txt') %>% 
    mutate(id = Sample) %>% 
    dplyr::select(id,READ_PAIR_DUPLICATES, PERCENT_DUPLICATION, ESTIMATED_LIBRARY_SIZE) %>% 
    dplyr::rename_with(tolower) 
  colnames(dedup)[2:ncol(dedup)] <- paste('dedup', colnames(dedup), sep = '_')[2:ncol(dedup)]

  meta <- meta %>% 
    left_join(
      dedup,
      by = c('id')
    ) 
}

meta %>% 
  write.table('read_metrics.csv', sep = ',', quote = F, row.names = F)

