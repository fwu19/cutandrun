#!/usr/bin/env Rscript

options(stringsAsFactors = F)
library(dplyr)

args <- commandArgs(T)
if (length(args) < 1){
  stop("Missing required arguments <path/to/sample_sheet.csv> ")
}

## sample sheet ####
ss <- read.csv(args[1])

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

bt2 <- bt2_target %>% 
  dplyr::select(-Sample) %>% 
  left_join(
      bt2_spikein %>% 
      dplyr::select(-Sample),
    by = 'id', suffix = c('_target', '_spikein')
  ) %>% 
  relocate(id)
colnames(bt2)[2:ncol(bt2)] <- paste('bt2', colnames(bt2), sep = '_')[2:ncol(bt2)]

dedup <- read.delim('multiqc_data/multiqc_picard_dups.txt') %>% 
  mutate(id = Sample) %>% 
  dplyr::select(id,READ_PAIR_DUPLICATES, PERCENT_DUPLICATION, ESTIMATED_LIBRARY_SIZE) %>% 
  dplyr::rename_with(tolower) 
colnames(dedup)[2:ncol(dedup)] <- paste('dedup', colnames(dedup), sep = '_')[2:ncol(dedup)]

meta <- bt2 %>% 
  left_join(
    dedup,
    by = c('id')
  ) %>% 
  left_join(
    ss %>% dplyr::select(id, group, replicate, control_id, sample_group, sample_replicate, target), 
    by = c('id')
)
write.table(meta, 'read_metrics.csv', sep = ',', quote = F, row.names = F)


## fragment lengths ####
flist <- list.files('fragment_length/', pattern = 'fragment_length', full.names = T, recursive = T)
bind_rows(lapply(
  flist,
  function(fname){
    if(file.size(fname) > 0 ){
      read.delim(fname, header = F, col.names = c('length','count'), colClasses = 'numeric') %>% 
        mutate(
          weight = count/sum(count), 
          id = gsub('.fragment_length.txt','',basename(fname)))
    }
  })) %>% 
  left_join(
    meta %>% select(id, sample_group, sample_replicate, target),
    by = 'id'
  ) %>% 
  saveRDS('fragment_length.rds')
