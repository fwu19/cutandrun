#!/usr/bin/env Rscript

options(stringsAsFactors = F)
library(dplyr)

args <- commandArgs(T)
if (length(args) < 1){
  stop("Please provide the following positional arguments, <path/to/sample_sheet.csv> <path/to/meta_table_ctrl.csv> <path/to/fragment_length>")
}

## sample sheet ####
ss <- read.csv(args[1])
mqc.dir <- paste(args[2], '04_reporting/multiqc/multiqc_data/', sep = '/')
target.genome <- args[3]
spikein.genome <- args[4]

## metadata ####
df <- read.delim(paste(mqc.dir, 'multiqc_bowtie2.txt', sep = '/')) %>% 
  mutate(
    id = gsub('.bowtie2|.spikein.bowtie2', '', Sample),
    total_aligned = paired_aligned_one + paired_aligned_multi,
    overall_alignment_rate = overall_alignment_rate/100 
    )
bt2 <- df %>% 
  filter(!grepl('spikein', Sample)) %>% 
  dplyr::select(-Sample) %>% 
  left_join(
    df %>% 
      filter(grepl('spikein', Sample)) %>% 
      dplyr::select(-Sample),
    by = 'id', suffix = c('_target', '_spikein')
  ) %>% 
  relocate(id)
colnames(bt2)[2:ncol(bt2)] <- paste('bt2', colnames(bt2), sep = '_')[2:ncol(bt2)]

dedup <- read.delim(paste(mqc.dir, 'multiqc_picard_dups.txt', sep = '/')) %>% 
  mutate(id = gsub('.target.filtered', '', Sample)) %>% 
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
flist <- list.files('./', pattern = 'fragment_length', full.names = T, recursive = T)
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
