#!/usr/bin/env Rscript

## Generate replicated peaks present in `min.reps` replicates
## Output singleton peaks if only 1 replicate is available.

options(stringsAsFactors = F)
options(scipen = 99)
options(warn = -1)
library(dplyr)
library(GenomicRanges)

args <- commandArgs(T) # path/to/original_peaks.rds, path/to/original_peak_metrics.csv
id <- args[1]
sample_group <- args[2]
tgt <- args[3]
peaks.txt <- list.files('peaks/', full.names = T)
hiconp.rds <- file.path('high_confidence_peaks', paste0(sample_group, '_', tgt, '.rds'))

if (file.exists(hiconp.rds)){
  reps <- list()
  hiconp <- readRDS(hiconp.rds)
  for (pkt in peaks.txt){
    peak <- read.delim(pkt, header = F) %>% 
      makeGRangesFromDataFrame(ignore.strand = T, seqnames.field = 'V1', start.field = 'V2', end.field = 'V3', starts.in.df.are.0based = T)
    
    if (grepl('narrow',pkt) & 'macs2narow' %in% names(hiconp)){
      reps$macs2narrow <- data.frame(
        file = basename(pkt),
        replicated_peaks = sum(countOverlaps(hiconp$macs2narrow, peak)>0),
        hiconf_peaks = length(hiconp$macs2narrow),
        caller = "MACS2narrow"
      )
    }
    
    if (grepl('broad',pkt) & 'macs2broad' %in% names(hiconp)){
      reps$macs2broad <- data.frame(
        file = basename(pkt),
        replicated_peaks = sum(countOverlaps(hiconp$macs2broad, peak)>0),
        hiconf_peaks = length(hiconp$macs2broad),
        caller = "MACS2broad"
      )
    }

    if (grepl('seacr',pkt) & 'seacr' %in% names(hiconp)){
      reps$seacr <- data.frame(
        file = basename(pkt),
        replicated_peaks = sum(countOverlaps(hiconp$seacr, peak)>0),
        hiconf_peaks = length(hiconp$seacr),
        caller = "SEACR"
      )
    }
    
    rm(peak)
    
  }
  
  bind_rows(reps) %>% 
    mutate(
      id = id
    ) %>% 
    relocate(id) %>% 
    write.table(paste0(id, '.replicated_peaks.csv'), sep = ',', quote = F,row.names = F)
}else{
  system(paste0('touch ', id, '.replicated_peaks.csv'))
}