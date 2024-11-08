#!/usr/bin/env Rscript

options(stringsAsFactors = F)
options(scipen = 99)
library(dplyr)
library(GenomicRanges)

args <- commandArgs(T)
id <- args[1]
peak.list <- list.files('peaks/', full.names = T)

## read original peaks ####
peaks <- lapply(
  peak.list,
  function(fname){
    if(file.size(fname) > 0){
      GenomicRanges::makeGRangesFromDataFrame(
        read.delim(fname, header = F)[1:3], seqnames.field = 'V1', start.field = 'V2', end.field = 'V3', starts.in.df.are.0based = T
      )
    }
  }); names(peaks) <- basename(peak.list)


## summarize original peaks ####
npeaks <- bind_rows(
  mapply(function(pk, fname){
    data.frame(file = fname, npeak = length(pk))
  }, peaks, names(peaks), SIMPLIFY = F
  )
) %>% 
  mutate(
    id = id,
    caller = case_when(
      grepl('seacr', basename(file)) ~ 'SEACR', 
      grepl('broad', basename(file)) ~ 'MACS2broad',
      grepl('narrow', basename(file)) ~ 'MACS2narrow',
      TRUE ~ 'others'
    ),
    toIgG = ifelse(grepl('noigg', file), "Target_only", "IgG_controlled"),
    filtered = ifelse(grepl('filtered', file),'filtered','unfiltered')
  )


## save results ####
write.table(npeaks, paste(id, 'original_peaks.csv', sep = '.'), sep = ',', quote = F, row.names = F)

