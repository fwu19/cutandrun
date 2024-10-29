#!/usr/bin/env Rscript

options(stringsAsFactors = F)
options(scipen = 99)
library(dplyr)
library(GenomicRanges)

args <- commandArgs(T)

## sample sheet ####
meta <- read.csv('read_metrics.csv')
targets <- sort(setdiff(unique(meta$target), 'IgG')) # targets to make QC plots

peak.list <- list.files('peaks/', full.names = T)
rip.list <- list.files('rip/', full.names = T)

# write.table(as.matrix(peak.list), 'peaks.txt', quote = F, row.names = F, col.names = F)
# write.table(as.matrix(rip.list), 'rip.txt', quote = F, row.names = F, col.names = F)


## read original peaks ####
peaks <- lapply(
  peak.list,
  function(fname){
    if(file.size(fname) > 0){
      GenomicRanges::makeGRangesFromDataFrame(
        read.delim(fname, header = F)[1:3], seqnames.field = 'V1', start.field = 'V2', end.field = 'V3', starts.in.df.are.0based = T
      )
    }
  }); names(peaks) <- peak.list


## summarize original peaks ####
npeaks <- bind_rows(
  mapply(function(pk, fname){
    data.frame(file = fname, npeak = length(pk))
  }, peaks, names(peaks), SIMPLIFY = F
  )
) %>% 
  mutate(id = gsub('.macs.*|.seacr.*','',basename(file))) %>% 
  left_join(
    meta %>% select(id, sample_group, target, sample_replicate, bt2_total_aligned_target), 
    by = c("id")
  ) %>% 
  mutate(
    caller = case_when(
      grepl('seacr', basename(file)) ~ 'SEACR', 
      grepl('broad', basename(file)) ~ 'MACS2broad',
      grepl('narrow', basename(file)) ~ 'MACS2narrow',
      TRUE ~ 'others'
    ),
    toIgG = ifelse(grepl('noigg', file), "Target_only", "IgG_controlled"),
    filtered = ifelse(grepl('filtered', file),'filtered','unfiltered')
  ) %>% 
  left_join(
    bind_rows(lapply(
      rip.list, 
      function(fname){
        if(file.size(fname) > 0){
          read.csv(fname, header = F, col.names = c('file', 'reads_in_peak'))
        }else{
          return(NULL)
        }
      }
    )) ,
    by = c('file' = 'file')
  ) %>% 
  mutate(frip = reads_in_peak/bt2_total_aligned_target)

## compute peak widths ####
wpeaks <- bind_rows(mapply(
  function(pk, file){
    if(length(pk)>0){
      tab <- table(width(pk))
      data.frame(
        file = file,
        length = as.integer(names(tab)),
        count = as.vector(tab)
      ) %>% 
        mutate(weight = count/ sum(count))
      
    }
  }, peaks, names(peaks), SIMPLIFY = F)
) %>% 
  left_join(
    npeaks %>% select(file, id, sample_group, target, sample_replicate, caller),
    by = c('file' = 'file')
  )

## save results ####
saveRDS(peaks, 'original_peaks.rds')
write.table(npeaks, 'original_peak_metrics.csv', sep = ',', quote = F, row.names = F)
saveRDS(wpeaks, 'original_peak_widths.rds')

