#!/usr/bin/env Rscript

options(stringsAsFactors = F)
library(dplyr)
library(patchwork)
library(ggplot2)
library(GenomicRanges)

args <- as.vector(commandArgs(T))

input.files <- c('sample_sheet.csv', 'read_metrics.csv')
if (sum(file.exists(input.files))<2){
  stop ('Provide sample_sheet.csv and read_metrics.csv!')
}

input.dirs <- c(
  'fragment_lengths',
  'original_peaks',
  'original_peak_widths',
  'reads_in_peak',
  'replicated_peaks',
  'saved_data'
)
  
dat <- list()

## read sample sheet ####
ss <- read.csv('sample_sheet.csv')
targets <- sort(setdiff(unique(ss$target), c('IgG', 'Input', 'input'))) # targets to make QC plots
sample_groups <- sort(unique(ss$sample_group))
dat$ss <- ss

## reads dat ####
{
## get read metrics ####
if (file.exists('read_metrics.csv')){
  dat$meta <- read.csv('read_metrics.csv') %>% 
    left_join(
      ss %>% dplyr::select(id, group, sample_group, target, sample_replicate), 
      by = 'id'
    )
  meta <- dat$meta # for plotting
  
}

## get fragment lengths ####
if (dir.exists('fragment_lengths')){
  length.list <- list.files('fragment_lengths/', pattern = 'fragment_lengths', full.names = T, recursive = T)
  
  dat$frag_lens <- bind_rows(lapply(
    length.list,
    function(fname){
      if(file.size(fname) > 1 ){
        read.delim(fname, header = F, col.names = c('length','count'), colClasses = 'numeric') %>%
          mutate(
            weight = count/sum(count),
            id = gsub('.fragment_lengths.txt','',basename(fname)))
      }
    })) %>% 
    left_join(
      dat$meta %>% dplyr::select(id,sample_group, target, sample_replicate),
      by = 'id'
    )
  
}
  
  
}

## peaks dat ####
{
## get metrics of the entire set of original peaks ####
if (dir.exists('original_peaks')){
  peak.metrics <- list.files('original_peaks', full.names = T, pattern = 'original_peaks')
  if (length(peak.metrics) > 0){
    dat$npeaks <- bind_rows(lapply(
      peak.metrics, 
      function(fname){
        if(file.size(fname) > 1){ 
          read.csv(fname)
        }else{
            return(NULL)
        }
      }
      )) %>% 
      left_join(
        dat$meta %>% dplyr::select(id,sample_group, target, sample_replicate),
        by = 'id'
      )
  }
  
  
}

## get metrics of the final set of original peak widths  ####
if (dir.exists('original_peak_widths')){
  peak.widths <- list.files('original_peak_widths', full.names = T, pattern = 'original_peak_widths')
  if (length(peak.widths) > 0){
    dat$wpeaks <- bind_rows(lapply(
      peak.widths,       
      function(fname){
        if(file.size(fname) > 1){ 
          read.csv(fname)
        }else{
          return(NULL)
        }
      }
      
      )) %>% 
      left_join(
        dat$npeaks %>% dplyr::select(file,id, sample_group, target, sample_replicate,caller),
        by = 'file'
      )
  }
  
  
}
## get reads in peak ####
if (dir.exists('reads_in_peak')){
  peak.reads <- list.files('reads_in_peak', full.names = T, pattern = 'reads_in_peak')
  dat$frip <- bind_rows(lapply(
    peak.reads,
    function(fname){
      if(file.size(fname) > 1 ){
        read.csv(fname, header = F, col.names = c('file','reads_in_peak')) %>%
          mutate(
            weight = reads_in_peak/sum(reads_in_peak),
            id = gsub('.reads_in_peak.csv','',basename(fname)))
      }
    })) %>% 
    left_join(
      dat$meta %>% dplyr::select(id,sample_group, target, sample_replicate,bt2_total_aligned_target),
      by = 'id'
    ) %>% 
    mutate(
      FRiP = reads_in_peak/bt2_total_aligned_target,
      caller = case_when(
        grepl('seacr', basename(file)) ~ 'SEACR', 
        grepl('broad', basename(file)) ~ 'MACS2broad',
        grepl('narrow', basename(file)) ~ 'MACS2narrow',
        TRUE ~ 'others'
      ),
      filtered = ifelse(grepl('filtered', file),'filtered','unfiltered')
    )
  
}


}

## get metrics of replicated peaks ####
if (dir.exists('replicated_peaks')){
  rep.metrics <- list.files('replicated_peaks', full.names = T, pattern = 'replicated_peaks')
  if (length(rep.metrics) > 0){
    k <- file.size(rep.metrics) > 0
    if (sum(k) > 0){
      dat$nreps <- bind_rows(lapply(rep.metrics[k], read.csv)) %>% 
        left_join(
          dat$npeaks %>% 
            dplyr::select(file, sample_group, target, sample_replicate) %>% 
            unique.data.frame(),
          by = 'file'
        )
    }
  }
}


## Save results and add to saved data if available ####
add_new <- function(df0, df1){
  if ('file' %in% colnames(df1)){
    bind_rows(
      df0 %>% filter(!file %in% df1$file),
      df1
    ) %>% 
      unique.data.frame()
    
  }else if ('id' %in% colnames(df1)){
    bind_rows(
      df0 %>% filter(!id %in% df1$file),
      df1
    ) %>% 
      unique.data.frame()
    
  }else{
    bind_rows(df0, df1) %>% 
      unique.data.frame()
  }
}

if (file.exists("saved_data/data.rds")){
  dat.old <- readRDS("saved_data/data.rds")
  
  common.names <- intersect(names(dat), names(dat.old))
  dat[common.names] <- mapply(
    add_new, dat.old[common.names], dat[common.names], SIMPLIFY = F
  )
  
  dat[setdiff(names(dat.old), names(dat))] <- dat.old[setdiff(names(dat.old), names(dat))]
  
}
saveRDS(dat, 'data.rds')
