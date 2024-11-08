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
) 

## save results ####
wpeaks %>% 
  write.table(paste(id, 'original_peak_widths.csv', sep = '.'), sep = ',', quote = F, row.names = F)

