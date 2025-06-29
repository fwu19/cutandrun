#!/usr/bin/env Rscript

options(stringsAsFactors = F)
options(scipen = 99)
options(warn = -1)
library(dplyr)
library(GenomicRanges)


args <- commandArgs(T) # path/to/original_peaks.rds, path/to/original_peak_metrics.csv
id <- args[1]
sample_group <- args[2]
tgt <- args[3]
beds <- list.files('peaks/', full.names = T)

## functions ####
bed2gr <- function(bed){
  require(dplyr)
  require(GenomicRanges)
  if(file.size(bed) > 0){
    read.delim(bed, header = F) %>% 
      makeGRangesFromDataFrame(ignore.strand = T, seqnames.field = 'V1', start.field = 'V2', end.field = 'V3', starts.in.df.are.0based = T)
  }else{
    GRanges()
  }
}

cmp2peaks <- function(que.bed, ref, ref.name){
  require(GenomicRanges)
  que <- bed2gr(que.bed)
  if (length(que) == 0){return(data.frame())}
  
  caller <- case_when(
    grepl('narrow',que.bed) ~ 'macs2_narrow',
    grepl('broad',que.bed) ~ 'macs2_broad',
    grepl('seacr',que.bed) ~ 'seacr'
  )
  if (length(ref[[caller]]) == 0){return(data.frame())}
  
  data.frame(
    file = basename(que.bed),
    ref = ref.name,
    caller = caller,
    no_que = length(que),
    no_que_in_ref = sum(countOverlaps(que, ref[[caller]])>0),
    no_ref = length(ref[[caller]]),
    no_ref_recalled = sum(countOverlaps(ref[[caller]], que)>0)
  )
}


#### analysis ####
if (dir.exists('reference')){
  ref.rds <- list.files(file.path('reference', tgt), full.names = T)
}
if (grepl('PE25', sample_group)){
  ref.rds <- grep('PE50', ref.rds, value = T, invert = T)
}else if (sample_group %in% 'PE50_digitonin'){
  ref.rds <- grep('PE50_triton', ref.rds, value = T, invert = T)
}

recall <- list()
for (rds in ref.rds){
  ref <- readRDS(rds)
  ref.name <- gsub('\\.rds$', '', basename(rds))
  if (class(ref) == 'GRanges'){
    ref <- list(
      macs2_narrow = ref,
      macs2_broad = ref,
      seacr = ref
    )
  }
  
  for (bed in beds){
    recall[[basename(bed)]][[ref.name]] <- cmp2peaks(bed, ref, ref.name)
  }
  
}

lapply(recall, bind_rows) %>% 
  bind_rows() %>% 
  write.table(paste0(id, '.replicated_peaks.csv'), sep = ',', quote = F,row.names = F)
