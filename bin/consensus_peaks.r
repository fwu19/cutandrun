#!/usr/bin/env Rscript

options(stringsAsFactors = F)
options(scipen = 99)
library(dplyr)
library(GenomicRanges)

args <- commandArgs(T) # path/to/replicated_peaks.rds, path/to/replicated_peak_metrics.csv

reps <- readRDS('replicated_peaks.rds')
ms <- read.csv('replicated_peak_metrics.csv')
targets <- sort(setdiff(unique(ms$target), 'IgG')) 

## generate consensus peaks of each target ####
generate_conp <- function(peak.list){
  require(GenomicRanges)
  
  peak.list <- peak.list[sapply(peak.list, length)>0]
  
  if(length(peak.list) > 1){
  ## merge peaks
    conp <- GRanges()
    for (pk in peak.list){
      conp <- c(conp, pk)
    }
    conp <- reduce(conp)
    
    ## identify sample_group that contribute to the conp
    shared.conp <- as.matrix(
      sapply(peak.list, function(pk){countOverlaps(conp,pk)})
    )
    colnames(shared.conp) <- sapply(strsplit(names(peak.list), split = ':'), function(v){v[2]})

    ## write out conp in bed format followed by sample contribution
    df <- as.data.frame(conp)[c(1:3)] %>% 
      mutate(
        start = start - 1,
        peak.id = paste0(seqnames, ':', start+1, '-', end), # add 1 back to start
        sample.groups = apply(shared.conp, 1, function(v){paste(sort(colnames(shared.conp)[v>0]), collapse = ',')}),
        peak.count = rowSums(shared.conp)
      )
    
    return(df)
  }else if (length(peak.list)==1){
    return(
      as.data.frame(peak.list[[1]])[c(1:3)] %>% 
        mutate(
          start = start - 1,
          peak.id = paste0(seqnames, ':', start+1, '-', end), # add 1 back to start
          sample.groups = strsplit(names(peak.list)[1], split = ':')[[1]][2],
          peak.count = 1
        )
      )
  }else{
    return(data.frame())
  }
    
}

conps <- lapply(
  split(reps[ms$rep.id], paste(ms$target, ms$caller, sep = ':')), 
  generate_conp
  )


out.bed <- paste0(
  gsub(':.*','.',names(conps)),
  dplyr::case_when(
    grepl('MACS2narrow', names(conps)) ~ 'narrow_peaks.bed',
    grepl('MACS2broad', names(conps)) ~ 'broad_peaks.bed',
    grepl('SEACR', names(conps)) ~ 'stringent.bed'
  )
)

mapply(
  function(x, fname){
    if(is.null(x) | nrow(x) == 0){
      system(paste('touch', fname))
    }else if (nrow(x) > 0){
      write.table(x, fname, sep = '\t', quote = F, row.names = F, col.names = F)
    }else{
      system(paste('touch', fname))
    }
  }, 
  conps, 
  out.bed,
  SIMPLIFY = F
)


## summarize conp ####
nconps <- as.data.frame(bind_rows(mapply(
  function(pk, id){
    data.frame(conp.id = id, nconp = nrow(pk))
  }, conps, names(conps), SIMPLIFY = F
)))
nconps[,c( 'target', 'caller')] <- do.call(rbind,strsplit(nconps$conp.id, split = ':'))


## save results ####
saveRDS(conps, 'consensus_peaks.rds')
write.table(nconps, 'consensus_peak_metrics.csv', sep = ',', quote = F, row.names = F)


