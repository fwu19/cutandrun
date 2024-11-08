#!/usr/bin/env Rscript

options(stringsAsFactors = F)
options(scipen = 99)
library(dplyr)
library(GenomicRanges)

args <- commandArgs(T) # path/to/original_peaks.rds, path/to/original_peak_metrics.csv
group.id <- args[1] # will be used as output prefix
min.reps <- as.integer(args[2])

merge_reps <- function(rep.peaks, min.reps){
  require(GenomicRanges)
  
  peaks2merge <- sapply(rep.peaks,length)>0 # at least one rep has peaks called
  
  if(sum(peaks2merge)>1){
    conp <- GRanges()
    
    for (pk in rep.peaks[peaks2merge]){
      conp <- c(conp, pk)
    }
    
    conp <- reduce(conp)
    
    keep.conp <- rowSums(
      as.matrix(sapply(rep.peaks, function(pk){
        if(length(pk)>0){countOverlaps(conp,pk)}else{rep(0,length(conp))}
      }))
      >0) >= min.reps # peaks shared by at least this number of replicates
    
    if(sum(keep.conp) > 0){
      return(conp[keep.conp])
    }
  }else if (sum(peaks2merge) == 1 ){
    return(rep.peaks[[which(peaks2merge)]])
  }else{
    return(NULL)
  }
}


## read the final set of peaks (filtered peaks or target-only peaks if control is not available)
peak.list <- list.files('peaks/', full.names = T)
peaks <- lapply(
  peak.list,
  function(fname){
    if(file.size(fname) > 0){
      GenomicRanges::makeGRangesFromDataFrame(
        read.delim(fname, header = F)[1:3], seqnames.field = 'V1', start.field = 'V2', end.field = 'V3', starts.in.df.are.0based = T
      )
    }
  }); names(peaks) <- basename(peak.list)

## exclude some samples if required
npeaks <- data.frame(
  file = names(peaks)
) %>% 
  mutate(
    id = gsub('.macs2.*|.seacr.*', '', file),
    caller = case_when(
      grepl('seacr', file) ~ 'SEACR',
      grepl('broad', file) ~ 'MACS2broad',
      grepl('narrow', file) ~ 'MACS2narrow',
      TRUE ~ 'others'
    )
  ) 

reps <- lapply(
  split(peaks, npeaks$caller),
  merge_reps, min.rep = min.reps
)

## summarize reproduced peaks ####
nreps <- bind_rows(mapply(
      function(pk, caller){
        data.frame(caller = caller, nrep = length(pk))
      }, reps, names(reps), SIMPLIFY = F
    )) %>% 
  mutate(
    group = group.id
  ) %>% 
  left_join(
    npeaks %>% 
      group_by(caller) %>% 
      summarise(
        caller = unique(caller),
        ids = paste(id, collapse = ';')
      ) %>% 
      dplyr::select(caller, ids),
    by = 'caller'
  )


## save results ####
saveRDS(reps, paste(group.id, 'replicated_peaks.rds', sep = '.'))
write.table(nreps, paste(group.id, 'replicated_peaks.csv', sep = '.'), sep = ',', quote = F, row.names = F)


## write out replicated peaks ####
rep.files <- paste(
      group.id, 
      plyr::mapvalues(names(reps), from = c('MACS2broad', 'MACS2narrow', 'SEACR'), to = c('macs2_broad_peaks.bed', 'macs2_narrow_peaks.bed', 'seacr_peaks.bed')), 
      sep = '.'
    )


mapply(
  function(pk, fname){
    if (length(pk) > 0){
      df <- as.data.frame(pk)[1:3]
      df[,2] <- df[,2] - 1
      write.table(df, fname, sep = '\t', quote = F, row.names = F, col.names = F)
    }
  }, 
  reps, 
  rep.files, 
  SIMPLIFY = F
)       

