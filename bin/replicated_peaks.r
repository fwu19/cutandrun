#!/usr/bin/env Rscript

## Generate replicated peaks present in `min.reps` replicates
## Output singleton peaks if only 1 replicate is available.

options(stringsAsFactors = F)
options(scipen = 99)
options(warn = -1)
library(dplyr)
library(GenomicRanges)

args <- commandArgs(T) # path/to/original_peaks.rds, path/to/original_peak_metrics.csv
group.id <- args[1] # will be used as output prefix
min.reps <- as.integer(args[2])

merge_reps <- function(peak.list, min.reps = 2){
  require(GenomicRanges)
  
  k <- sapply(peak.list,length)>0 # at least one rep has peaks called

  if(sum(k) == 0){return(NULL)}
  
  peak.list[!k] <- NULL
  lst <- list()
  if (length(peak.list) == 1 ){ # only 1 replicate
    lst$replicated <- GRanges()
    lst$singleton <- peak.list[[1]]
    
  }else { # more than 2 replicates
    conp <- GRanges()
    
    for (pk in peak.list){
      conp <- c(conp, pk)
    }
    
    conp <- reduce(conp)
    
    keep.conp <- rowSums(as.matrix(sapply(peak.list, function(pk){countOverlaps(conp,pk)}))>0) >= min.reps # peaks shared by at least this number of replicates
    
    if(sum(keep.conp) > 0){
      lst$replicated <- conp[which(keep.conp)]
    }else{
      lst$replicated <- GRanges()
    }
    lst$singleton <- GRanges()
  }
  
  return(lst)
}


## read the final set of peaks (filtered peaks or target-only peaks if control is not available)
peak.list <- list.files('peaks/', full.names = T)
peaks <- lapply(
  peak.list,
  function(fname){
    if(file.size(fname) > 1){
      GenomicRanges::makeGRangesFromDataFrame(
        read.delim(fname, header = F)[1:3], seqnames.field = 'V1', start.field = 'V2', end.field = 'V3', starts.in.df.are.0based = T
      )
    }else{
      GRanges()
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
      function(lst, caller){
        if (is.null(lst)){
          data.frame(
            caller = caller, 
            peak.count = 0
          )
          
        }else {
          data.frame(
            caller = caller, 
            peak.count = length(lst$replicated)
          )
        }

      }, 
      reps, names(reps), SIMPLIFY = F
    )) %>% 
  mutate(
    group = group.id,
    output_file = paste(
      group.id, 
      plyr::mapvalues(caller, from = c('MACS2broad', 'MACS2narrow', 'SEACR'), to = c('macs2_broad_peaks.bed', 'macs2_narrow_peaks.bed', 'seacr_peaks.bed')), 
      sep = '.')
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
  ) %>% 
  mutate(
   output_file = paste(ifelse(peak.count > 0, 'multiple_replicates', 'single_replicate'), output_file, sep = '/') 
  )


## save results ####
saveRDS(reps, paste(group.id, 'replicated_peaks.rds', sep = '.'))
write.table(nreps, paste(group.id, 'replicated_peaks.csv', sep = '.'), sep = ',', quote = F, row.names = F)


## write out replicated peaks ####
reps[sapply(reps, is.null)] <- NULL

rep.files <- paste(
      group.id, 
      plyr::mapvalues(names(reps), from = c('MACS2broad', 'MACS2narrow', 'SEACR'), to = c('macs2_broad_peaks.bed', 'macs2_narrow_peaks.bed', 'seacr_peaks.bed')), 
      sep = '.'
    )

## replicated peaks
od <- 'multiple_replicates'
if (!dir.exists(od)){dir.create(od, recursive = T)}
mapply(
  function(lst, fname){
    if (length(lst$replicated) > 0){
      df <- as.data.frame(lst$replicated)[1:3]
      df[,2] <- df[,2] - 1
      write.table(df, file.path(od, fname), sep = '\t', quote = F, row.names = F, col.names = F)
    }
    
  }, 
  reps, 
  rep.files, 
  SIMPLIFY = F
)       

## singleton
od <- 'single_replicate'
if (!dir.exists(od)){dir.create(od, recursive = T)}
mapply(
  function(lst, fname){
    if (length(lst$singleton) > 0){
      df <- as.data.frame(lst$singleton)[1:3]
      df[,2] <- df[,2] - 1
      write.table(df, file.path(od, fname), sep = '\t', quote = F, row.names = F, col.names = F)
    }
    
  }, 
  reps, 
  rep.files, 
  SIMPLIFY = F
)       
