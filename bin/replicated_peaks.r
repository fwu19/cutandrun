#!/usr/bin/env Rscript

options(stringsAsFactors = F)
options(scipen = 99)
library(dplyr)
library(GenomicRanges)

args <- commandArgs(T) # path/to/original_peaks.rds, path/to/original_peak_metrics.csv
ss <- read.csv(args[1])
min.reps <- as.integer(args[2])
peaks <- readRDS('original_peaks.rds')
npeaks <- read.csv('original_peak_metrics.csv')
targets <- sort(setdiff(unique(npeaks$target), 'IgG')) 

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

## exclude some samples if required
if ('exclude_rep_peaks' %in% colnames(ss)){
  npeaks <- npeaks %>% 
    left_join(
      ss %>% dplyr::select(sample_id, exclude_rep_peaks),
      by = 'sample_id'
    )
}else{
  npeaks$exclude_rep_peaks <- 'false'
}
npeaks <- npeaks %>% 
  filter(exclude_rep_peaks == 'false')

if ('filtered' %in% colnames(npeaks) & sum(npeaks$filtered %in% 'filtered') > 0){
  k <- npeaks$filtered %in% 'filtered'
}else{
  k <- 1:nrow(npeaks)
}
reps <- lapply(
  split(peaks[npeaks$file[k]], paste(npeaks$caller, npeaks$sample_group, npeaks$target, sep = ':')[k]),
  merge_reps, min.rep = min.reps
)

## write out replicated peaks
rep.files <- sapply(
  strsplit(names(reps), split = ":"),
  function(v){
    paste0(
      v[2], '_', v[3], '.', 
      dplyr::case_when(
        grepl('narrow', v[1]) ~ 'narrow_peaks.bed', 
        grepl('broad', v[1]) ~ 'broad_peaks.bed',
        TRUE ~ 'stringent.bed')
    )
  }
)

mapply(
  function(pk, fname){
    if (length(pk) > 0){
      df <- as.data.frame(pk)[1:3]
      df[,2] <- df[,2] - 1
      
      write.table(df, fname, sep = '\t', quote = F, row.names = F, col.names = F)
    }else{
      system(paste('touch', fname))
    }
  }, reps, rep.files, SIMPLIFY = F
)       


## summarize reproduced peaks ####
# nreps <- as.data.frame(bind_rows(mapply(
#   function(pk, id){
#     data.frame(rep.id = id, nrep = length(pk))
#   }, reps, names(reps), SIMPLIFY = F
# )))
# nreps[,c( 'caller', 'sample_group', 'target')] <- do.call(rbind,strsplit(nreps$rep.id, split = ':'))




df <- as.data.frame(data.table::rbindlist(mapply(
  function(pk, id){
    data.frame(rep.id = id, nrep = length(pk))
  }, reps, names(reps), SIMPLIFY = F
), use.names = T, fill = T))
df[,c( 'caller', 'sample_group', 'target')] <- do.call(rbind, strsplit(df$rep.id, split = ':'))
nreps <- npeaks %>% 
  filter(filtered == 'filtered') %>% 
  mutate(sample_replicate = paste0('rep', sample_replicate, '.peaks')) %>% 
  dplyr::select(sample_group, target, caller, sample_replicate, npeak) %>% 
  tidyr::pivot_wider(names_from = sample_replicate, values_from = npeak) 
nreps$valid.replicates <- rowSums(nreps[, grep('^rep.*peaks$', colnames(nreps))]>0, na.rm = T)

nreps <- nreps %>%
  left_join(
    df %>% dplyr::rename("replicated.peaks" = "nrep"),
    by = c('sample_group', 'target', 'caller')
  ) %>% 
  mutate(
    replicated.peaks = ifelse(valid.replicates < 2, NA, replicated.peaks)
  ) %>%
  arrange(sample_group, target, caller) %>% 
  dplyr::select(-valid.replicates) %>% 
  relocate(rep.id, .after = last_col())





## save results ####
saveRDS(reps, 'replicated_peaks.rds')
write.table(nreps, 'replicated_peak_metrics.csv', sep = ',', quote = F, row.names = F)



