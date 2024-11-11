#!/usr/bin/env Rscript

options(stringsAsFactors = F)
options(scipen = 99)
library(dplyr)
library(GenomicRanges)

args <- commandArgs(T) # path/to/replicated_peaks.rds, path/to/replicated_peak_metrics.csv
ss <- read.csv(args[1])
tgt <- args[2]

## parse sample sheet and keep only groups used for conp generation
ssg <- ss %>% 
  filter(target %in% tgt & call_con_peak == 'true') %>% 
  pull(group) %>% 
  unique()

## read peaks 
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
    group = gsub('.macs2.*|.seacr.*', '', file),
    caller = case_when(
      grepl('seacr', file) ~ 'SEACR',
      grepl('broad', file) ~ 'MACS2broad',
      grepl('narrow', file) ~ 'MACS2narrow',
      TRUE ~ 'others'
    )
  ) %>% 
  filter(
    group %in% ssg
  )

if (nrow(npeaks) == 0){
  quit()
}




## generate consensus peaks of each target ####
generate_conp <- function(peak.list){
  require(GenomicRanges)
  
  peak.list <- peak.list[sapply(peak.list, length)>0]
  
  if (length(peak.list) > 0){
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
    colnames(shared.conp) <- gsub('.macs2.*|.seacr.*', '', names(peak.list))
    colnames(shared.conp) <- gsub(paste0('_',tgt), '', colnames(shared.conp))
    
    ## write out conp in bed format followed by sample contribution
    df <- as.data.frame(conp)[c(1:3)] %>% 
      mutate(
        start = start - 1,
        peak.id = paste0(seqnames, ':', start+1, '-', end), # add 1 back to start
        sample.groups = apply(shared.conp, 1, function(v){paste(sort(colnames(shared.conp)[v>0]), collapse = ',')}),
        peak.count = rowSums(shared.conp)
      )
    
    return(df)
  }else{
    return(data.frame())
  }
    
}

peaks <- peaks[match(npeaks$file, names(peaks))]
conps <- lapply(
  split(peaks, npeaks$caller), 
  generate_conp
  )

## summarize conp ####
nconps <- bind_rows(mapply(
      function(pk, caller){
        data.frame(caller = caller, nconp = nrow(pk))
      }, conps, names(conps), SIMPLIFY = F
    )) %>% 
  mutate(
    target = tgt
  ) %>% 
  left_join(
    npeaks %>% 
      group_by(caller) %>% 
      summarise(
        groups = paste(group, collapse = ';')
      ) %>% 
      dplyr::select(caller, groups),
    by = 'caller'
  )


## save results ####
saveRDS(conps, paste(tgt, 'consensus_peaks.rds', sep = '.'))
write.table(nconps, paste(tgt, 'consensus_peaks.csv', sep = '.'), sep = ',', quote = F, row.names = F)


## write out consensus peaks ####
conp.files <- paste(
  tgt, 
  plyr::mapvalues(names(conps), from = c('MACS2broad', 'MACS2narrow', 'SEACR'), to = c('macs2_broad_peaks.bed', 'macs2_narrow_peaks.bed', 'seacr_peaks.bed')), 
  sep = '.'
)

mapply(
  function(x, fname){
    if (nrow(x) > 0){
      write.table(x, fname, sep = '\t', quote = F, row.names = F, col.names = F)
    }
  },
  conps,
  conp.files,
  SIMPLIFY = F
)
