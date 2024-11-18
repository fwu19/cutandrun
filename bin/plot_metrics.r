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
  'consensus_peaks',
  'consensus_beds',
  'consensus_annotation',
  'differential_peaks'
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

## conp dat ####
{
## get metrics of replicated peaks ####
if (dir.exists('replicated_peaks')){
  rep.metrics <- list.files('replicated_peaks', full.names = T, pattern = 'replicated_peaks')
  if (length(rep.metrics) > 0){
    dat$nreps <- bind_rows(lapply(rep.metrics, read.csv)) %>% 
      left_join(
        dat$meta %>% 
          dplyr::select(group,sample_group, target) %>% 
          unique.data.frame(),
        by = 'group'
      )
  }
}
## get metrics of consensus peaks ####
if (dir.exists('consensus_peaks')){
  
  conp.metrics <- list.files('consensus_peaks', full.names = T, pattern = 'consensus_peaks')
  if (length(conp.metrics) > 0){
    dat$nconps <- bind_rows(lapply(conp.metrics, read.csv))
  }
}

## get consenus peaks ####
if (dir.exists('consensus_beds')){
  conp.bed <- list.files('consensus_beds/', full.names = T, pattern = '.bed')

  ## Compare consensus peaks across sample groups ####  
  conps <- lapply(
    conp.bed, function(fname){
      read.delim(fname, header = F, col.names = c('chrom', 'start', 'end', 'conp.id','sample.groups', 'npeak'))
    }
  ); names(conps) <- basename(conp.bed)
  
  if(length(sample_groups)==1){
    dat$rep2conp <- lapply(
      conps, function(pk){
        mtx <- matrix(rep(1, length(pk)), ncol = 1)
        dimnames(mtx) <- list(pk$conp.id, sample_groups)
      }
    )
  }else{
    dat$rep2conp <- lapply(
      conps, function(pk){
        mtx <- t(sapply(
          strsplit(pk$sample.groups, split = ','),
          function(v){
            as.integer(sample_groups %in% v)
          }
        ))
        dimnames(mtx) <- list(pk$conp.id, sample_groups)
        return(mtx)
      }
    )
    
  }
  
  ## Compare MACS2 and SEACR by consensus peaks ####
  conps.gr <- lapply(
    conps, function(x){
      if(nrow(x)>0){
        makeGRangesFromDataFrame(x, ignore.strand = T, seqnames.field = 'chrom', start.field = 'start', end.field = 'end', starts.in.df.are.0based = T)
      }
    }
  )
  
  dat$seacr2macs <- lapply(
    split(names(conps.gr), gsub('.macs2_.*|.seacr_.*', '', names(conps.gr))),
    function(i){
      require(GenomicRanges)
      peak.list <- conps.gr[i]
      
      if(length(peak.list) == 1){return(NULL)} # only one peak set is available
      
      peak.list[sapply(peak.list, length) == 0] <- NULL
      
      if(length(peak.list) < 2){return(NULL)} # less than two peak sets are available
      
      tgt <- gsub('.macs2_.*|.seacr_.*', '', i[1])
      
      conp <- GRanges()
      for (pk in peak.list){
        conp <- c(conp, pk)
      }
      conp <- reduce(conp)
      
      overlap2conp <- as.data.frame(sapply(
        peak.list, function(pk){countOverlaps(conp, pk)>0}
      ))
      colnames(overlap2conp) <- gsub('_peaks', '', stringr::str_extract(colnames(overlap2conp), "macs.*peaks|seacr.*peaks"))
      
      return(overlap2conp)
    })
  
}

## genomic distribution of consensus peaks ####
if (dir.exists('consensus_annotation')){
  ann.txt <- list.files('consensus_annotation', pattern = 'annotation.txt', full.names = T, recursive = T)
  
  
  dat$conp_ann <- lapply(ann.txt, function(fname){
    read.delim(fname, header = T) %>%   
      subset(!is.na(Annotation)) %>% 
      mutate(
        conp.id = .[,1],
        genomic.location = gsub(' .*', '', Annotation)
      ) %>% 
      dplyr::select(conp.id, genomic.location)
    
  }); names(dat$conp_ann) <- gsub('.annotation.txt','',basename(ann.txt))
  
}


}

## Save results ####
saveRDS(dat, 'data.rds')

## cat dp.rds ####
if (dir.exists('differential_peaks')){
  dp <- do.call(c, lapply(
    list.files('differential_peaks', full.names = T), readRDS
  ))
  dp %>% saveRDS('dp.rds')
}

## prepare report.Rmd ####
file.copy('report/report.setup.Rmd', 'report.Rmd')
if (file.exists('read_metrics.csv')){file.append('report.Rmd', 'report/report.read_qc.Rmd')}
if (dir.exists('original_peaks')){file.append('report.Rmd', 'report/report.orig_qc.Rmd')}
if (dir.exists('consensus_peaks')){file.append('report.Rmd', 'report/report.conp_qc.Rmd')}
if (dir.exists('differential_peaks')){
  file.append('report.Rmd', 'report/report.diff_peaks.Rmd')
  }
file.append('report.Rmd', 'report/report.deliverables.Rmd')
file.append('report.Rmd', 'report/report.methods.Rmd')
