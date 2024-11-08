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
  'consensus_annotation'
)

dat <- list()
funcs <- list()
figs <- list()

## read sample sheet ####
ss <- read.csv('sample_sheet.csv')
targets <- sort(setdiff(unique(ss$target), c('IgG', 'Input', 'input'))) # targets to make QC plots

## get read metrics ####
if (file.exists('read_metrics.csv')){
 dat$meta <- read.csv('read_metrics.csv') %>% 
   left_join(
     ss, by = 'id'
   )
 meta <- dat$meta # for plotting
}

## plot accordingly
{
  ## input reads ####
  qc <- 'seq_depth'
  
  funcs[[qc]] <- function(df, tgt, var.x = 'bt2_total_reads_target', var.y = 'sample_group', xlab = 'Total reads (in million)', ylab = '', color = 'Replicate', plot.title = 'Sequenced reads'){
    require(ggplot2)
    
    as.data.frame(df) %>% 
      mutate(sample_replicate = ifelse(is.na(sample_replicate), 1, sample_replicate)) %>% 
      mutate(x = df[,var.x]/1e6, y = df[,var.y], color = factor(sample_replicate)) %>% 
      subset(target %in% tgt) %>% 
      ggplot(mapping = aes(x = x, y = y, color = color))+
      geom_jitter(height = 0.1, width = 0)+
      labs(x = xlab, y = ylab, color = color, title = plot.title)+
      theme_bw()+
      theme(
        text = element_text(size = 8),
        axis.text.x = element_text(angle = 90, hjust = 1),
        strip.text.y = element_text(angle = 0),
        legend.position = 'top'
      )
    
  }
  
  figs[[qc]] <- lapply(c(targets,'IgG'), funcs[[qc]], df = meta)
  names(figs[[qc]]) <- c(targets,'IgG')
  
  
  ## Aligned reads to the target genome ####
  qc <- 'aligned_reads'
  
  funcs[[qc]] <- function(df, tgt, var.x = 'bt2_total_aligned_target', var.y = 'sample_group', xlab = 'Total aligned reads (in million) to the target genome', ylab = '', color = 'Replicate', plot.title = 'Reads aligned to target genome'){
    require(ggplot2)
    
    as.data.frame(df) %>% 
      mutate(sample_replicate = ifelse(is.na(sample_replicate), 1, sample_replicate)) %>% 
      mutate(x = df[,var.x]/1e6, y = df[,var.y], color = factor(sample_replicate)) %>% 
      subset(target %in% tgt) %>%
      ggplot(mapping = aes(x = x, y = y, color = color))+
      geom_jitter(height = 0.1, width = 0)+
      labs(
        x = xlab, y = ylab, 
        color = color, title = plot.title)+
      theme_bw()+
      theme(
        text = element_text(size = 8),
        axis.text.x = element_text(angle = 90, hjust = 1),
        strip.text.y = element_text(angle = 0),
        legend.position = 'top'
      )
    
  }
  
  figs[[qc]] <- lapply(c(targets,'IgG'), funcs[[qc]], df = meta)
  names(figs[[qc]]) <- c(targets,'IgG')
  
  
  ## Alignment rate to the target genome ####
  qc <- 'aligned_pct' 
  funcs[[qc]] <- function(df, tgt, var.x = 'bt2_overall_alignment_rate_target', var.y = 'sample_group', xlab = 'Fraction of reads aligned to the target genome', ylab = '', color = 'Replicate', plot.title = 'Alignment rate to target genome'){
    require(ggplot2)
    
    as.data.frame(df) %>% 
      mutate(sample_replicate = ifelse(is.na(sample_replicate), 1, sample_replicate)) %>% 
      mutate(x = df[,var.x], y = df[,var.y], color = factor(sample_replicate)) %>% 
      subset(target %in% tgt) %>%
      ggplot(mapping = aes(x = x, y = y, color = color))+
      geom_jitter(height = 0.1, width = 0)+
      labs(
        x = xlab, y = ylab, 
        color = color, title = plot.title)+
      theme_bw()+
      theme(
        text = element_text(size = 8),
        axis.text.x = element_text(angle = 90, hjust = 1),
        strip.text.y = element_text(angle = 0),
        legend.position = 'top'
      )
    
  }
  
  figs[[qc]] <- lapply(c(targets,'IgG'), funcs[[qc]], df = meta)
  names(figs[[qc]]) <- c(targets,'IgG')
  
  
  ## Aligned reads to spike-in genome ####
  qc <- 'aligned_reads_spikein'
  
  funcs[[qc]] <- function(df, tgt, var.x = 'bt2_overall_alignment_rate_spikein', var.y = 'sample_group', xlab = 'Total aligned reads (in thousand) to the spike-in genome', ylab = '', color = 'Replicate', plot.title = 'Reads aligned to spike-in genome'){
    require(ggplot2)
    
    as.data.frame(df) %>% 
      mutate(sample_replicate = ifelse(is.na(sample_replicate), 1, sample_replicate)) %>% 
      mutate(x = df[,var.x]/1e3, y = df[,var.y], color = factor(sample_replicate)) %>% 
      subset(target %in% tgt) %>%
      ggplot(mapping = aes(x = x, y = y, color = color))+
      geom_point()+
      labs(
        x = xlab, y = ylab, 
        color = color, title = plot.title)+
      theme_bw()+
      theme(
        text = element_text(size = 8),
        axis.text.x = element_text(angle = 90, hjust = 1),
        strip.text.y = element_text(angle = 0),
        legend.position = 'top'
      )
    
  }
  
  figs[[qc]] <- lapply(c(targets,'IgG'), funcs[[qc]], df = meta)
  names(figs[[qc]]) <- c(targets,'IgG')
  
  
  ## Duplication rate ####
  qc <- 'dup_rate'
  
  funcs[[qc]] <- function(df, tgt, var.x = 'dedup_percent_duplication', var.y = 'sample_group', xlab = 'Duplication Rate', ylab = '', color = 'Replicate', plot.title = 'Duplication rate'){
    require(ggplot2)
    
    as.data.frame(df) %>% 
      mutate(sample_replicate = ifelse(is.na(sample_replicate), 1, sample_replicate)) %>% 
      mutate(x = df[,var.x], y = df[,var.y], color = factor(sample_replicate)) %>% 
      subset(target %in% tgt) %>%
      ggplot(mapping = aes(x = x, y = y, color = color))+
      geom_jitter(height = 0.1, width = 0)+
      labs(
        x = xlab, y = ylab, 
        color = color, title = plot.title)+
      theme_bw()+
      theme(
        text = element_text(size = 8),
        axis.text.x = element_text(angle = 90, hjust = 1),
        strip.text.y = element_text(angle = 0),
        legend.position = 'top'
      )
  }
  
  figs[[qc]] <- lapply(c(targets,'IgG'), funcs[[qc]], df = meta)
  names(figs[[qc]]) <- c(targets,'IgG')
  
  
  
  ## Estimated library size ####
  qc <- 'est_lib_size' 
  funcs[[qc]] <- function(df, tgt, var.x = 'dedup_estimated_library_size', var.y = 'sample_group', xlab = 'Esitmated library size (in million)', ylab = '', color = 'Replicate', plot.title = 'Estimated library size'){
    require(ggplot2)
    
    as.data.frame(df) %>% 
      mutate(sample_replicate = ifelse(is.na(sample_replicate), 1, sample_replicate)) %>% 
      mutate(x = df[,var.x]/1e6, y = df[,var.y], color = factor(sample_replicate)) %>% 
      subset(target %in% tgt) %>%
      ggplot(mapping = aes(x = x, y = y, color = color))+
      geom_jitter(height = 0.1, width = 0)+
      labs(
        x = xlab, y = ylab,
        color = color, title = plot.title)+
      theme_bw()+
      theme(
        text = element_text(size = 8),
        strip.text.y = element_text(angle = 0),
        legend.position = 'top'
      )
  }
  
  figs[[qc]] <- lapply(c(targets,'IgG'), funcs[[qc]], df = meta)
  names(figs[[qc]]) <- c(targets,'IgG')
  
  
  
  
}

## get fragment lengths ####
if (dir.exists('fragment_lengths')){
  length.list <- list.files('fragment_lengths/', pattern = 'fragment_lengths', full.names = T, recursive = T)
}
dat$frag_lens <- bind_rows(lapply(
  length.list,
  function(fname){
    if(file.size(fname) > 0 ){
      read.delim(fname, header = F, col.names = c('length','count'), colClasses = 'numeric') %>%
        mutate(
          weight = count/sum(count),
          id = gsub('.fragment_length.txt','',basename(fname)))
    }
  })) %>% 
  left_join(
    dat$meta %>% dplyr::select(id,sample_group, target, sample_replicate),
    by = 'id'
  )
funcs[['frag_lens_dens']] <- function(df, tgt, var.x = 'length', var.y = 'count', var.group = 'id', facet.row = 'sample_group', xlab = 'Fragment length (bp)', ylab = 'Occurrences', color = 'Replicate', plot.title = 'Fragment length'){
  require(ggplot2)
  
  df <- as.data.frame(df)
  p <- df %>%
    mutate(x = df[,var.x], y = df[,var.y], group = df[,var.group], color = factor(sample_replicate), facet.row = df[,facet.row]) %>%
    subset(target %in% tgt) %>%
    ggplot(mapping = aes(x=x, y=y, group = group, color = color))+
    geom_line(size=0.5)+
    labs(x=xlab, y=ylab, color=color, title = plot.title)+
    theme_bw()+
    theme(
      text = element_text(size = 8),
      strip.text.y = element_text(angle = 0),
      legend.position = 'top'
    )
  
  if (length(unique(df$sample_group)) > 1){
    p+ 
      facet_grid(facet.row ~ .)
  }else{
    p
  }
}


figs[['frag_lens_dens']] <- lapply(
  c(targets, 'IgG'), funcs[['frag_lens_dens']], df = dat$frag_lens
); names(figs[['frag_lens_dens']]) <- c(targets, 'IgG')


## get metrics of the entire set of original peaks ####
peak.metrics <- list.files('original_peaks', full.names = T, pattern = 'original_peaks')
if (length(peak.metrics) > 0){
  dat$npeaks <- bind_rows(lapply(peak.metrics, read.csv)) %>% 
    left_join(
      dat$meta %>% dplyr::select(id,sample_group, target, sample_replicate),
      by = 'id'
    )
}

qc <- 'count_peaks' 
funcs[[qc]] <- function(df, tgt, var.x = 'npeak', var.y = 'sample_group', var.color = 'toIgG', add_facet = facet_grid(~caller, scales = 'free'), var.shape = 'filtered', xlab = 'Total Peaks (in thousand)', ylab = '', scale_color = scale_color_manual('', values = c("IgG_controlled"="indianred", "Target_only"="steelblue")), shape = '', plot.title = 'Peaks from each sample'){
  require(ggplot2)
  
  df <- as.data.frame(df) %>% 
    mutate(sample_replicate = ifelse(is.na(sample_replicate), 1, sample_replicate))
  df %>% 
    mutate(x = df[,var.x]/1e3, y = df[,var.y], color = df[,var.color], shape = df[,var.shape]) %>% 
    subset(target %in% tgt) %>%
    ggplot(mapping = aes(x = x, y = y, color = color, shape = shape))+
    geom_jitter(height = 0.1, width = 0)+
    add_facet+
    scale_color+
    scale_shape_manual(values = c(filtered=19, unfiltered=1))+
    labs(
      x = xlab, y = ylab, 
      shape = shape, title = plot.title)+
    theme_bw()+
    theme(
      text = element_text(size = 8),
      axis.text.x = element_text(angle = 90, hjust = 1),
      legend.position = 'top'
    )
}

df <- dat$npeaks 

figs[[qc]] <- lapply(targets, funcs[[qc]], df = df, var.y = 'id')
names(figs[[qc]]) <- targets


## get metrics of the final set of original peak widths  ####
peak.widths <- list.files('original_peak_widths', full.names = T, pattern = 'original_peak_widths')
if (length(peak.widths) > 0){
  dat$wpeaks <- bind_rows(lapply(peak.widths, read.csv)) %>% 
    left_join(
      dat$npeaks %>% dplyr::select(file,id, sample_group, target, sample_replicate,caller),
      by = 'file'
    )
}

qc <- 'peak_width'

funcs[[qc]] <- function(df, tgt, var.x = 'length', var.y = 'id', var.color = 'caller', xlab = 'Peak Width (bp)', ylab = '', scale_color = scale_color_manual('', values = c(SEACR='darkblue', MACS2narrow='red', MACS2broad='brown')), plot.title = 'Peak width'){
  require(ggplot2)
  
  ## determine data range
  length.max <- as.data.frame(df) %>% 
    subset(target %in% tgt) %>% 
    group_by(sample_group, caller) %>% 
    reframe(
      qt = quantile(length, 0.75)
    ) %>% 
    reframe(
      max = max(qt)
    ) %>% 
    as.matrix() %>% 
    as.vector()
  
  ## make plots
  as.data.frame(df) %>% 
    mutate(x = .[,var.x], y = .[,var.y], color = .[,var.color]) %>% 
    subset(target %in% tgt) %>%
    ggplot(mapping = aes(x = x, y = y, color = color, weight = weight))+
    geom_violin(bw=5, trim = T, draw_quantiles = 0.5, orientation = 'y')+
    scale_color+
    coord_cartesian(xlim = c(0, length.max))+
    labs(x=xlab, y=ylab, title = plot.title)+
    theme_bw()+
    theme(
      text = element_text(size = 8),
      axis.text.x = element_text(angle = 45, hjust = 1),
      legend.position = 'top'
    )
}

df <- dat$wpeaks 

figs[[qc]] <- lapply(targets, funcs[[qc]], df = df)
names(figs[[qc]]) <- targets


## get reads in peak ####
peak.reads <- list.files('reads_in_peak', full.names = T, pattern = 'reads_in_peak')
dat$frip <- bind_rows(lapply(
  peak.reads,
  function(fname){
    if(file.size(fname) > 0 ){
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

qc <- 'peak_frip'

funcs[[qc]] <- function(df, tgt, var.x = 'FRiP', var.y = 'id', var.color = 'caller', var.shape = 'filtered', add_facet = NULL, xlab = 'Fraction of reads in peak', ylab = '', scale_color = scale_color_manual('', values = c(SEACR='darkblue', MACS2narrow='red', MACS2broad='brown')), scale_shape = scale_shape_manual('', values = c(filtered = 1, unfiltered = 2)), plot.title = 'Fraction of reads in peak'){
  require(ggplot2)
  
  as.data.frame(df) %>% 
    mutate(x = .[,var.x], y = .[,var.y], color = .[,var.color], shape = .[,var.shape]) %>% 
    subset(target %in% tgt) %>%
    ggplot(mapping = aes(x = x, y = y, color = color, shape = shape))+
    geom_jitter(height = 0.1, width = 0)+
    add_facet+
    scale_color+
    scale_shape+
    labs(
      x = xlab, y = ylab, 
      title = plot.title
    )+
    theme_bw()+
    theme(
      text = element_text(size = 8),
      axis.text.x = element_text(angle = 90, hjust = 1),
      plot.title = element_text(size = 10),
      legend.position = 'top'
    )
  
  
}

df <- dat$frip

figs[[qc]] <- lapply(targets, funcs[[qc]], df = df)
names(figs[[qc]]) <- targets


## get metrics of replicated peaks ####
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

## get metrics of consensus peaks ####
conp.metrics <- list.files('consensus_peaks', full.names = T, pattern = 'consensus_peaks')
if (length(conp.metrics) > 0){
  dat$nconps <- bind_rows(lapply(conp.metrics, read.csv))
}

## get consenus peaks ####
conp.bed <- list.files('consensus_beds/', full.names = T, pattern = '.bed')

conps <- lapply(
  conp.bed, function(fname){
    read.delim(fname, header = F, col.names = c('chrom', 'start', 'end', 'conp.id','sample.groups', 'npeak'))
  }
); names(conps) <- basename(conp.bed)

conps.gr <- lapply(
  conps, function(x){
    if(nrow(x)>0){
      makeGRangesFromDataFrame(x, ignore.strand = T, seqnames.field = 'chrom', start.field = 'start', end.field = 'end', starts.in.df.are.0based = T)
    }
  }
)

## make plots
{
  ## compare consensus peaks across sample groups ####
  qc <- 'rep2conp'
  
  funcs[[qc]] <- function(df, plot.title){
    df <- df %>% 
      mutate(
        shared = ifelse(grepl(',', sample.groups), 'shared', 'unique')
      ) %>% 
      rowwise() %>% 
      reframe(
        conp.id = conp.id,
        shared = shared,
        sample.group = unlist(strsplit(sample.groups, split = ','))
      )
    
    df %>% 
      ggplot(aes(y = sample.group, fill = factor(shared, levels = c('unique', 'shared'))))+
      geom_bar()+
      labs(y = '', x = 'peak count', fill = '', title = plot.title)+
      theme_bw(base_size = 8)
    
  }
  
  figs[[qc]] <- lapply(
    paste(rep(targets, each = 3), c('macs2_narrow_peaks.bed', 'macs2_broad_peaks.bed', 'seacr_peaks.bed'), sep = '.'),
    function(i){
      if (i %in% names(conps)){
        funcs[[qc]](conps[[i]], i)
      }else{
        plot_spacer()
      }
    }
  )
  
  ## Compare MACS2 and SEACR by consensus peaks ####
  qc <- 'seacr2macs'
  
  funcs[[qc]] <- function(peak.list){
    require(GenomicRanges)
    
    tgt <- gsub('.macs2_.*|.seacr_.*', '', names(peak.list)[1])
    
    if(sum(sapply(peak.list, length)>0) < 2){return(NULL)} # only one peak set is available
    
    peak.list <- peak.list[sapply(peak.list, length)>0]

    conp <- GRanges()
    for (pk in peak.list){
      conp <- c(conp, pk)
    }
    conp <- reduce(conp)
    
    overlap2conp <- as.data.frame(sapply(
      peak.list, function(pk){countOverlaps(conp, pk)>0}
    ))
    colnames(overlap2conp) <- gsub('_peaks', '', stringr::str_extract(colnames(overlap2conp), "macs.*peaks|seacr.*peaks"))
    
    venn::venn(
      overlap2conp,
      ggplot = T,
      box = F
    )+
      labs(title = tgt)+
      theme(
        text = element_text(size = 10),
        plot.margin = unit(rep(0.5, 4), 'line'),
        plot.title = element_text(size = 12)
      )
    
  }
  
  figs[[qc]] <- lapply(
    split(conps.gr, gsub('.macs2_.*|.seacr_.*', '', names(conps.gr))),
    funcs[[qc]]
  )
  
  
}

## genomic distribution of consensus peaks ####
if (dir.exists('consensus_annotation')){
  ann.txt <- list.files('consensus_annotation', pattern = 'annotation.txt', full.names = T, recursive = T)
}

dat$conp_ann <- lapply(ann.txt, function(fname){
  read.delim(fname, header = T) %>%   
    subset(!is.na(Annotation)) %>% 
    mutate(
      genomic.location = gsub(' .*', '', Annotation)
    ) 
  
}); names(dat$conp_ann) <- gsub('.annotation.txt','',basename(ann.txt))


qc <- 'conp_ann'
funcs[[qc]] <- function(df, plot.title = NULL){
  require(dplyr)
  require(ggplot2)
  
  df %>% 
    ggplot(
      aes(y = genomic.location)
    )+
    geom_bar()+
    scale_x_continuous(expand = c(0,0))+
    labs(x = 'Consensus peak count', y = '', title = plot.title)+
    theme_bw()+
    theme(
      text = element_text(size = 8),
      plot.title = element_text(size = 8)
    )
  
}

figs[[qc]] <- mapply(funcs[[qc]], dat$conp_ann, names(dat$conp_ann), SIMPLIFY = F)
names(figs[[qc]]) <- names(dat$conp_ann)


## Save results ####
saveRDS(funcs, 'funcs.rds')
saveRDS(figs, 'figs.rds')
saveRDS(dat, 'data.rds')

