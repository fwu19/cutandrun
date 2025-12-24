#!/usr/bin/env Rscript

# Author: @fwu19

options(stringsAsFactors = F)
library(dplyr)

## functions ####
add_genomic_locations <- function(bed, ann_dir = 'peak_annotations'){
    peaks <- read.delim(bed, header = F)[1:5]
    colnames(peaks)[1:5] <- c('#chrom', 'start', 'end', 'peak_id', 'sample_groups')
    peaks <- peaks %>% 
        mutate(
            score = 1000,
            strand = '.'
        ) %>% 
        relocate(score, strand, .after = 'peak_id')
    
    txt <- file.path(ann_dir, gsub('.bed$', '.annotation.txt', basename(bed)))
    if (!file.exists(txt) | file.size(txt)==0){ return(peaks) }
    ann <- read.delim(txt)
    colnames(ann)[1] <- 'peak_id'
    
    peaks %>% 
        left_join(
            ann %>% 
                mutate(
                    genomic_location = gsub(' .*', '', Annotation)
                ) %>% 
                dplyr::select(peak_id, genomic_location),
            by = 'peak_id'
        )
        
}
    
add_differential_peaks <- function(peaks, dp){
    if (is.null(dp)){ return(peaks) }
    
    for (i in names(dp)){
        df <- dp[[i]]$df %>% 
            dplyr::select(conp.id, is.sig) 
        colnames(df) <- c('peak_id', i)
        peaks <- peaks %>% 
            left_join(df, by = 'peak_id')
        peaks[,i] <- plyr::mapvalues(
            peaks[,i],
            from = c(1, 0, -1, NA),
            to = c('Up', 'Not_DE', 'Down', 'Not_tested')
        )
    }
    
    
    return(peaks)
}    

## read arguments ####
args <- as.vector(commandArgs(T))
beds <- list.files('peaks', full.names = T)

## add peak annotations ####
if (dir.exists('peak_annotations')){
    peaks <- lapply(beds, add_genomic_locations, ann_dir = 'peak_annotations')
    names(peaks) <- basename(beds)
}else{
    peaks <- beds
}

## add differential peaks ####
if (dir.exists('differential_peaks')){
    rds <- list.files('differential_peaks/', full.names = T)
    if (length(rds) > 0){
        dp <- c()
        for (i in rds){
        dp <- c(dp, readRDS(i))
        }
        peaks <- mapply(add_differential_peaks, peaks, dp[names(peaks)], SIMPLIFY = F)
    }
}

## write out results ####
mapply(
    function(pk, out_file){
        pk %>% 
            write.table(out_file, sep = '\t', quote = F, row.names = F)
    }, peaks, names(peaks), SIMPLIFY = F
)
