#!/usr/bin/env Rscript

# Author: @fwu19

options(stringsAsFactors = F)
library(dplyr)

## functions ####
add_metadata <- function(ss, meta_csv){
    ## update with metadata if provided
    # metadata contains a required columns id and optional columns: sample_group, sample_replicate, target, control, call_peak, call_rep_peak, call_con_peak
    
    if (grepl('dummy', meta_csv)){
        return(ss)
    }
    
    if (file_test('-f', meta_csv) & grepl('.csv$', meta_csv)){
        meta <- read.csv(meta_csv)
        ss <- ss %>% 
            inner_join(
                meta, by = 'id', suffix = c(".x", "")
            ) %>% 
            dplyr::select(!ends_with(".x")) 
    }
    
    return(ss)
    
}

## read arguments ####
args <- as.vector(commandArgs(T))
in_csv <- args[1]
out_csv <- args[2]
use_control <- as.logical(args[3]) # not used
meta_csv <- ifelse(length(args) > 3, args[4], '')


## generate sample sheet ####
ss <- read.csv(in_csv)

## check single end fastq
if (!'id' %in% colnames(ss)){
    stop ( 'Missing column id!' )
}else if ( ! 'fastq_1' %in% colnames(ss)){
    stop ( 'Missing column fastq_1!' )
}else if ( ! 'fastq_2' %in% colnames(ss)){
    ss$fastq_2 <- ''
}
ss$fastq_2 <- ifelse(is.na(ss$fastq_2), "", ss$fastq_2)
ss$single_end <- ifelse(ss$fastq_2 == "", 'true', 'false')

## add metadata if available
ss <- add_metadata(ss, meta_csv)

## identify controls
controls <- unique(ss$control)
ss <- ss %>% 
    mutate(
        group = paste(sample_group, target, sep = '_'), # to be compatible with nf-core
        replicate = sample_replicate, # to be compatible with nf-core
        control = ifelse(is.na(control), "", control), # convert NA to ""
        is_control = ifelse(id %in% controls, 'true', 'false'),
        control_group = ifelse(is_control, id, control),
        call_peak = ifelse(is_control == 'true', 'false', call_peak),
        call_rep_peak = ifelse(is_control == 'true', 'false', call_rep_peak),
        call_con_peak = ifelse(is_control == 'true', 'false', call_rep_peak)
    ) %>%
    dplyr::relocate(id, group, replicate, single_end, is_control, control_group, control, fastq_1, fastq_2, target, sample_group, sample_replicate, call_peak, call_rep_peak, call_con_peak)
 
## write sample sheet ####
ss %>% 
    write.table(out_csv, sep = ',', quote = F, row.names = F)

