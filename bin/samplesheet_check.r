#!/usr/bin/env Rscript

# Author: @fwu19

options(stringsAsFactors = F)
library(dplyr)

## functions ####
get_fastq <- function(flowcells){
    flowcells <- unlist(strsplit(flowcells, split=' '))
    
    ## get paths to fastq files
    fqs <- grep('undetermined', list.files(flowcells, recursive = T, full.names = T, pattern = "fastq.gz"), invert = T, value = T, ignore.case = T)
    fqs1 <- sort(grep("_S[0-9]+(_L[0-9]+)?_R1_", fqs, value = T))
    fqs2 <- sort(grep("_S[0-9]+(_L[0-9]+)?_R2_", fqs, value = T))
    
    nfq1 <- length(fqs1)
    nfq2 <- length(fqs2)
    if (nfq1 == 0){
        stop ('No read 1 files found!')
    }
    
    if (nfq1 > 0 & nfq2 > 0 & nfq1 != nfq2){
        stop ('Found different numbers of read 1 and read 2 files!')
    }
    
    ## generate a sample sheet without metadata
    ss <- data.frame(
        fastq_1 = fqs1,
        single_end = ifelse(nfq2 == 0, 'true', 'false')
    ) %>% 
        mutate(
            id = gsub("_S[0-9]+(_L[0-9]+)?_R1_.*", "", basename(fastq_1)),
            fastq_2 = ifelse(single_end, "", fqs2)
        ) %>% 
        filter(
            !grepl("^PC_.*K562", id) # exclude processing controls
        ) %>% 
        mutate(
            sample_group = id,
            sample_replicate = 1,
            target = "",
            control = "",
            call_peak = 'true',
            call_rep_peak = 'true',
            call_con_peak = 'true'
        )
    
    ## 
    return(ss)
}

add_metadata <- function(ss, meta.csv=NULL){
    ## update with metadata if provided
    # metadata contains a required columns id and optional columns: sample_group, sample_replicate, target, control, call_peak, call_rep_peak, call_con_peak
    
    if (file_test('-f', meta.csv) & grepl('.csv$', meta.csv)){
        meta <- read.csv(meta.csv)
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
input <- args[1]
out_csv <- args[2]
use_control <- as.logical(args[3]) # not used
meta.csv <- ifelse(length(args) > 3, args[4], NULL)


## generate sample sheet ####
if (file_test('-f', input) & grepl('.csv$', input)){
    ss <- read.csv(input)
}else if (file_test('-d', input) ){
    ss <- get_fastq(input)    
}else {
    stop (paste(input, 'should be a .csv file or a path to fastq files!'))
}

ss <- add_metadata(ss, meta.csv)

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
    
ss %>% 
    write.table(out_csv, sep = ',', quote = F, row.names = F)

