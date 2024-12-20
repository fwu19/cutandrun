#!/usr/bin/env Rscript

# Author: @fwu19

options(stringsAsFactors = F)
library(dplyr)

## functions ####
get_fastqs <- function(fq_dirs, workflow){
    ## get paths to fastq files
    fqs <- normalizePath(grep('undetermined', list.files(fq_dirs, recursive = T, full.names = T, pattern = "fastq.gz"), invert = T, value = T, ignore.case = T))
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
        mutate(
            sample_group = id,
            sample_replicate = 1,
            target = id,
            control = "",
            call_peak = 'true',
            call_rep_peak = 'true',
            call_con_peak = 'true'
        )
    
    if (workflow == 'process_controls'){
        ss <- ss %>% 
            filter(
                grepl("^PC_.*K562", id) # keep processing controls only
            )
        
    }else{
        ss <- ss %>% 
            filter(
                !grepl("^PC_.*K562", id) # exclude processing controls
            )
    }
    
    ## 
    return(ss)
}

## read arguments ####
args <- as.vector(commandArgs(T))
workflow <- args[1]
if(file_test('-f', args[2])){
    input_dirs <- scan(args[2], what = 'character')
}else if (file_test('-d', args[2])){
    input_dirs <- args[2]
}else{
    stop("--input_dir takes either path/to/fastq/dir or a file containing paths/to/fastq/dir (one path in each row)")
}

ss <- get_fastqs(input_dirs, workflow)

ss %>% 
    write.table('input.csv', sep = ',', quote = F, row.names = F)

