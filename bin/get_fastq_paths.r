#!/usr/bin/env Rscript

# Author: @fwu19

options(stringsAsFactors = F)
library(dplyr)

## functions ####
get_fastqs <- function(
    fq_dirs, workflow, 
    r1_pattern = "_S[0-9]+(_L[0-9]+)?_R1_",
    r2_pattern = "_S[0-9]+(_L[0-9]+)?_R2_"
    ){
    ## get paths to fastq files
    fqs <- normalizePath(grep('undetermined', list.files(fq_dirs, recursive = T, full.names = T, pattern = "fastq.gz"), invert = T, value = T, ignore.case = T))
    fqs1 <- sort(grep(r1_pattern, fqs, value = T))
    fqs2 <- sort(grep(r2_pattern, fqs, value = T))
    
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
lst <- strsplit(args, split = '=')
for (x in lst){
    assign(x[1],x[2])
} # read arguments: workflow r1_pattern r2_pattern
rm(lst)

if(!exists('r1_pattern')){r1_pattern <- "_S[0-9]+(_L[0-9]+)?_R1_"}
if(!exists('r2_pattern')){r2_pattern <- "_S[0-9]+(_L[0-9]+)?_R2_"}

fqs <- list.files('fastq', full.names = T)
if(file_test('-f', fqs[1])){
    input_dirs <- scan(fqs[1], what = 'character')
}else if (file_test('-d', fqs[1])){
    input_dirs <- fqs[1]
}else{
    stop("--input_dir takes either path/to/fastq/dir or a file containing paths/to/fastq/dir (one path in each row)")
}

ss <- get_fastqs(input_dirs, workflow)

ss %>% 
    write.table('input.csv', sep = ',', quote = F, row.names = F)

