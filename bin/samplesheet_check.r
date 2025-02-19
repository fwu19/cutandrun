#!/usr/bin/env Rscript

# Author: @fwu19

options(stringsAsFactors = F)
library(dplyr)

## functions ####
add_col <- function(df, col, default.value){
    if(!col %in% colnames(df)){
        df[,col] <- default.value
    }
    return(df)
}

add_metadata <- function(ss, meta_csv){
    ## update with metadata if provided
    # metadata contains a required columns id and optional columns: sample_group, sample_replicate, target, control, call_peak, call_rep_peak, call_con_peak
    
    if (file_test('-f', meta_csv) & grepl('.csv$', meta_csv) & !grepl('dummy', meta_csv)){
        meta <- read.csv(meta_csv)
        ss <- ss %>% 
            inner_join(
                meta, by = 'id', suffix = c(".x", "")
            ) %>% 
            dplyr::select(!ends_with(".x")) 
    }
    
    ## add missing columns
    ss <- ss %>% 
        add_col('sample_group', ss$id) %>% 
        add_col('target', ss$id) %>% 
        add_col('sample_replicate', 1) %>% 
        add_col('control', "") %>% 
        add_col('call_peak', 'true') %>% 
        add_col('call_rep_peak', 'true') %>% 
        add_col('call_con_peak', 'true') %>% 
        mutate(
            group = paste(sample_group, target, sep = '_'), # to be compatible with nf-core
            replicate = sample_replicate
        )
    
    if (identical(ss$target, ss$id)){
        ss$group <- ss$id
    }else{
        ss$group <- paste(ss$sample_group, ss$target, sep = '_')
    }
    ss$replicate <- ss$sample_replicate
    
    ## identify controls and modify accordingly
    controls <- unique(ss$control)
    ss <- ss %>% 
        mutate(
            control = ifelse(is.na(control), "", control), # convert NA to ""
            is_control = ifelse(id %in% controls, 'true', 'false'),
            control_group = ifelse(is_control == 'true', id, control),
            call_peak = ifelse(is_control == 'true', 'false', call_peak),
            call_rep_peak = ifelse(is_control == 'true', 'false', call_rep_peak),
            call_con_peak = ifelse(is_control == 'true', 'false', call_rep_peak)
        ) %>%
        dplyr::relocate(id, group, replicate, single_end, is_control, control_group, control, fastq_1, fastq_2, target, sample_group, sample_replicate, call_peak, call_rep_peak, call_con_peak)
    
    
    return(ss)
    
}

add_metadata_process_controls <- function(ss){
    ss$flowcell <- basename(gsub('.Unaligned.*', '', ss$fastq_1))
    ss$sample_group <- ifelse(as.integer(gsub('_.*', '', ss$flowcell)) > 211110, 'PE50', 'PE25')
    ss$sample_label <- gsub('^PC_|K562_','',basename(ss$id))
    ss$target <- case_when(
        grepl('H3K27|H3K23|K27', ss$id, ignore.case = T) ~ 'H3K27me3',
        grepl('Pol', ss$id, ignore.case = T) ~ 'Pol2Ser5',
        grepl('Myc', ss$id, ignore.case = T) ~ 'cMyc',
        grepl('CTCF', ss$id, ignore.case = T) ~ 'CTCF',
        grepl('IgG', ss$id, ignore.case = T) ~ 'IgG',
        TRUE ~ 'TBD'
    )
    ss <- ss %>% 
        group_by(flowcell,target) %>% 
        mutate(
            i = 1:n(),
            n = n()
        ) %>% 
        mutate(
            sample_replicate = paste0(
                sapply(
                    strsplit(flowcell, split = '_'), 
                    function(v){
                        paste(c(stringr::str_sub(v[1], -6), v[length(v)]), collapse = '')
                    }),
                ifelse(n == 1, "", letters[i]))
        ) %>% 
        mutate(
            id = paste(target, sample_replicate, sep = '_')
        )
    
    ## add missing columns
    ss <- ss %>% 
        add_col('control', "") %>% 
        add_col('call_peak', 'true') %>% 
        add_col('call_rep_peak', 'true') %>% 
        add_col('call_con_peak', 'true') %>% 
        mutate(
            group = paste(sample_group, target, sep = '_'), # to be compatible with nf-core
            replicate = sample_replicate
        )
    
    ## identify controls and modify accordingly
    controls <- unique(ss$control)
    ss <- ss %>% 
        mutate(
            control = ifelse(is.na(control), "", control), # convert NA to ""
            is_control = ifelse(id %in% controls, 'true', 'false'),
            control_group = ifelse(is_control == 'true', id, control),
            call_peak = ifelse(is_control == 'true', 'false', call_peak),
            call_rep_peak = ifelse(is_control == 'true', 'false', call_rep_peak),
            call_con_peak = ifelse(is_control == 'true', 'false', call_rep_peak)
        ) %>%
        dplyr::select(-c(i,n)) %>% 
        dplyr::relocate(id, group, replicate, single_end, is_control, control_group, control, fastq_1, fastq_2, target, sample_group, sample_replicate, call_peak, call_rep_peak, call_con_peak, sample_label)
    
    
    return(ss)
    
}

## read arguments ####
args <- as.vector(commandArgs(T))
in_csv <- args[1]
out_csv <- args[2]
workflow <- args[3]
meta_csv <- ifelse(length(args) > 3, args[4], '')


## generate sample sheet ####
ss <- read.csv(in_csv, colClasses = 'character')

## check single end fastq
if (!'id' %in% colnames(ss)){
    stop ( 'Missing column id!' )
}else if ( ! 'fastq_1' %in% colnames(ss)){
    stop ( 'Missing column fastq_1!' )
}else if ( ! 'fastq_2' %in% colnames(ss)){
    ss$fastq_2 <- ''
}
ss <- ss %>% 
    add_col('fastq_2', "") %>% 
    mutate(single_end = ifelse(fastq_2 == "", 'true', 'false'))


## add metadata if available
if (workflow == "cutandrun"){
    ss <- add_metadata(ss, meta_csv)
}else if (workflow == "process_controls"){
    ss <- add_metadata_process_controls(ss)
}


## write sample sheet ####
ss %>% 
    write.table(out_csv, sep = ',', quote = F, row.names = F)

