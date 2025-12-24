#!/usr/bin/env Rscript

options(stringsAsFactors = F)
options(scipen = 99)
library(dplyr)
library(ggplot2)
library(patchwork)

## functions ####
generate_count_matrix <- function(conp.bed, infiles, ids = NULL){

    conp <- read.delim(conp.bed, header = F, col.names = c('chrom', 'start', 'end', 'conp.id', 'sample.groups', 'npeak')) %>%
        mutate(length = end - start)

    cts <- do.call(cbind, lapply(
        infiles, function(fname){
            df <- read.delim(fname, header = T, comment.char = '#')
            stopifnot(identical(conp$conp.id, df$Geneid))
            df[7]
        }
    ))
    if (!is.null(ids)){
        colnames(cts) <- ids
    }

    return(cbind(conp, cts))
}


wrapper_one_conp <- function(ss, tgt, conp.bed, count.txts){
    cts <- generate_count_matrix(conp.bed, count.txts, ids = gsub('.*macs2_narrow_peaks.|.*macs2_broad_peaks.|.*seacr_peaks.|.fragmentCounts.txt', '', basename(count.txts)))
    
    cts %>%
        write.table(
            gsub('\\.bed$', '.raw_counts.txt', basename(conp.bed)),
            sep = '\t', quote = F, row.names = F
        )

}


## read arguments ####
args <- as.vector(commandArgs(T)) # ss, cmp, path/to/fragmentCounts.txt, path/to/conp.bed

lst <- strsplit(args, split = '=')
for (x in lst){
    assign(x[1],x[2])
} # read arguments: ss, comparison, rds, fdr, fc, fdr2, fc2
rm(lst)

ss <- read.csv(ss_csv)

conp.bed.all <- list.files('conp/', full.names = T)
count.txts.all <- list.files('counts/', full.names = T)

for (conp.bed in conp.bed.all){
    count.txts <- grep(gsub('.bed$', '', basename(conp.bed)), count.txts.all, value = T)
    wrapper_one_conp(ss, tgt, conp.bed, count.txts)

}


