#!/usr/bin/env Rscript

# Author: @fwu19

options(stringsAsFactors = F)
library(dplyr)

args <- as.vector(commandArgs(T))
in_csv <- args[1]
out_csv <- args[2]
use_control <- as.logical(args[3])

ss <- read.csv(in_csv) %>%
    mutate(
        control = ifelse(use_control, control, ""),
        is_control = ifelse(use_control, is_control, false)
    )

ss %>%
    dplyr::select(id, group, replicate, control, single_end, fastq_1, fastq_2, is_control) %>%
    write.table(out_csv, sep = ',', quote = F, row.names = F)
