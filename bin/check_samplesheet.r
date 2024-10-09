#!/usr/bin/env Rscript

# Author: @fwu19

options(stringsAsFactors = F)

args <- as.vector(commandArgs(T))
in_csv <- args[1]
out_csv <- args[2]
use_control <- as.logical(args[3])

ss <- read.csv(in_csv)
if (!use_control){
    ss$control <- ""
    ss$is_control <- 0
}

write.table(ss, out_csv, sep = ',', quote = F, row.names = F)