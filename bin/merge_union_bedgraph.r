#!/usr/bin/env Rscript

# Author: @fwu19

options(stringsAsFactors = F)
options(scipen = 99)
library(dplyr)

## read arguments ####
args <- as.vector(commandArgs(T))
in_bdg <- args[1]
out_bdg <- args[2]

df <- read.delim(in_bdg, header = F)
df$sum <- rowSums(df[,4:ncol(df)])
df[c('V1', 'V2', 'V3', 'sum')] %>% 
    write.table(out_bdg, sep = '\t', quote = F, row.names = F, col.names = F)