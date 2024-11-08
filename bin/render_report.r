#!/usr/bin/env Rscript

options(stringsAsFactor=F)
library(rmarkdown)
args <- as.vector(commandArgs(T))

rmarkdown::render(args[1], params = list(dat.dir = args[2]))
