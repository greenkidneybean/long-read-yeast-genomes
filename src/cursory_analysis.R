#!/usr/bin/env Rscript

library(data.table)
library(ggplot2)
library(foreach)
library(doMC)
library(RColorBrewer)
library(zoo)
library(ggthemes)


filenames <- list.files('data/imputed', pattern="*.txt", full.names=TRUE)

o <- foreach(f = filenames, .combine='rbind') %do% {
    dat <- fread(f)[chrID==4]
    dat[]
}