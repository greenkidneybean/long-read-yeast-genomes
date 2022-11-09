#!/usr/bin/env Rscript

library(data.table)
library(ggplot2)
library(foreach)
args <- commandArgs(trailingOnly=TRUE)
filename <- args[1]

dat <- fread(filename)
setnames(dat, c('sample','mapped','unmapped'))
dat[, N := sum(mapped, unmapped)]
dat[, RPM := 1e6 * (mapped / N)]
dat[, sample2 := sample]

#ggplot(dat, aes(x=1, y=RPM)) + geom_violin()

fwrite(dat[,c('sample','sample2','RPM')], file="phenotypes.txt", quote=F, row.names=F, col.names=F, sep="\t")