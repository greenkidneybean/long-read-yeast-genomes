#!/usr/bin/env bash
library(data.table)
library(ggplot2)
library(foreach)

filenames <- list.files(path='data/', pattern='*.blast', full.names=TRUE)

blast <- foreach(file=filenames, .combine='rbind') %do% {
    fread(file)
}

blast[, c('V6', 'V7', 'V8', 'V9', 'V10', 'V11', 'V12') := NULL]
setnames(blast, c('query','db', 'number', 'target', 'PID', 'evalue', 'bits'))

filenames <- list.files(path='data/', pattern='*summary.txt', full.names=TRUE)
nreads <- foreach(file=filenames, .combine='rbind') %do% {
    fread(file)
}
setnames(nreads, c('query', 'nReads'))

setkey(blast, query)
setkey(nreads,query)

dat <- merge(blast, nreads)
dat[, db := NULL]

dat[, plate := tstrsplit(query, split="_")[2]]
dat[, well := tstrsplit(query, split="_")[3]]
dat[, well := as.numeric(well)]

dat.ag <- dat[, list(.N), by=list(query, target, plate, well, nReads)]

ggplot(dat.ag, aes(x=nReads, y=N, color=target)) + geom_point(shape=21, alpha=0.85) +
labs(x="total reads genome-wide", y="reads aligned to target via BLAST")

