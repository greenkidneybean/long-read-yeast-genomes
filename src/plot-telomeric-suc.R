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

dat[, plate := tstrsplit(query, split="_")[2]]
dat[, well := tstrsplit(query, split="_")[3]]
dat[, well := as.numeric(well)]

dat.ag <- dat[, list(.N), by=list(query, db, target, plate, well, nReads)]

## UNFINISHED FROM HERE

dat.ag[, pct_mapped := (N/nReads), by=list(query, target)]
dat.wide <- dcast(dat.ag, query~target, value.var="N")




dat.wide <- dcast(dat.ag, plate+well~match, value.var="N")

dat.final <- melt(dat.wide, measure.vars=c("telo1", "telo2"), variable.name="telo-specific-region", value.name="reads")


ggplot(dat.final, aes(x=intact, y=reads, color=`telo-specific-region`)) + geom_point() +
labs(x="reads mapping to whole 