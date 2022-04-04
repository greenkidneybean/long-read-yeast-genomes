#!/usr/bin/env bash
library(data.table)
library(ggplot2)
library(foreach)

readThreshold <- 50000

filenames <- list.files(path='data/', pattern='*.blast', full.names=TRUE)

blast <- foreach(file=filenames, .combine='rbind') %do% {
    fread(file)
}


blast[, c('V7', 'V8', 'V9', 'V10') := NULL]
setnames(blast, c('query','db', 'number', 'target', 'PID', 'length', 'start', 'stop', 'evalue', 'bits'))

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

dat.ag <- dat[nReads >= readThreshold & PID >= 99 & length >= 120, list(.N), by=list(query, target, plate, well, nReads)]
# find strains with > 1 read mapping to (both unique sequences)
allCandidates <- unique(dat.ag[, query])
twoSucs <- unique(dat.ag[target=='unique1'][N > 1][,query])
oneSuc <- unique(dat.ag[! query %in% twoSucs][, query])
allStrains <- unique(dat[,query])

dat.ag[query %in% twoSucs, classification := 'twoSucs']
dat.ag[query %in% oneSuc, classification := 'oneSuc']
dat.ag[! query %in% oneSuc][! query %in% twoSucs, classification := 'indeterminate']


dat[query %in% twoSucs, classification := 'twoSucs']
dat[query %in% oneSuc, classification := 'oneSuc']
dat[is.na(classification), classification := 'indeterminate']


ggplot(dat[PID >= 99 & length >= 120, list(.N), by=list(query, target, plate, well, nReads, classification)], aes(x=nReads, y=N, color=classification)) + geom_point(shape=21, alpha=0.50) +
labs(x="total reads genome-wide", y="reads aligned to target via BLAST", title="nReads >= 50k, PID >= 99, length >= 120") +
facet_wrap(~target, ncol=2, scales='free')

classifications <- dat[!duplicated(dat[, c('query', 'plate', 'well', 'classification')])]

fwrite(classifications, file='suc-classificiation.txt', quote=F, row.names=F, col.names=T, sep='\t')

