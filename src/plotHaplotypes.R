#!/usr/bin/env R
library(data.table)
library(ggplot2)
library(foreach)

dat.all <- fread('data/3003_G1.blast')

# only use perfect blast for now
dat.all[, sample := tstrsplit(V1, split="_")[3]]
dat.all[, sample := as.numeric(sample)]
dat.all[, c('V1','V2') := NULL]
setnames(dat.all, 'V4', 'contig')
setnames(dat.all, 'V3', 'readID')
setnames(dat.all, 'V5', 'PID')
setnames(dat.all, 'V6', 'matchLength')
dat.all[, c('V7','V8','V9','V10') := NULL]
setnames(dat.all, 'V11', 'contigStart')
setnames(dat.all, 'V12', 'contigEnd')
setnames(dat.all, 'V13', 'Evalue')
setnames(dat.all, 'V14', 'bits')

minMatchLength <- 150
minPID <- 98.0
windowSize <- 30000

o <- foreach(i=1:96, .combine='rbind') %do% {
    dat <- dat.all[sample == i][matchLength >= minMatchLength][PID >= minPID]
    dat.perf <- dat[dat[, .I[bits==max(bits)], by=readID]$V1]
    dat.perf <- dat.perf[readID %in% dat.perf[, .N, by=readID][N==1, readID]]
    #uniques <- dat.perf[, .N, by=readID][N==1, readID]
    #dat.perf <- dat.perf[readID %in% uniques]
    dat.perf[, genomeWindow := cut(contigStart, breaks=seq(0,1.6e6, windowSize))]
    dat.perf[, genomeWindow := as.numeric(genomeWindow)]
    dat.wide <- dcast(dat.perf[, .N, by=list(contig, genomeWindow)], genomeWindow~contig, fill=0)
    setnames(dat.wide, c('genomeWindow', 'contig_A', 'contig_B'))
    dat.wide[is.na(contig_A), contig_A := 0]
    dat.wide[is.na(contig_B), contig_B := 0]
    dat.wide[, total := contig_A + contig_B]
    dat.wide[, fractionA := contig_A / (contig_A + contig_B)]
    dat.wide[, sample := i]
}

g <- ggplot(o[total > 5], aes(x=genomeWindow, y=sample, fill=fractionA)) + geom_tile()

ggsave(g, file="haplotypes-test.png", width=24, height=18, units='cm')