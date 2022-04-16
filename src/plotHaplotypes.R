#!/usr/bin/env R
library(data.table)
library(ggplot2)
library(foreach)

dat.all <- fread('data/3003_G1.blast')

# only use perfect blast for now
dat.all <- dat.all[V5 == 100]

dat.all[, sample := tstrsplit(V1, split="_")[3]]
dat.all[, sample := as.numeric(sample)]

o <- foreach(i=1:96, .combine='rbind') %do% {
    dat <- dat.all[sample == i]
    dat.perf <- dat[dat[, .I[V5==max(V5)], by=V3]$V1]
    uniques <- dat.perf[, .N, by=V3][N==1, V3]
    dat.uniques <- dat.perf[V3 %in% uniques]
    dat.uniques[, tenKbBin := cut(V11, breaks=seq(0,1.6e6, 10000))]
    dat.uniques[, tenKbBin := as.numeric(tenKbBin)]
    dat.wide <- dcast(dat.uniques[, .N, by=list(V4, tenKbBin)], tenKbBin~V4)
    dat.wide[is.na(tig00000001), tig00000001 := 0]
    dat.wide[is.na(tig00000087j113rc), tig00000087j113rc := 0]
    dat.wide[, pctContig001 := tig00000001 / (tig00000001 + tig00000087j113rc)]
    dat.out <- dat.wide[, c("tenKbBin","pctContig001")]
    dat.out[, sample := i]
}

ggplot(o, aes(x=tenKbBin, y=sample, fill=pctContig001)) + geom_tile()