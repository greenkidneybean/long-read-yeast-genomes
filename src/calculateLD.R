#!/usr/bin/env Rscript

library(data.table)
library(ggplot2)
library(foreach)
setwd("data/vcf/chrIV")
window_size <- 2000 


files <- list.files(pattern='.vcf')
firstfile <- files[1]

assignGTs <- function(DT, windowSize){
    DT[, bin := cut(POS, breaks=seq(0,windowSize+max(DT$POS), windowSize))]
    DT <- DT[, list('nRef'=sum(GT==0), 'nAlt'=sum(GT==1)), by=bin]
    DT[, tot := nRef + nAlt]
    DT[, fracRef := nRef/tot]
    DT[, fracAlt := nAlt/tot]
    DT[fracRef > 0.9, GT := 0]
    DT[fracAlt > 0.9, GT := 1]
    DT <- DT[, c('bin','GT')]
    return(DT[])
}

# get first file
keep <- fread(firstfile)
setnames(keep, c('CHROM','POS',' ID','REF','ALT','QUAL','FILTER','INFO','FORMAT','GT'))
keep <- assignGTs(keep, window_size)
keep <- keep[, c('bin')]
setkey(keep, bin)
#ggplot(keep, aes(x=bin, y=GT)) + geom_point()






# iterate over other files
for(filename in files) {
    dat <- fread(filename)
    setnames(dat, c('CHROM','POS',' ID','REF','ALT','QUAL','FILTER','INFO','FORMAT','GT'))
    dat <- assignGTs(dat, window_size)
    dat <- dat[, c('bin','GT')]
    sample <- paste0(strsplit(filename, split="_")[[1]][1:3], collapse='_')
    setnames(dat, 'GT', sample)
    setkey(dat, bin)
    keep <- merge(keep, dat, all=TRUE)
}

sampleNames <- colnames(keep)[2:961]

# Check mapping bias
keep[, nRef := apply(.SD, 1, function(x) sum(x == 0, na.rm=TRUE)), .SDcols=sampleNames]
keep[, nAlt := apply(.SD, 1, function(x) sum(x == 1, na.rm=TRUE)), .SDcols=sampleNames]
keep[, p := nRef / (nRef + nAlt)]
keep[, q := nAlt / (nRef + nAlt)]
# keep.long <- melt(keep[nRef > 250 & nAlt > 250 & nRef < 500 & nAlt < 500][,c('bin','nRef','nAlt')], measure.vars=c('nRef','nAlt'))
keep.long <- melt(keep[,c('bin','nRef','nAlt')], measure.vars=c('nRef','nAlt'))
ggplot(keep.long, aes(x=bin, y=value, color=variable)) + geom_point()

# disable scientific notation
options(scipen=999)

bins <- as.numeric(sapply(strsplit(tstrsplit(keep[,bin], split=",")[[1]], split="\\("), "[[", 2))

binNames <- paste0('bin_', as.character(bins))

keep.t <- as.data.table(t(as.matrix(keep[,.SD, .SDcols=sampleNames])))
setnames(keep.t, binNames)

# get bin 
keep[, bin := binNames]
setkey(keep, bin)
keep[.('bin_12000'), nRef]
keep[.('bin_12000'), nAlt]
rownames(keep.t) <- sampleNames

fwrite(keep.t, file="genotypes.mat", quote=F, row.names=T, col.names=T, sep="\t", na="NA")

getLD <- function(DT, bin1, bin2, alleleCountDT) {
    p12 <- alleleCountDT[.(bin1), p] * alleleCountDT[.(bin2), q]
    p21 <- alleleCountDT[.(bin2), p] * alleleCountDT[.(bin1), q]
    pA_x_pB <- p12 * p21
    DT.sub <- data.table('bin1'=keep.t[,get(bin1)], 'bin2'=keep.t[,get(bin2)])
    DT.sub <- DT.sub[! is.na(bin1) & ! is.na(bin2)]

    pAB <- nrow(DT.sub[bin1==0 & bin2==0]) / nrow(DT.sub)
    D <- pAB - pA_x_pB
    data.table(bin1, bin2, p12, p21, pA_x_pB, pAB, D)[]
}

o <- foreach(firstBin=binNames[1:10], .combine='rbind') %do% {
    foreach(secondBin=binNames[1:10], .combine='rbind') %do% {
        getLD(keep.t, firstBin, secondBin, keep)
    }
}

getLD(keep.t, 'bin_1486000', 'bin_1488000', keep)

keep.t[]
recombination_rates <- new.env()
for(chromosome in chromosomes) {
    recombination_rates[[chromosome]] <- sum(bed[chr==chromosome]$M)       # convert c (cM per Megabase) to Morgans
}


fwrite(keep.t, file="genotypes.mat", quote=F, row.names=T, col.names=T, sep="\t" na="NA")

foreach(i=1:960, .combine=

ggplot(data=keep, aes(x=
#bin    #samplename,GT

samples <- colnames(keep)
samples <- samples[2:length(samples)]

o <- foreach(sample=samples, .combine='rbind') %do% {
    genotypes <- keep[[sample]]
    nRef <- sum(genotypes==0, na.rm=T)
    nAlt <- sum(genotypes==1, na.rm=T)
    total <- length(genotypes)
    nNA = total - nRef - nAlt
    data.table('sample'=sample, 'nRef'=nRef, 'nAlt'=nAlt, 'nNA'=nNA)
}

refBias <- sum(o$nRef) / (sum(o$nRef) + sum(o$nAlt))