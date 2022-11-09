#!/usr/bin/env Rscript

library(data.table)
library(ggplot2)
library(foreach)

window_size <- 1000 



files <- list.files(pattern="*.call")

assignGTs <- function(DT, windowSize){
    sampleName <- colnames(DT)[5]
    setnames(DT, c("CHROM", "POS", "REF", "ALT", "GT"))
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

firstfile <- files[1]


# get first file
keep <- fread(firstfile)
keep <- assignGTs(keep, window_size)
keep <- keep[, c('bin')]
setkey(keep, bin)






# iterate over other files
for(filename in files) {
    dat <- fread(filename)
    setnames(dat, c('CHROM','POS','REF','ALT','GT'))
    dat <- assignGTs(dat, window_size)
    dat <- dat[, c('bin','GT')]
    sample <- paste(strsplit(filename, "_|\\.")[[1]][1:3], collapse='_')
    setnames(dat, 'GT', sample)
    setkey(dat, bin)
    keep <- merge(keep, dat, all=TRUE)
}

keep <- keep[!duplicated(bin)]

sampleNames <- colnames(keep)[2:length(colnames(keep))]

# disable scientific notation
options(scipen=999)

bins <- as.numeric(sapply(strsplit(tstrsplit(keep[,bin], split=",")[[1]], split="\\("), "[[", 2))

binNames <- paste0('bin_', as.character(bins))

keep.t <- as.data.table(t(as.matrix(keep[,.SD, .SDcols=sampleNames])))
setnames(keep.t, binNames)

# get bin 
# keep[, bin := binNames]
# setkey(keep, bin)
# keep[.('bin_12000'), nRef]
# keep[.('bin_12000'), nAlt]
rownames(keep.t) <- sampleNames

fwrite(keep.t, file="genotypes.mat", quote=F, row.names=T, col.names=T, sep="\t", na="NA")

## END HERE
