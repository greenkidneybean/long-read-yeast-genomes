#!/usr/bin/env Rscript

library(data.table)
library(ggplot2)
library(foreach)

#window_size <- 2000 



files <- list.files(pattern="*.call")


firstfile <- files[1]


# get first file
keep <- fread(firstfile)
setnames(keep, "#CHROM", "CHROM")
setkey(keep, CHROM, POS, REF, ALT)

sample <- paste(strsplit(firstfile, "_|\\.")[[1]][2:3], collapse='_')







# iterate over other files
for(filename in files) {
    dat <- fread(filename)
    setnames(dat, c('CHROM','POS','REF','ALT','GT'))
    sample <- paste(strsplit(filename, "_|\\.")[[1]][2:3], collapse='_')
    setnames(dat, 'GT', sample)
    setkey(dat, CHROM, POS, REF, ALT)
    keep <- merge(keep, dat, all=TRUE)
}

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
