#!/usr/bin/env Rscript

library(ggplot2)
library(foreach)
library(doMC)
library(data.table)
library(ggthemes)

tsvfiles <- list.files(path='MUMandCo_output', pattern='*.tsv', full.names=TRUE)

o <- foreach(f=tsvfiles, .combine='rbind') %do% {
    fread(cmd=paste0('cut -f 1-8 ', f))

}
o[, accession := tstrsplit(ref_chr, split='\\|')[2]]
setkey(o, accession)

ids <- fread('chromosome_ids.txt', header=FALSE)
setnames(ids, c('accession','chr'))
setkey(ids, accession)

dat <- merge(o, ids)
dat[, accession := NULL]
dat[, ref_chr := NULL]
dat[, strain := tstrsplit(query_chr, split='_')[1]]
setnames(dat, 'chr', 'ref_chr')

## OK

setcolorder(dat, c('ref_chr','ref_start','ref_stop','size',
                   'query_chr','query_start','query_stop',
                   'SV_type'))

SVs <- dat[, .N, by=list(strain, ref_chr, SV_type)]


chrLvls <- c(
'chrI','chrII','chrIII','chrIV','chrV',
'chrVI','chrVII','chrVIII','chrIX','chrX',
'chrXI','chrXII','chrXIII','chrXIV','chrXV',
'chrXVI','mitochondrion')

SVs[, ref_chr := factor(ref_chr, levels=chrLvls)]

ggplot(data=SVs, aes(x=ref_chr, y=N, fill=SV_type)) +
    geom_bar(stat='identity') +
    facet_grid(strain~.) +
    labs(x='Chromosome', y='Count', fill='Structural Variant') +
    theme_few(12)