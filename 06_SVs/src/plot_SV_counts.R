#!/usr/bin/env Rscript

library(ggplot2)
library(foreach)
library(doMC)
library(data.table)
library(ggthemes)

tsvfiles <- list.files(path='reference_alignment', pattern='*SVs.tsv', full.names=TRUE)

o <- foreach(f=tsvfiles, .combine='rbind') %do% {
    fread(f)
}

o[, accession := tstrsplit(ref_chr, split='\\|')[2]]
setkey(o, accession)

ids <- fread('chromosome_ids.txt', header=FALSE)
setnames(ids, c('accession','chr'))
setkey(ids, accession)

dat <- merge(o, ids)
dat[, accession := NULL]
dat[, ref_chr := NULL]
dat[, c('strain', 'query_chr') := tstrsplit(query_chr, split='_')]
setnames(dat, 'chr', 'ref_chr')

dat2 <- copy(dat)
dat2[, bin := cut(

# get CHR lengths



SVs <- dat[, .N, by=list(strain, query_chr, SV_type)]


chrLvls <- c(
'chrI','chrII','chrIII','chrIV','chrV',
'chrVI','chrVII','chrVIII','chrIX','chrX',
'chrXI','chrXII','chrXIII','chrXIV','chrXV',
'chrXVI','mitochondrion')

strainLvls <- c(
'BY','273614','CBS2888','CLIB219','CLIB413','I14','M22','PW5',
'RM','Y10','YJM145','YJM454','YJM978','YJM981','YPS1009',
'YPS163'
)


SVs[, query_chr := factor(query_chr, levels=chrLvls)]
SVs[, strain := factor(strain, levels=strainLvls)]

SVs[SV_type == 'deletion_mobile', SV_type := 'Deletion (mobile)']
SVs[SV_type == 'deletion_novel', SV_type := 'Deletion (novel)']
SVs[SV_type == 'insertion_mobile', SV_type := 'Insertion (mobile)']
SVs[SV_type == 'insertion_novel', SV_type := 'Insertion (novel)']
SVs[SV_type == 'duplication', SV_type := 'Duplication']
SVs[SV_type == 'inversion', SV_type := 'Inversion']
SVs[SV_type == 'contraction', SV_type := 'Contraction']
SVs[SV_type == 'transloc', SV_type := 'Translocation']

SVlvls <- c(
'Deletion (mobile)',
'Deletion (novel)',
'Insertion (mobile)',
'Insertion (novel)',
'Duplication',
'Inversion',
'Contraction',
'Translocation'
)
SVs[, SV_type := factor(SV_type, levels=SVlvls)]

fwrite(SVs, file='reference_alignment/SV_counts.tsv', quote=F, row.names=F, col.names=T, sep='\t')

g <- ggplot(data=SVs, aes(x=query_chr, y=N, fill=SV_type)) +
    geom_bar(stat='identity') +
    facet_grid(strain~.) +
    labs(x='Assembled Chromosome', y='Count', fill='Structural Variant') +
    theme_few(12) +
    theme(legend.position = 'bottom')

ggsave(g, file='reference_alignment/SV_counts.png', width=40, height=40, units='cm')
ggsave(g, file='reference_alignment/SV_counts.svg', width=40, height=40, units='cm')
