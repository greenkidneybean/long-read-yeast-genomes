#!/usr/bin/env Rscript

library(ggplot2)
library(foreach)
library(doMC)
library(data.table)
library(ggthemes)

tsvfiles <- list.files(path='reference_alignment', pattern='*SVs.tsv', full.names=TRUE)

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
setkey(dat2, strain, query_chr)

chr_lengths <- fread('assembly_lengths.tsv', header=TRUE)
setnames(chr_lengths, c('strain','query_chr','query_chr_length'))
setkey(chr_lengths, strain, query_chr)
binSize <- 20000

dat2 <- merge(dat2, chr_lengths)
dat2[, dist_from_R_end := query_chr_length - query_stop]
dat2[, dist_from_L_end := query_start]
dat2[, idx := 1:.N]
dat2[, min_dist_to_end := min(dist_from_R_end, dist_from_L_end), by=idx]

o2 <- foreach(minSVsize=c(0,100,500,1000,5000,10000), .combine='rbind') %do% {
    foreach(binSize=c(10000,20000,50000), .combine='rbind') %do% {
        dat.tmp <- copy(dat2)
        dat.tmp[, bin := trunc(min_dist_to_end /binSize) ]
        dat.tmp.ag <- dat.tmp[size >= minSVsize, .N, by=list(SV_type,bin)]
        dat.tmp.ag[, 'binSize' := binSize]
        dat.tmp.ag[, 'minSVsize' := minSVsize][]
        return(dat.tmp.ag)
    }
}

o2[SV_type == 'deletion_mobile', SV_type := 'Deletion (mobile)']
o2[SV_type == 'deletion_novel', SV_type := 'Deletion (novel)']
o2[SV_type == 'insertion_mobile', SV_type := 'Insertion (mobile)']
o2[SV_type == 'insertion_novel', SV_type := 'Insertion (novel)']
o2[SV_type == 'duplication', SV_type := 'Duplication']
o2[SV_type == 'inversion', SV_type := 'Inversion']
o2[SV_type == 'contraction', SV_type := 'Contraction']
o2[SV_type == 'transloc', SV_type := 'Translocation']

o2[, SV_type := factor(SV_type, levels=SVlvls)]

g.bins <- ggplot(o2[binSize==20000 & minSVsize==500], aes(x=bin, y=N)) + geom_point() +
    facet_grid(SV_type~.) +
    scale_x_continuous(breaks=c(0,5,10,15,20,25,30,35,40), 
                       labels=c(0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8))+
    theme_few() +
    labs(x='Distance from nearest end of chromosome (Mb)',
         y='Structural Variants (\u2265 500 bp) within 20 kilobase bin')

ggsave(g.bins, file='reference_alignment/SV_bins.png', width=15, height=40, units='cm')
ggsave(g.bins, file='reference_alignment/SV_bins.pdf', width=15, height=40, units='cm', dpi=300)

SVs <- dat[size >= 500, .N, by=list(strain, query_chr, SV_type)]


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


SVs[, SV_type := factor(SV_type, levels=SVlvls)]

fwrite(SVs, 
    file='reference_alignment/chr-specific-SV-counts-500bp-min.tsv',
    quote=F, row.names=F, col.names=T, sep='\t')



g <- ggplot(data=SVs, aes(x=query_chr, y=N, fill=SV_type)) +
    geom_bar(stat='identity') +
    facet_grid(strain~.) +
    labs(x='Assembled Chromosome', y='Structural Variants (\u2265 500 bp)', fill='Structural Variant') +
    theme_few(12) +
    theme(legend.position = 'bottom')

ggsave(g, file='reference_alignment/SV_counts.png', width=40, height=40, units='cm')
ggsave(g, file='reference_alignment/SV_counts.pdf', width=40, height=40, units='cm', dpi=300)




## GENOME-WIDE SV COUNT FIGURES
genomeWide <- copy(SVs)
genomeWide <- genomeWide[, list('N'=sum(N)), by=list(strain, SV_type)]


# get all pairwise combos to fill in 0s which get dropped otherwise
allCombos <- CJ('strain'=strainLvls, 'SV_type'=SVlvls)
setkey(allCombos, strain, SV_type)
setkey(genomeWide, strain, SV_type)
allGenomeWideCounts <- merge(genomeWide, allCombos, all=TRUE)
allGenomeWideCounts[is.na(N), N := 0]


# save table
fwrite(allGenomeWideCounts, 
    file='reference_alignment/genome-wide-SV-counts-500bp-min.tsv',
    quote=F, row.names=F, col.names=T, sep='\t')


# re-factor
allGenomeWideCounts[, strain := factor(strain, levels=strainLvls)]
allGenomeWideCounts[, SV_type := factor(SV_type, levels=SVlvls)]


# group by SV type
g.genomewide1 <- ggplot(data=allGenomeWideCounts) +
    geom_bar(stat='identity', width = 0.8, position = position_dodge(width = 1), aes(x='1', y=N, fill=strain)) +
    labs(y='Structural Variants (\u2265 500 bp), Genome-Wide', fill='Isolate') +
    theme_few(12) +
    facet_grid(.~SV_type, switch='x') +
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
    theme(axis.title.x=element_blank()) +
    scale_x_discrete(breaks=NULL, labels=NULL)

ggsave(g.genomewide1, file='reference_alignment/SV_counts_genomewide1.png', width=40, height=15, units='cm')
ggsave(g.genomewide1, file='reference_alignment/SV_counts_genomewide1.pdf', width=40, height=15, units='cm', dpi=300)


# group by strain
g.genomewide2 <- ggplot(data=allGenomeWideCounts) +
    geom_bar(stat='identity', width = 0.8, position = position_dodge(width = 1), aes(x='1', y=N, fill=SV_type)) +
    labs(y='Structural Variants (\u2265 500 bp), Genome-Wide', fill='Structural Variant') +
    theme_few(12) +
    facet_grid(.~strain, switch='x') +
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
    theme(axis.title.x=element_blank()) +
    scale_x_discrete(breaks=NULL, labels=NULL)

ggsave(g.genomewide2, file='reference_alignment/SV_counts_genomewide2.png', width=40, height=15, units='cm')
ggsave(g.genomewide2, file='reference_alignment/SV_counts_genomewide2.pdf', width=40, height=15, units='cm', dpi=300)