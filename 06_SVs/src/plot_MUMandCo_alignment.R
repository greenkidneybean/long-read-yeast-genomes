#!/usr/bin/env Rscript

library(ggplot2)
library(foreach)
library(doMC)
library(data.table)
library(ggthemes)

minLen <- 10000

coordsfiles <- list.files(path='reference_alignment/', pattern='*.coords', full.names=TRUE)

o <- foreach(f=coordsfiles, .combine='rbind') %do% {
    fread(f, skip='[S1]')
    # warnings from fread because final two columns (V14 and V15) are unnamed, this is fine
    # and new columns are added automatically
}

# fix column names
setnames(o, c('start_ref','end_ref','start_query','end_query','len_ref','len_query','pct_id',
                'ref_chr_len','query_chr_len','cov_r','cov_q','ref_strand','query_strand',
                'ref_contig_tmp','query_chr'))

# remove unneeded columns
o[, 'cov_r' := NULL]
o[, 'cov_q' := NULL]
o[, 'ref_strand' := NULL]
o[, 'query_strand' := NULL]

o[, accession := tstrsplit(ref_contig_tmp, split='\\|')[2]]
o[, 'ref_contig_tmp' := NULL]
setkey(o, accession)

ids <- fread('chromosome_ids.txt', header=FALSE)
setnames(ids, c('accession','chr'))
setkey(ids, accession)

dat <- merge(o, ids)
setnames(dat, 'chr', 'ref_chr')
dat[, accession := NULL]
dat[, strain := tstrsplit(query_chr, split='_')[1]]


chrLvls <- c(
'chrI','chrII','chrIII','chrIV','chrV',
'chrVI','chrVII','chrVIII','chrIX','chrX',
'chrXI','chrXII','chrXIII','chrXIV','chrXV',
'chrXVI','mitochondrion')




## OK

# Note A=Query, B=Ref


dat[, query_chr := tstrsplit(query_chr, split='_')[2]]
dat[, chrA := factor(query_chr, levels=chrLvls)]
dat[, chrB := factor(ref_chr, levels=chrLvls)]



dat[, sequentialID := as.numeric(chrB)]
chrLengths <- unique(dat[, c('ref_chr_len', 'sequentialID')])
chrLengths[, cumulativeLength := cumsum(ref_chr_len)]
chrLengths[, cumulativeLength := shift(cumulativeLength, n=1L, type='lag')]
chrLengths[is.na(cumulativeLength), cumulativeLength := 0]

setkey(chrLengths, sequentialID)
setkey(dat, sequentialID)

dat2 <- merge(chrLengths, dat)


## ?ok?

dat2[,startB := start_ref + cumulativeLength]
dat2[,stopB := end_ref + cumulativeLength]



boxes <- data.table(x1=c(0,unique(dat2[cumulativeLength != 0, cumulativeLength])))
boxes[, x2 := shift(x1, n=1L, type='lead', fill=max(dat2$stopB))]
boxes[, midpoint := (x1+x2)/2]
boxes[, chr := chrLvls]
boxes[chr=='mitochondrion', chr := 'm']

dat2[, spanA := end_query - start_query]
dat2[, spanB := end_ref - start_ref]

strainLvls <- c(
'BY','273614','CBS2888','CLIB219','CLIB413','I14','M22','PW5',
'RM','Y10','YJM145','YJM454','YJM978','YJM981','YPS1009',
'YPS163'
)

dat2[, strain := factor(strain, levels=strainLvls)]

plot_all <- function() {
    g <- ggplot(data=dat2) + 
        geom_segment(aes(x=startB, xend=stopB, y=start_query, yend=end_query, color=chrA)) +
        theme_few() +
        scale_x_continuous(breaks=boxes$midpoint, labels=boxes$chr, expand = c(0, 0), limits = c(0, NA)) +
        scale_y_continuous(breaks=NULL, labels=NULL, expand = c(0, 0), limits = c(0, NA)) +
        geom_vline(xintercept=boxes$x1, alpha=0.15, linetype='F1') +
        facet_grid(strain~.) +
        labs(x='S288C Genome Coordinate', y='Aligned Chromosome Coordinate', color='Aligned Chromosome') +
        theme(legend.position = 'bottom') +
        guides(colour = guide_legend(override.aes = list(size = 3)))
    return(g)
}

g.all <- plot_all()
ggsave(g.all, file='reference_alignment/all-aligned-S288C.png', width=40, height=40, units='cm')
ggsave(g.all, file='reference_alignment/all-aligned-S288C.svg', width=40, height=40, units='cm')

quit()































############ STOP







setcolorder(dat, c('ref_chr','ref_start','ref_stop','size',
                   'query_chr','query_start','query_stop',
                   'SV_type'))

SVs <- dat[, .N, by=list(strain, ref_chr, SV_type)]


chrLvls <- c(
'chrI','chrII','chrIII','chrIV','chrV',
'chrVI','chrVII','chrVIII','chrIX','chrX',
'chrXI','chrXII','chrXIII','chrXIV','chrXV',
'chrXVI','mitochondrion')

dat[, ref_chr := factor(ref_chr, levels=chrLvls)]

ggplot(data=SVs, aes(x=ref_chr, y=N, fill=SV_type)) +
    geom_bar(stat='identity') +
    facet_grid(strain~.) +
    labs(x='Chromosome', y='Count', fill='Structural Variant') +
    theme_few(12)