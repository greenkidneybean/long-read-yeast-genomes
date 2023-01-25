#!/usr/bin/env Rscript

library(data.table)
library(optparse)
library(ggthemes)
library(foreach)
library(doMC)
library(parallel)
library(cowplot)
registerDoMC(cores=2)

# minimum segment length to plot
minLen <- 10000


# library(optparse)
# parser <- add_option(parser, 
#                         c("-r", "--reference"), 
#                         action="store_true", 
#                         default=TRUE,
#                         help="Print extra output [default]"
# )


import_paf <- function(filename) {

    dat <- fread(filename, sep='\t', fill=TRUE)
    # Note A=Query, B=Ref
    dat <- dat[,1:9]
    setnames(dat, c('seqA','lenA','startA','stopA','strand',
                    'seqB','lenB','startB','stopB'))

    dat[, c('strainA','chrA') := tstrsplit(seqA, split='_')]
    dat[, c('tmpA','tmpB') := tstrsplit(seqB, split='\\|')]
    dat[, tmpA := NULL]
    dat[, seqB := NULL]
    setnames(dat, 'tmpB', 'NC')

    library(ggplot2)
    library(ggrepel)

    dat[, spanA := 1 + stopA-startA]
    dat[, spanB := 1 + stopB-startB]
    return(dat[])
}


paf_files <- list.files(path='reference_alignment/', pattern='*.paf', full.names=TRUE)
chr_id_file <- 'chromosome_ids.txt'

print(chr_id_file)
print(paf_files)


dat.all <- foreach(f=paf_files, .combine='rbind') %do% {
    import_paf(f)
}


chrLvls <- c(
'chrI','chrII','chrIII','chrIV','chrV',
'chrVI','chrVII','chrVIII','chrIX','chrX',
'chrXI','chrXII','chrXIII','chrXIV','chrXV',
'chrXVI','mitochondrion')


id_table <- fread(chr_id_file, header=F)
setnames(id_table, c('NC', 'chrB'))
setkey(id_table, NC)
setkey(dat.all, NC)
dat.merge <- merge(dat.all, id_table)
dat.merge[, chrA := factor(chrA, levels=chrLvls)]
dat.merge[, chrB := factor(chrB, levels=chrLvls)]


dat.merge[, sequentialID := as.numeric(chrB)]

chrLengths <- unique(dat.merge[, c('lenB', 'sequentialID')])
chrLengths[, cumulativeLength := cumsum(lenB)]
chrLengths[, cumulativeLength := shift(cumulativeLength, n=1L, type='lag')]
chrLengths[is.na(cumulativeLength), cumulativeLength := 0]
chrLengths[, lenB := NULL]
setkey(chrLengths, sequentialID)
setkey(dat.merge, sequentialID)

dat <- merge(chrLengths, dat.merge)

dat[,startB := startB + cumulativeLength]
dat[,stopB := stopB + cumulativeLength]
dat[, colorBackground := sequentialID%%2]   # odds -> 1, evens -> 0



boxes <- data.table(x1=c(0,unique(dat[cumulativeLength != 0, cumulativeLength])))
boxes[, x2 := shift(x1, n=1L, type='lead', fill=max(dat$stopB))]
boxes[, midpoint := (x1+x2)/2]
boxes[, chr := chrLvls]
boxes[chr=='mitochondrion', chr := 'm']

plot_all <- function() {
    g <- ggplot(data=dat[spanA >= minLen & spanB >= minLen]) + 
        geom_segment(aes(x=startB, xend=stopB, y=startA, yend=stopA, color=chrA)) +
        theme_few() +
        scale_x_continuous(breaks=boxes$midpoint, labels=boxes$chr, expand = c(0, 0), limits = c(0, NA)) +
        scale_y_continuous(breaks=NULL, labels=NULL, expand = c(0, 0), limits = c(0, NA)) +
        geom_vline(xintercept=boxes$x1, alpha=0.15, linetype='F1') +
        facet_grid(strainA~.) +
        labs(x='S288C Genome Coordinate', y='Aligned Chromosome Coordinate', color='Aligned Chromosome') +
        theme(legend.position = 'bottom') +
        guides(colour = guide_legend(override.aes = list(size = 3)))
    return(g)
}

g.all <- plot_all()
ggsave(g.all, file='reference_alignment/all.png', width=40, height=40, units='cm')

quit()

