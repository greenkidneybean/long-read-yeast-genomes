library(data.table)
library(ggthemes)

args <- commandArgs(trailingOnly=TRUE)
pafFile <- args[1]

dat <- fread(pafFile, header=FALSE, fill=TRUE)
dat <- dat[,1:9]

setnames(dat, c('seqA','lenA','startA','stopA','strand',
                'seqB','lenB','startB','stopB'))

dat[, c('strainA','chrA') := tstrsplit(seqA, split='_')]
dat[, c('strainB','chrB') := tstrsplit(seqB, split='_')]

chrLvls <- c(
'chrI','chrII','chrIII','chrIV','chrV',
'chrVI','chrVII','chrVIII','chrIX','chrX',
'chrXI','chrXII','chrXIII','chrXIV','chrXV',
'chrXVI','mitochondrion')

dat[, chrA := factor(chrA, levels=chrLvls)]
dat[, chrB := factor(chrB, levels=chrLvls)]

library(ggplot2)
library(ggrepel)

dat[, spanA := 1 + stopA-startA]
dat[, spanB := 1 + stopB-startB]

strainA <- unique(dat$strainA)
strainB <- unique(dat$strainB)

minLen <- 10000
ggplot(dat[spanA >= minLen & spanB >= minLen], aes(x=startA, xend=stopA, y=startB, yend=stopB, color=chrB, label=chrB)) +
geom_segment() +
geom_text(nudge_x=1e5) +
facet_wrap(~seqA, scales='free') +
theme_few() +
labs(x=strainA, y=strainB)