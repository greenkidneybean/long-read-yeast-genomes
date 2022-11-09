#!/usr/bin/env R


library(data.table)
library(ggplot2)
library(foreach)
library(ggthemes)

files <- list.files(pattern='*.counts.txt')

dat <- foreach(file=files, .combine='rbind') %do% {
    fread(file)
}
~/pacbio-yeast-genomes/SRR9330831_39-tf_G1.counts.txt
# dat <- fread('SRR9330831_39-tf_R1.counts.txt')
setnames(dat, c('sample_gene','length','count','count2'))

dat[, c('sample','gene') := tstrsplit(sample_gene, split=' ')]
dat[, sample_gene := NULL]
setcolorder(dat, c('sample', 'gene','length','count','count2'))
dat[gene=='*', gene := 'unmapped']
dat[gene=='unmapped', count := count2]
dat[, count2 := NULL]

dat[, total := sum(count), by=sample]
dat[, RPM := 1e6 * count / total ]
dat <- dat[gene != 'unmapped']

medians <- dat[, list('medRPM'=median(RPM)), by=gene]

#dummy <- data.table(gene = c('39-tf1','39-tf3','39-tf4'), Z = c(89.996172, 7.869651, 107.417401))

ggplot(data=dat, aes(x=1, y=RPM)) + geom_violin() +
facet_grid(~gene) +
geom_hline(data=medians, aes(yintercept=medRPM), alpha=0.6, linetype='dashed', color='red')

dat.merge <- merge(dat,medians)

dat.merge[RPM >= medRPM, genotype := 1]
dat.merge[RPM < medRPM, genotype := 0]
dat.merge[, gene := factor(gene, levels=c('39-tf1','39-tf3','39-tf4'))]
dat.merge[, plate := tstrsplit(sample, split='_')[[2]]]
dat.merge[, N := tstrsplit(sample, split='_')[[3]]]
dat.merge[, N := as.numeric(N)]

ggplot(dat.merge, aes(x=gene, y=N, fill=genotype)) + geom_tile(color='gray40') + facet_wrap(~plate, nrow=1) +
theme_few() +
scale_fill_manual(values=c('gray20','gray90')) +
theme(panel.border = element_blank(),
axis.text.y = element_blank(),
axis.ticks.y = element_blank()) +
guides(fill=guide_legend(title='Gene Present'))



dat.tmp <- dcast(dat.merge, sample~gene, value.var='genotype')
fwrite(dat.tmp, file='MALR.mat', quote=F, row.names=F, col.names=T, sep='\t')


dat.tmp[, combined := paste0(`39-tf1`,`39-tf3`,`39-tf4`)]
dat.tmp[, .N, by=list(combined)]

Genotypes
000 424
111 418
110  36
010   8
100  12
011  18
001  30
101  14


