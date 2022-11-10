#!/usr/bin/env Rscript

library(data.table)
library(foreach)
library(ggplot2)
library(ggthemes)
library(viridis)
library(cowplot)

args <- commandArgs(trailingOnly=TRUE)

assoc_file <- args[1]
ld_file <- args[2]

plot_fig <- function(assoc_filename, ld_filename) {
    dat.assoc <- fread(assoc_filename)
    f2 <- strsplit(strsplit(basename(assoc_filename),'\\.')[[1]][1], split="_")[[1]]
    strain1 <- f2[1]
    gene <- f2[2]
    strain2 <- f2[3]
    chromosome <- f2[4]

    dat.ld <- fread(ld_filename)

    max1 <- max(dat.assoc[,BP])
    max2 <- max(dat.ld[,BP_B])
    overallmax <- max(max1, max2)

    g.bar <- ggplot(dat.assoc, aes(x=BP, y=-1*log10(P))) + geom_bar(stat='identity',position='dodge') +
    theme_few() +
    labs(title=paste0(strain1, ' x ', strain2, ' ', chromosome, ', association to ', gene), x='') +
    scale_x_continuous(labels = seq(0, 50, 1), breaks=seq(0, 5e6, 1e5), limits=c(0,overallmax))

    g.ld <- ggplot(dat.ld, aes(x=BP_A, y=BP_B, fill=R2)) + geom_tile() +
    theme(panel.background = element_blank(),
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank()
            ) +
    theme_few() +
    scale_fill_viridis() +
    labs(   x='Distance along chromosome (100 kb)',
            y='Distance along chromosome (100 kb)') + 
    scale_x_continuous(labels = seq(0, 50, 1), breaks=seq(0, 5e6, 1e5), limits=c(0,overallmax)) +
    scale_y_continuous(labels = seq(0, 50, 1), breaks=seq(0, 5e6, 1e5), limits=c(0,overallmax)) +
    theme(legend.position = c(0.9, 0.2))

    g.out <- plot_grid(g.bar, g.ld, align = "hv", nrow=2, rel_heights=c(0.2, 1))
    ggsave(g.out, 
            file=paste0('plots/', strain1, '_', gene, '_', strain2, '_', chromosome, '.pdf'),
            units='cm',
            width=25,
            height=30,
            dpi=300)
}

plot_fig(assoc_file, ld_file)
