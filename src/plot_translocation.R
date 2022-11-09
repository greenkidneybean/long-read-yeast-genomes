#!/usr/bin/env Rscript

library(data.table)
library(foreach)
library(ggplot2)
library(ggthemes)
library(viridis)
library(cowplot)

args <- commandArgs(trailingOnly=TRUE)
ld_file <- args[1]

#ld_files <- list.files('../data/output/', pattern='*_translocation.ld.txt', full.names=TRUE)

plot_translocation <- function(f) {
    dat.ld <- fread(f)
    f2 <- strsplit(strsplit(basename(f),'\\.')[[1]][1], split="_")[[1]]
    strain1 <- f2[1]
    strain2 <- f2[2]
    chrom1 <- f2[3]
    chrom2 <- f2[4]

    dat.ld[CHR_A==1, CHR1 := chrom1]
    dat.ld[CHR_A==2, CHR1 := chrom2]
    dat.ld[CHR_B==1, CHR2 := chrom1]
    dat.ld[CHR_B==2, CHR2 := chrom2]
    dat.ld[, c("CHR_A","CHR_B","SNP_A","SNP_B") := NULL]
    # plot X to X, XV to XV, X to XV
    g.out <- ggplot(dat.ld, aes(x=BP_A, y=BP_B, fill=R2)) + geom_tile() +
    theme(panel.background = element_blank(),
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank()
            ) +
    theme_few() +
    scale_fill_viridis() +
    labs(   x='Distance along chromosome (100 kb)',
            y='Distance along chromosome (100 kb)',
            title=paste0(strain1, ' x ', strain2)) + 
    scale_x_continuous(labels = seq(0, 50, 1), breaks=seq(0, 5e6, 1e5)) +
    scale_y_continuous(labels = seq(0, 50, 1), breaks=seq(0, 5e6, 1e5)) +
    theme(legend.position = c(0.9, 0.2)) +
    facet_grid(CHR2~CHR1)

    ggsave(g.out, 
            file=paste0(strain1, '_', 
            strain2, '_', 
            chrom1, '_',
            chrom2, '.pdf'),
            units='cm',
            width=32,
            height=30,
            dpi=300)
}

plot_translocation(ld_file)
