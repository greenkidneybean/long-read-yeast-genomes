library(data.table)
library(ggplot2)
library(foreach)
library(doMC)
library(RColorBrewer)
library(zoo)
library(ggthemes)

registerDoMC(cores=2)

args <- commandArgs(trailingOnly=TRUE)
parent1 <- args[1]
parent2 <- args[2]

get_params_from_filename <- function(blast_file) {
    filename <- unlist(strsplit(blast_file, split='.blast.gz'))[1]
    filestem <- unlist(strsplit(filename, split='data/blasts/'))[2]
    filestem_split <- unlist(strsplit(filestem, split='_'))
    cross <- filestem_split[1]
    plate <- filestem_split[2]
    well <- filestem_split[3]
    return(c(cross,plate,well))
}

format_outfile <- function(cross, plate, well) {
    f <- paste0("data/imputed/", cross, "_", plate, "_", well, ".txt")
    return(f)
}

format_outimage <- function(cross, plate, well) {
    f <- paste0("data/imputed/", cross, "_", plate, "_", well, ".png")
    return(f)
}

import_blast <- function(blast_filename) {
    dt <- fread(blast_filename, header=F)
    setnames(dt, c("readID", "parent_chr", "PID", "match_length", "V5", "V6", 
                    "V7", "V8", "match_start", "match_end", "evalue", "bits"))
    dt[, c("V5", "V6", "V7", "V8") := NULL]
    dt[, c("strain_tmp", "chr") := tstrsplit(parent_chr, split="_")]
    dt[strain_tmp==parent1, strain := 'genotype1']
    dt[strain_tmp==parent2, strain := 'genotype2']
    dt[, strain_tmp := NULL]
    dt[, parent_chr := NULL]
    return(dt[])
}

filter_blast <- function(dt, min_PID=100, min_evalue=1e-60, min_match_length=130) {
    dt <- dt[evalue <= min_evalue & match_length >= min_match_length & PID >= min_PID][]
    dt[, .SD[evalue == min(evalue)], by=readID]
}

convert_to_wide <- function(dt, window_size = 35000) {
    dt[, window := window_size * floor(match_start / window_size)]
    dt.counts <- dt[, .N, by=list(chr, strain, window)]
    dcast(dt.counts, window+chr~strain, value.var='N', fill=0)
}

binom_test <- function(x) {
    binom.test(x, p=0.5, alternative="two.sided")$p.value
}

calculate_log_odds <- function(dt) {
    dt[, logOdds := log((genotype1 + 1) / (genotype2 + 1))]
}

calculate_binomial_p <- function(dt, parent1, parent2) {
    dt[, p := sapply(Map(c, dt$genotype1, dt$genotype2), binom_test)]
}

fill_empty_windows <- function(dt, window_size = 35000) {
    tmpTable <- foreach(chr.i = unique(dt[,chr]), .combine='rbind') %do% {
        chr.min <- 0
        chr.max <- max(dt[chr==chr.i, window])
        vals <- seq(0, chr.max, window_size)
        data.table('window' = vals, 'chr' = chr.i)
    }
    setkey(tmpTable, chr, window)
    setkey(dt, chr, window)
    return(merge(dt, tmpTable, all=TRUE))
}

fill_haplotypes <- function(dt) {
    dt[logOdds > 0, haplotype := parent1]
    dt[logOdds < 0, haplotype := parent2]
    dt[p > 0.05, haplotype := NA]
    dt[, nextHap := na.locf(haplotype, na.rm=FALSE)]
    dt[, prevHap := na.locf(haplotype, na.rm=FALSE, fromLast=TRUE)]
    dt[window == min(window), prevHap := NA, by=chr]
    dt[window == max(window), nextHap := NA, by=chr]

    # if NA but previous and next haplotype are the same, fill:
    dt[is.na(haplotype) & nextHap == prevHap, haplotype := nextHap]

    # if NA, fill with adjacent non-NA haplotype
    dt[is.na(haplotype) & is.na(prevHap), haplotype := nextHap]
    dt[is.na(haplotype) & is.na(nextHap), haplotype := prevHap]
    dt[, prevHap := NULL]
    dt[, nextHap := NULL]
    dt[is.na(haplotype), haplotype := 'indeterminate']
    return(dt[])
}

collapse_windows <- function(dt, window_size = 35000) {
    setkey(dt, chr, window)
    dt[, id := rleid(haplotype), by=chr]

    haplotypeMap <- dt[, list(haplotype, 'start'=min(window)+1, 'stop' = max(window)+window_size), by=list(id,chr)]
    haplotypeMap[, id := NULL]
    haplotypeMap <- haplotypeMap[!duplicated(haplotypeMap)]
    haplotypeMap[, chr := factor(chr, levels=c('chrI','chrII','chrIII','chrIV','chrV','chrVI',
                                                'chrVII','chrVIII','chrIX','chrX','chrXI','chrXII',
                                                'chrXIII','chrXIV','chrXV','chrXVI','mitochondrion'))]

    haplotypeMap[, chrID := as.numeric((chr))]
    haplotypeMap[is.na(haplotype), haplotype := "Indeterminate"]
    setkey(haplotypeMap, chrID, start, stop)
    return(haplotypeMap[])
}



plot_haplotypes <- function(dt) {
    ggplot(data=dt, aes(xmin=start, xmax=stop, ymin = 0, ymax=1, fill=haplotype)) +
    geom_rect() +
    facet_grid(chr~.) +
    theme_few() +
    theme(axis.ticks.y = element_blank(),axis.text.y = element_blank()) +
    theme(strip.text.y = element_text(angle=0, hjust=0))
}


process_sample <- function(filename) {
    blast <- convert_to_wide(filter_blast(import_blast(filename)))

    params <- get_params_from_filename(filename)
    cross <- params[1]
    plate <- params[2]
    well  <- params[3]

    outfile_name <- format_outfile(cross, plate, well)
    outimage_name <- format_outimage(cross, plate, well)
    calculate_log_odds(blast)
    calculate_binomial_p(blast)
    blast <- fill_empty_windows(blast)
    blast <- fill_haplotypes(blast)

    haplotypes <- collapse_windows(blast)
    haplotypes[, 'cross' := cross]
    haplotypes[, 'plate' := plate]
    haplotypes[, 'well' := well]
    ggsave(plot_haplotypes(haplotypes), file=outimage_name, width=15, height=15, units='cm')
    fwrite(haplotypes, file=outfile_name, quote=F, row.names=F, col.names=T, sep='\t')
}



main <- function() {
    blast_filenames <- list.files('data/blasts', full.names=TRUE)
    for(i in blast_filenames) {
        print(i)
        try(process_sample(i))
    }
}


main()