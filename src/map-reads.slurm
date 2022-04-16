#!/usr/bin/env bash


function makeFasta() {
    cross=${1}
    plate=${2}
    well=${3}
    repoPath=${4}
    sample="${plate}_${well}"
    reference=${5}

    samtools fastq \
    -1 ${cross}_${sample}.fastq \
    -2 ${cross}_${sample}.fastq \
    -0 ${cross}_${sample}.fastq \
    -s ${cross}_${sample}.fastq \
    <(cat ${repoPath}/data/header-${cross}.txt <(unzip -p ${repoPath}/data/${cross}.sam.zip ${cross}_${sample}.sam))

    src/bwa.sif bwa index data/MSY37-36-chrIV.fasta
    src/bwa.sif bwa mem data/MSY37-36-chrIV.fasta  ${cross}_${sample}.fastq | \
    src/bwa.sif samtools view -b | src/bwa.sif samtools sort -o ${cross}_${sample}.bam -
    src/bwa.sif samtools index ${cross}_${sample}.bam

}
export -f makeFasta

makeFasta $1 $2 $3 $PWD $5

src/minimap.sif minimap2 -ax sr -B 25 data/MSY37-36-chrIV.fasta 3003_G1_01.fastq > mapping.sam
samtools sort -o mapping.bam mapping.sam
samtools depth mapping.bam > mapping.depth

singularity exec ../pseudodiploidy/src/pseudodiploidy.sif R              

library(ggplot2)
library(data.table)
dat <- fread('mapping.depth')
setnames(dat, c('contig', 'pos', 'depth'))

dat[, posBin := cut(pos, breaks=seq(0,1.6e6, 1000))]
dat[, posBin := as.numeric(posBin)]
dat.ag <- dat[, list('medianDepth'=median(depth)), by=list(contig,posBin)]
dat.ag[contig == 'tig00000001', medianDepth := -1 * medianDepth]
ggplot(dat.ag, aes(x=posBin, y=medianDepth, color=contig)) + geom_point()
