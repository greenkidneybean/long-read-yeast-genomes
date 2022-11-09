#!/usr/bin/env bash

call_phenos() {
    local strain1=${1}
    local strain2=${2}
    local run=${3}
    local fasta=${4}
    sbatch ../src/call_phenotypes.sh \
        --strain1   ${strain1} \
        --strain2   ${strain2} \
        --fasta     fastas/${fasta} \
        --fastqs    ../data/shortreads/${run}.fastq.zip \
        --gitdir    /home/wellerca/pacbio-yeast-genomes 
}
export -f call_phenos

call_genos() {
    local strain1=${1}
    local strain2=${2}
    local chr=${3}
    local run=${4}
    sbatch ../src/call_genotypes.sh \
        --strain1   ${strain1} \
        --fasta1    ../data/assemblies/${strain1}_${chr}.fasta \
        --strain2   ${strain2} \
        --fasta2    ../data/assemblies/${strain2}_${chr}.fasta \
        --chr       ${chr} \
        --fastqs    ../data/shortreads/${run}.fastq.zip \
        --gitdir    /home/wellerca/pacbio-yeast-genomes 
}
export -f call_genos

qtl_map() {
    local strain1=${1}
    local strain2=${2}
    local fasta=${3}
    local chr=${4}
    sbatch ../src/qtl_map.sh \
        --strain1   ${strain1} \
        --strain2   ${strain2} \
        --fasta     fastas/${fasta} \
        --chr       ${chr} \
        --gitdir    /home/wellerca/pacbio-yeast-genomes 
}
export -f qtl_map

plot_association() {
    module load R/3.6.3
    local strain1=${1}
    local strain2=${2}
    local fasta=${3}
    local chr=${4}
    mkdir -p ./plots
    assoc_file=$(realpath "../data/output/${strain1}_${fasta%.fasta}_${strain2}_${chr}.qassoc.txt")
    ld_file=$(realpath "../data/output/${strain1}_${strain2}_${chr}.ld.txt")
    Rscript ../src/plot_association.R ${assoc_file} ${ld_file}
}
export -f plot_association

