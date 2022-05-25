#!/usr/bin/env bash

module load bwa
#module load R/3.6.3

fastq_file=${1}
ref_strain1=${2}
ref_strain2=${3}

sample=$(basename ${fastq_file%.fastq})

echo $sample


# index parental ref_strain 1
if [ ! -f "data/pacbio/pacbio_${ref_strain1}.fasta.bwt" ]; then
    bwa index "data/pacbio/pacbio_${ref_strain1}.fasta"
fi

# index parental ref_strain 2
if [ ! -f "data/pacbio/pacbio_${ref_strain2}.fasta.bwt" ]; then
    bwa index "data/pacbio/pacbio_${ref_strain2}.fasta"
fi

mkdir -p data/{bam,pileup}

function bwa_mapping() {
    sample=${1}
    ref_strain=${2}

    bwa mem data/pacbio/pacbio_${ref_strain}.fasta data/fastq/${sample}.fastq | \
        samtools view -hb - | samtools sort - > data/bam/${sample}_${ref_strain}.bam
    
    bwa index data/bam/${sample}_${ref_strain}.bam
}

export -f bwa_mapping

parallel -j 1 bwa_mapping ::: ${sample} ::: ${ref_strain1} ${ref_strain2}






# Rscript - <<EOF
#     library(data.table)
#     dat <- fread('test_273614.pileup')
#     dat[, V5 := NULL]
#     dat[, V6 := NULL]
#     setnames(dat, c("chr","pos","ref","count"))

#     dat[, window := cut(pos, breaks=seq(0,10e6, 5000))]
#     dat[, window := as.numeric(window)]

#     counts <- dat[, list("N"=sum(count)), by=list(chr, window)]

#     print(counts)

# EOF