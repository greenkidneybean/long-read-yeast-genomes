#!/usr/bin/env bash

module load bwa
module load samtools
#module load R/3.6.3

fastq_file=${1}
ref_strain=${2}

sample=$(basename ${fastq_file%.fastq})

echo $sample


# index parental ref_strain 1
if [ ! -f "data/pacbio/${ref_strain}.fasta.bwt" ]; then
    bwa index "data/pacbio/${ref_strain}.fasta"
fi


mkdir -p data/{bam,pileup}

function bwa_mapping() {
    sample=${1}
    ref_strain=${2}

    bwa mem data/pacbio/${ref_strain}.fasta data/fastq/${sample}.fastq | \
        samtools view -hb - | samtools sort - > data/bam/${sample}_${ref_strain}.bam
    
    bwa index data/bam/${sample}_${ref_strain}.bam
}

export -f bwa_mapping

parallel -j 1 bwa_mapping ::: ${sample} ::: ${ref_strain}
