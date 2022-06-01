#!/usr/bin/env bash

module load bwa
module load samtools
module load bcftools
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
    if [ ! -f data/bam/${sample}_${ref_strain}.bam ]; then
        bwa mem data/pacbio/${ref_strain}.fasta data/fastq/${sample}.fastq | \
            samtools view -hb - | samtools sort - > data/bam/${sample}_${ref_strain}.bam
    fi

    if [ ! -f data/bam/${sample}_${ref_strain}.bam.bai ]; then
        samtools index data/bam/${sample}_${ref_strain}.bam
    fi
    
    if [ ! -f data/vcf/${sample}.vcf ]; then
        bcftools mpileup \
            --targets-file data/pacbio/${ref_strain}.targets.txt \
            --fasta-ref data/pacbio/${ref_strain}.fasta \
            data/bam/${sample}_${ref_strain}.bam | \
            bcftools call \
            --ploidy 1 -m -Ob | \
            bcftools view | \
            sed 's/1:.*$/1/g' | \
            grep -v "^##" > data/vcf/${sample}.vcf
    fi



}

export -f bwa_mapping

parallel -j 1 bwa_mapping ::: ${sample} ::: ${ref_strain}
