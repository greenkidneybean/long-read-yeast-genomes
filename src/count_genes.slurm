#!/usr/bin/env bash
#SBATCH --ntasks=1
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=1:00:00
#SBATCH --gres=lscratch:10


module load blast
module load samtools
module load bwa

cross=${1}
gene=${2}
plate=${3}

echo CROSS ${cross}
echo GENE ${gene}
echo PLATE ${plate}


TMP="/lscratch/${SLURM_JOB_ID}"
mkdir -p ${TMP} && cd ${TMP}

cp ~/pacbio-yeast-genomes/${gene}.fasta .


if [ ! -f "${gene}.fasta.bwt" ]; then
    bwa index "${gene}.fasta"
fi

makeblastdb -dbtype nucl -in ${gene}.fasta

for sample in ${cross}_${plate}_{01..96}; do
    echo $sample
    unzip -p /data/SBGE/cory/pacbio/${cross}.fastq.zip ${sample}.fastq > ${sample}.fastq
    bwa mem -B 40 -O 60 -E 10 -L 100 ${gene}.fasta ${sample}.fastq | \
            samtools view -hb - | samtools sort -o - > ${sample}.bam
    samtools index ${sample}.bam
    samtools idxstats ${sample}.bam | awk -v s=${sample} '{print s,$0}' > ${sample}.counts.txt
    rm ${sample}.{fastq,bam,bam.bai}
done

cat *.counts.txt > /data/wellerca/projects/pacbio-yeast-genomes/${cross}_${gene}_${plate}.counts.txt

cd

