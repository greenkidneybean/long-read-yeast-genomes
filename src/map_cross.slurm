#!/usr/bin/env bash
#SBATCH --ntasks=1
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --time=1:00:00

RUN=${1}
STRAIN1=${2}
STRAIN2=${3}
CHR=${4}

run_mummer() {
    module load mummer
    local strain1=${1}
    local strain2=${2}
    local chr=${3}
    local cross="${strain1}_${strain2}"
    nucmer -p ${cross}_${chr} ${strain1}_${chr}.fasta ${strain2}_${chr}.fasta && \
    dnadiff -p ${cross}_${chr} -d ${cross}_${chr}.delta
}
export -f run_mummer

filter_snps() {
    local strain1=${1}
    local strain2=${2}
    local chr=${3}
    local cross="${strain1}_${strain2}"
   python ~/pacbio-yeast-genomes/src/filter_snps.py ${cross}_${chr}.snps 20 > ${cross}_${chr}.mask.filter.bed
}
export -f filter_snps

print_vcf_targets() {
    local strain1=${1}
    local strain2=${2}
    local chr=${3}
    local cross="${strain1}_${strain2}"
    local filter_dist='20'
    local input="${cross}_${chr}.snps"
    local output="${cross}_${chr}.targets.txt"
    python ~/pacbio-yeast-genomes/src/print_target_snps.py ${input} ${filter_dist} > ${output}
}
export -f print_vcf_targets

mask_ref_fasta() {
    local strain1=${1}
    local strain2=${2}
    local chr=${3}
    local cross="${strain1}_${strain2}"
    local input_fasta="${strain1}_${chr}.fasta"
    local input_bed="${cross}_${chr}.mask.filter.bed"
    local output="${cross}_${chr}.mask.fasta"
    python ~/pacbio-yeast-genomes/src/mask_fasta.py ${input_fasta} ${input_bed} --out ${output}
}
export -f mask_ref_fasta

index_masked_ref_fasta() {
    module load bwa
    local strain1=${1}
    local strain2=${2}
    local chr=${3}
    local cross="${strain1}_${strain2}"
    local input="${cross}_${chr}.mask.fasta"
    local output="${cross}_${chr}.mask.fasta.bwt"
    bwa index ${input}
    # check output?
}
export -f index_masked_ref_fasta

map_and_call() {
    module load bwa
    module load samtools
    local run=${1}
    local strain1=${2}
    local strain2=${3}
    local chr=${4}
    local cross="${strain1}_${strain2}"
    local ref_fasta=${strain1}_${strain2}_${chr}.mask.fasta
    local plate=${5}
    local well=${6}
    local id="${run}_${plate}_${well}"

    bwa mem ${ref_fasta} <((unzip -p ${run}.fastq.zip ${id}.fastq)) | \
        samtools view -hb - | \
        samtools sort - > ${id}.bam

    samtools index ${id}.bam
    bcftools mpileup \
        --targets-file ${cross}_${chr}.targets.txt \
        --fasta-ref ${ref_fasta} \
        ${id}.bam | \
        bcftools call \
        --ploidy 1 -m -Ob | \
        bcftools view | \
        sed 's/1:.*$/1/g' | \
        grep -v "^##" | awk '{print $1,$2,$4,$5,$10}' | sed 's/\.bam$//g' > ${id}.call
    
    rm ${id}.bam && \
    rm ${id}.bam.bai
}
export -f map_and_call

prep_for_bwa() {
    local strain1=${1}
    local strain2=${2}
    local chr=${3}
    run_mummer ${strain1} ${strain2} ${chr}
    filter_snps ${strain1} ${strain2} ${chr}
    print_vcf_targets ${strain1} ${strain2} ${chr}
    mask_ref_fasta ${strain1} ${strain2} ${chr}
    index_masked_ref_fasta ${strain1} ${strain2} ${chr}
}

# define temporary folder on ramdisk
TMP="/dev/shm/${USER}/${SLURM_JOB_ID}/"
mkdir -p ${TMP}

# move requisite files to ramdisk
cp /data/SBGE/cory/pacbio/${RUN}.fastq.zip ${TMP}
cp /data/SBGE/cory/pacbio/${STRAIN1}_${CHR}.fasta ${TMP}
cp /data/SBGE/cory/pacbio/${STRAIN2}_${CHR}.fasta ${TMP}

cd ${TMP}

prep_for_bwa ${STRAIN1} ${STRAIN2} ${CHR}

parallel -j 8 map_and_call {} \
    ::: ${RUN} \
    ::: ${STRAIN1} \
    ::: ${STRAIN2} \
    ::: ${CHR} \
    ::: {G,R}{1,2,3,4,5} \
    ::: {01..96}

# combine and export 
ls *.call | tar -czv --files-from - -f /data/SBGE/cory/pacbio/${RUN}_${STRAIN1}_${STRAIN2}_${CHR}.call.tar.gz

cd && rm -rf "${TMP}"