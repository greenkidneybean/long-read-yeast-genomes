#!/usr/bin/env bash
#SBATCH --ntasks=1
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --time=1:00:00

RUN="SRR9330831"
STRAIN1="CBS2888"
STRAIN2="CLIB219"
CHR="chrXV"

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
        --targets-file targets.txt \
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


# define temporary folder on ramdisk
TMP="/dev/shm/${USER}/${SLURM_JOB_ID}/"
mkdir -p ${TMP}

# move requisite files to ramdisk
cp /data/SBGE/cory/pacbio/SRR9330831.fastq.zip ${TMP}
cp /data/SBGE/cory/pacbio/CBS2888_chrXV.fasta ${TMP}
cp /data/SBGE/cory/pacbio/CLIB219_chrX.fasta ${TMP}
cp /data/SBGE/cory/pacbio/CLIB219_chrXV.fasta ${TMP}


cd ${TMP}

cat CLIB219_chrX.fasta CLIB219_chrXV.fasta > CLIB219_merged.fasta

module load mummer
module load python/3.9

nucmer CBS2888_chrXV.fasta CLIB219_merged.fasta && \
dnadiff -d out.delta

awk '$11~"CBS2888_pacbio_chrXV$"' out.snps > out.fixed.snps
rm out.snps
mv out.fixed.snps out.snps

python ~/pacbio-yeast-genomes/src/filter_snps.py out.snps 20 > mask.filter.bed

python ~/pacbio-yeast-genomes/src/print_target_snps.py out.snps 20 > targets.txt

python ~/pacbio-yeast-genomes/src/mask_fasta.py CBS2888_chrXV.fasta mask.filter.bed --out CBS2888_CLIB219_chrXV.mask.fasta

bwa index CBS2888_CLIB219_chrXV.mask.fasta

parallel -j 8 map_and_call {} \
    ::: ${RUN} \
    ::: ${STRAIN1} \
    ::: ${STRAIN2} \
    ::: ${CHR} \
    ::: {G,R}{1,2,3,4,5} \
    ::: {01..96}

# combine and export 
ls *.call | tar -czv --files-from - -f /data/SBGE/cory/pacbio/${RUN}_${STRAIN1}_${STRAIN2}_${CHR}translocation2.call.tar.gz

cd && rm -rf "${TMP}"