#!/usr/bin/env bash
#SBATCH --ntasks=1
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --time=3:59:59
#SBATCH --gres=lscratch:200
#SBATCH --partition=quick,norm

export STRAIN1=${1}
export STRAIN2=${2}
export CHR1=${3}
export CHR2=${4}
export SCRATCH="/lscratch/${SLURM_JOB_ID}/"
export DATADIR="/data/SBGE/cory/pacbio/"
export GITDIR="/home/wellerca/pacbio-yeast-genomes/"
cd ${SCRATCH}

get_sra_run() {
    s1=${1}
    s2=${2}

    matchedrun=$(cat ${GITDIR}/runs.txt | grep ${STRAIN1} | grep ${STRAIN2} | awk '{print $4}')
    echo ${matchedrun}    
}

# determine SRA RUN ID based on STRAIN1 and STRAIN2 (see `runs.txt`)
export RUN=$(get_sra_run ${STRAIN1} ${STRAIN2})

echo "STRAIN1 = ${STRAIN1}"
echo "STRAIN2 = ${STRAIN2}"
echo "CHR1 = ${CHR1}"
echo "CHR2 = ${CHR2}"
echo "RUN = ${RUN}"

# load modules
module load bwa/0.7.17
module load samtools/1.16.1
module load mummer/4.0.0beta2
module load python/3.9

# move ref fasta seqs to local scratch
unzip ${DATADIR}/pacbio-assemblies.zip ${STRAIN1}_${CHR1}.fasta
unzip ${DATADIR}/pacbio-assemblies.zip ${STRAIN2}_${CHR1}.fasta
unzip ${DATADIR}/pacbio-assemblies.zip ${STRAIN2}_${CHR2}.fasta

# move fastq zip to local scratch
cp ${DATADIR}/${RUN}.fastq.zip .

map_and_call() {
    local plate=${1}
    local well=${2}
    local cross="${STRAIN1}_${STRAIN2}"
    local id="${RUN}_${plate}_${well}"

    bwa mem ${REF_FASTA} <((unzip -p ${RUN}.fastq.zip ${id}.fastq)) | \
        samtools view -hb - | \
        samtools sort - > ${id}.bam

    samtools index ${id}.bam
    bcftools mpileup \
        --targets-file targets.txt \
        --fasta-ref ${REF_FASTA} \
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

# merge both chromosomes for STRAIN2 into a single fasta
cat ${STRAIN2}_${CHR1}.fasta ${STRAIN2}_${CHR2}.fasta > ${STRAIN2}_merged.fasta

# perform CHR1 to (merged CHR1+CHR2) alignment
nucmer ${STRAIN1}_${CHR1}.fasta ${STRAIN2}_merged.fasta && \
dnadiff -d out.delta

# ensure everything from perspective of first strain
awk -v s1=${STRAIN1} '$11~s1' out.snps > out.fixed.snps
rm out.snps
mv out.fixed.snps out.snps

# generate filters
python ${GITDIR}/src/filter_snps.py out.snps 20 > mask.filter.bed
python ${GITDIR}/src/print_target_snps.py out.snps 20 > targets.txt

# mask reference fasta before mapping
export REF_FASTA="${STRAIN1}_${STRAIN2}_${CHR1}.mask.fasta"
python ${GITDIR}/src/mask_fasta.py ${STRAIN1}_${CHR1}.fasta mask.filter.bed --out ${REF_FASTA}

# index new masked reference
bwa index ${REF_FASTA}

# do the actual mapping/genotyping
parallel -j 8 map_and_call {} \
    ::: {G,R}{1,2,3,4,5} \
    ::: {01..96}

# combine and export 
ls *.call | tar -czv --files-from - -f ${DATADIR}/${STRAIN1}_${STRAIN2}_${CHR1}_${CHR2}_translocation.call.tar.gz

cd