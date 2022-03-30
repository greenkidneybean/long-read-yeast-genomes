#!/usr/bin/env bash
#SBATCH --ntasks=1
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --time=1:00:00

set -e -o pipefail

export TMPDIR="/dev/shm/${USER}/${SLURM_JOB_ID}"
function finish {
    rm -rf ${TMPDIR}
}
trap finish EXIT
mkdir -p ${TMPDIR}

module load samtools

# Extract plate ID from slurm array ID
let i=${SLURM_ARRAY_TASK_ID}-1 || true      # offset array ID by 1 (bash arrays begin at 0); || true bypases exit code 0
allPlates=(G1 G2 G3 G4 G5 R1 R2 R3 R4 R5)   # array of all plates
plate=${allPlates[$i]}                      # take i-th plate

cross=${1}
evalue=${2} # e.g. 1e-10
repoPath=${3}
dbFasta=${4}

function makeFastaDoBlast() {
    cross=${1}
    plate=${2}
    well=${3}
    repoPath=${4}
    evalue=${5}
    dbFasta=${6}
    tmp=${7}
    sample="${plate}_${well}"

    mkdir -p ${tmp} && cd ${tmp}

    samtools fasta \
    -1 ${cross}_${sample}.fasta \
    -2 ${cross}_${sample}.fasta \
    -0 ${cross}_${sample}.fasta \
    -s ${cross}_${sample}.fasta \
    <(cat ${repoPath}/data/header-${cross}.txt <(unzip -p ${repoPath}/data/${cross}.sam.zip ${cross}_${sample}.sam))

    # run blastn
    ${repoPath}/blast/bin/blastn -outfmt 6 -query ${cross}_${sample}.fasta -db ${repoPath}/data/${dbFasta} | awk -v query="${cross}_${sample}" -v db="${dbFasta}" '{OFS="\t"; print query,db,$0}' >> ${tmp}/${cross}_${sample}.blast

    # count number of lines
    nReads=$(wc -l ${cross}_${sample}.fasta | awk '{print $1}')
    let nReads=${nReads}/2

    echo -e "${cross}_${sample}\t${nReads}" >> ${tmp}/${cross}_${sample}.txt

    # clean up unneeded files
    rm ${cross}_${sample}.fasta
}
export -f makeFastaDoBlast

parallel -j 8 makeFastaDoBlast ::: \
    ${cross} ::: \
    ${plate} ::: \
    $(seq -w 1 96) ::: \
    ${repoPath} ::: \
    ${evalue} ::: \
    ${dbFasta} ::: \
    ${TMPDIR}  &&
cat ${TMPDIR}/${cross}_${plate}_*.blast > ${repoPath}/data/${cross}_${plate}.blast &&
cat ${TMPDIR}/${cross}_${plate}_*.txt > ${repoPath}/data/${cross}_${plate}-summary.txt