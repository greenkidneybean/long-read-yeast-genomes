#!/usr/bin/env bash
#SBATCH --ntasks=1
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=120G
#SBATCH --time=6:00:00
#SBATCH --gres=lscratch:200


WORKDIR="/data/SBGE/cory/pacbio/"
RUN=${1}
TMPDIR="/lscratch/${SLURM_JOB_ID}/"


cd $TMPDIR

function get_samfile_list() {
    local RUN=${1}
    local samfiles=$(zip -sf ${WORKDIR}/${RUN}.sam.zip | tail -n +2 | head -n -1)
    echo ${samfiles}
}

export -f get_samfile_list

function build_fastq() {
    local RUN=${1}
    local ZIP_FILENAME=${RUN}.sam.zip
    local SAM_FILENAME=${2}
    # local FILESTEM=$(echo ${SAM_FILENAME%.sam} | sed 's/^A/A_/g')
    local FILESTEM=$(echo ${SAM_FILENAME%.sam})
    local CROSS=$(echo $FILESTEM | cut -d '_' -f 1)
    local PLATE=$(echo $FILESTEM | cut -d '_' -f 2)
    local WELL=$(echo $FILESTEM | cut -d '_' -f 3)
    samtools fastq \
    -1 ${RUN}_${PLATE}_${WELL}.fastq \
    -2 ${RUN}_${PLATE}_${WELL}.fastq \
    -0 ${RUN}_${PLATE}_${WELL}.fastq \
    -s ${RUN}_${PLATE}_${WELL}.fastq \
    <(cat ${WORKDIR}/${RUN}.header.txt <(unzip -p  ${WORKDIR}/${ZIP_FILENAME} ${SAM_FILENAME}))
}

export -f build_fastq

# function build_header() {
#     local RUN=${1}
#     local GZ_FILE=${WORKDIR}/${RUN}.sam.gz
#     local ZIP_FILE=${WORKDIR}/${RUN}.sam.zip
#     # get cross ID
#     local CROSS=$(zip -sf ${ZIP_FILE} | head -n 2 | tail -n 1 | xargs echo -n | cut -d "_" -f 1)
#     # Extract full header
#     zcat ${GZ_FILE} | head -n 30000 | grep "^@" > ${RUN}.header.full.txt

#     # Include all @HD lines
#     grep "^@HD" ${RUN}.header.full.txt > ${RUN}.header.txt

#     # Include all @SQ lines
#     grep "^@SQ" ${RUN}.header.full.txt >> ${RUN}.header.txt

#     # Include all @RG lines that include $RUN
#     grep "^@RG" ${RUN}.header.full.txt | grep "ID:$CROSS"  >> ${RUN}.header.txt

#     # Include all @PG lines that include $RUN
#     grep "^@PG" ${RUN}.header.full.txt | grep "$CROSS" >> ${RUN}.header.txt

#     # remove initial tmp header
#     rm ${RUN}.header.full.txt
# }

# export -f build_header

echo "beginning run ${RUN}"
# echo 'building header...'
# build_header ${RUN}

module load samtools

samfile_list=($(get_samfile_list $RUN))

echo 'generating fastqs'
for samfile in ${samfile_list[@]}; do
    echo $samfile
    build_fastq ${RUN} ${samfile}
done

echo 'building zip'
ls ${RUN}*.fastq | zip -@ -j ${RUN}.fastq.zip

echo "copying ${RUN}.fastq.zip to ${WORKDIR}/"
cp ${RUN}.fastq.zip ${WORKDIR}/
