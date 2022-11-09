#!/usr/bin/env bash
#SBATCH --ntasks=1
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=200G
#SBATCH --time=8:00:00
#SBATCH --gres=lscratch:20


export -f build_header
# Split large sam file contents into separate files while on biowulf ramdisk
TMPDIR="/dev/shm/${USER}/${SLURM_JOB_ID}"           # name tmp directory on /dev/shm/
mkdir -p ${TMPDIR} && cd ${TMPDIR}                  # cd into tmp directory
run=${1}
cramid=${2}

permdir="/data/SBGE/cory/pacbio/"
srafile=${permdir}/${run}.sra

if [ ! -f ${permdir}/${run}.sam.gz ]; then
    # if sam zip file does not yet exist

    if [ ! -f ${srafile} ]; then
        echo "downloading sra file"
        wget -O ${run}.sra https://sra-pub-run-odp.s3.amazonaws.com/sra/${run}/${run}
    else
        echo "copying sra file from /data/..."
        cp ${srafile} .
    fi


    module load sratoolkit

    echo "dumping lines from sam files..."
    sam-dump ${run}.sra  | awk '{split($12,a,":"); s=a[3]".sam"; print > s}'
    echo "done splitting sam files"

    echo "zipping sam archive"
    ls *.sam | zip -@ ${permdir}/${run}.sam.zip 
    echo "done zipping sam archive"
else
    # when sam zip file already exists
    echo "${permdir}/${run}.sam.gz already exists! "
fi

# Generate header if does not exist
if [ ! -f ${permdir}/${run}.header.txt ]; then
    sam-dump ${run}.sra | awk '{print; if (NR > 30000) exit}' | grep "^@" > ${run}.header.full.txt
    grep "^@HD" ${run}.header.full.txt > ${run}.header.txt
    grep "^@SQ" ${run}.header.full.txt >> ${run}.header.txt
    grep "^@RG" ${run}.header.full.txt | grep "ID:${cramid}"  >> ${run}.header.txt
    # currently do not include @PG lines that include
    # grep "^@PG" ${run}.header.full.txt | grep "$cramID" >> ${run}.header.txt
    cp ${run}.header.txt ${permdir}/${run}.header.txt
    # remove initial tmp header
fi

cd