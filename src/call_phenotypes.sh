#!/usr/bin/env bash
#SBATCH --ntasks-per-node=8
#SBATCH --nodes=1
#SBATCH --mem=48G
#SBATCH --time=3:59:59
#SBATCH --gres=lscratch:100
#SBATCH --partition=quick,norm
#SBATCH --output %j.slurm.out
#SBATCH --error %j.slurm.out


# Function Definitions
## copy_fasta ensures extension is .fa or .fasta 
copy_fasta() {
    local filename=${1}
    local destination=${2}
    case ${filename#*.} in
        'fasta.gz' | 'fa.gz' )
            echo "Supplied gzipped fasta file ${filename} but it must not be gzipped!"
            exit 1
            ;;
        'fa' | 'fasta' )
            echo 'regular fasta!'
            cp ${filename} ${destination}
            ;;
        'zip' | 'tar' | 'tar.gz' ) echo 'archives not supported, use extracted fasta file'
            ;;
        * )
            echo 'cannot determine input file type'
            ;;
    esac
}

## map_to_ref maps fastq files to defined reference and counts mapped/unmapped reads
map_to_ref() {
    local zipfile=${1}
    local fastqfile=${2}
    local id=${fastqfile%.fastq}
    echo $id
    unzip -p ${zipfile} ${fastqfile} | \
    bwa mem -B 40 -O 60 -E 10 -L 100 ref.fasta - | \
            samtools view -hb - | \
            samtools sort -o - > ${id}.bam
    samtools index ${id}.bam
    samtools idxstats ${id}.bam | awk -v s=${id} '{print s,$0}' | \
    awk 'NR==1 {s=$1; i=$4}; NR==2 {j=$5}; END {print s,i,j}' >> ${id}.txt
    rm ${id}.{bam,bam.bai}
}
export -f map_to_ref


## convenience functions for arg parsing
usage_error () { echo >&2 "$(basename $0):  $1"; exit 2; }
assert_argument () { test "$1" != "$EOL" || usage_error "$2 requires an argument"; }


# Parse Arguments
## import $@
if [ "$#" != 0 ]; then
    EOL=$(printf '\1\3\3\7')
    set -- "$@" "$EOL"
    while [ "$1" != "$EOL" ]; do
        opt="$1"; shift
        case "$opt" in

            # Your options go here.
            -h|--help) HELP='true';;
            --strain1) assert_argument "$1" "$opt"; STRAIN1="$1"; shift;;
            --strain2) assert_argument "$1" "$opt"; STRAIN2="$1"; shift;;
            --fasta) assert_argument "$1" "$opt"; FASTA="$1"; shift;;
            --chr) assert_argument "$1" "$opt"; CHR="$1"; shift;;
            --fastqs) assert_argument "$1" "$opt"; FASTQS="$1"; shift;;
            --gitdir) assert_argument "$1" "$opt"; GITDIR="$1"; shift;;
      
            # Arguments processing. You may remove any unneeded line after the 1st.
            -|''|[!-]*) set -- "$@" "$opt";;                                          # positional argument, rotate to the end
            --*=*)      set -- "${opt%%=*}" "${opt#*=}" "$@";;                        # convert '--name=arg' to '--name' 'arg'
            -[!-]?*)    set -- $(echo "${opt#-}" | sed 's/\(.\)/ -\1/g') "$@";;       # convert '-abc' to '-a' '-b' '-c'
            --)         while [ "$1" != "$EOL" ]; do set -- "$@" "$1"; shift; done;;  # process remaining arguments as positional
            -*)         usage_error "unknown option: '$opt'";;                        # catch misspelled options
            *)          usage_error "this should NEVER happen ($opt)";;               # sanity test for previous patterns
    
        esac
    done
    shift  # $EOL
fi


## Print help message if requested
if [ "${HELP}" == 'true' ]; then
cat << EndOfHelp
    --help | -h         show this message
    --strain1           first strain in cross (required)
    --strain2           second strain in cross (required)
    --fasta             diagnostic fasta for counting mapped/unmapped reads (required)
    --fastqs            zip file containing ONLY relevant sample fastqs for cross (required)
EndOfHelp
exit 0
fi


## Check for required args
argexit='false'
if [ "${STRAIN1}" == '' ]; then echo "ERROR: --strain1 is required"; argexit='true'
else echo "--strain1 is ${STRAIN1}"; fi

if [ "${STRAIN2}" == '' ]; then echo "ERROR: --strain2 is required"; argexit='true'
else echo "--strain2 is ${STRAIN2}"; fi

if [ "${FASTA}" == '' ]; then echo "ERROR: --fasta is required"; argexit='true'
else echo "--fasta is ${FASTA}"; fi

if [ "${FASTQS}" == '' ]; then echo "ERROR: --fastqs is required"; argexit='true'
else echo "--fastqs is ${FASTQS}"; fi

if [ ${argexit} == 'true' ]; then echo 'exiting due to missing arguments'; exit 1; fi


## Define SCRATCH temporary working directory
if [ -w "/lscratch/${SLURM_JOB_ID}" ]; then
    export SCRATCH="/lscratch/${SLURM_JOB_ID}"                  # use /lscratch/%j if writable
elif [ -w "/tmp" ]; then
    export SCRATCH="/tmp/${USER}/$$" && mkdir -p ${SCRATCH}     # else use /tmp if writable
else
    export SCRATCH="/home/${USER}/$$" && mkdir -p ${SCRATCH}    # else use tmp dir in home directory
fi

export FILESTEM=$(basename ${FASTA%.fasta})
## Ensure paths are not relative
export FASTA=$(realpath ${FASTA})
export FASTQS=$(realpath ${FASTQS})
export GITDIR="${GITDIR:=${PWD}}"
export GITDIR=$(realpath ${GITDIR})
echo "gitdir is ${GITDIR}"


## check for final output and exit if it exists
if [ -f ${GITDIR}/data/output/${STRAIN1}_${FILESTEM}_${STRAIN2}.counts.txt ]; then
    echo "final output already exists!"
    echo "${GITDIR}/data/output/${STRAIN1}_${FILESTEM}_${STRAIN2}.counts.txt"
    echo "Exiting"
    exit 0
fi

# RUN 


## Copy required fastas to temp working directory
echo "copying ${FASTA} to ${SCRATCH} as ref.fasta"
copy_fasta ${FASTA} ${SCRATCH}/ref.fasta
echo "copying ${FASTQS} to ${SCRATCH} as fastqs.zip"
cp ${FASTQS} ${SCRATCH}/fastqs.zip
cd ${SCRATCH}


## Load Modules
echo "Loading modules"
module load bwa/0.7.17
module load samtools/1.16.1
module load mummer/4.0.0beta2
module load python/3.9


## Index first strain
echo "Indexing ${STRAIN1} ${CHR} fasta"
bwa index ref.fasta


## Get list of fastq files within zip
files=($(unzip -l ${FASTQS} | tail -n +4 | head -n -2 | awk '{print $4}'))


## Iterate over all files within fastqs.zip
echo "Iterating over ${#files[@]} fastq files within ${FASTQS}"
parallel -j 8 map_to_ref {} ::: ${FASTQS} ::: ${files[@]}


## combine and export 
echo "Archiving output to directory ${GITDIR}/data/output"
echo "as file name ${STRAIN1}_${FILESTEM}_${STRAIN2}.counts.txt"
mkdir -p ${GITDIR}/data/output
cat *.txt > ${GITDIR}/data/output/${STRAIN1}_${FILESTEM}_${STRAIN2}.counts.txt

cd
exit 0