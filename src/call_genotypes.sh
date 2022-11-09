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

## map_and_call streams fastq file from zip archive and pipes through
## bwa mem and bcftools to call genotypes
map_and_call() {
    local zipfile=${1}
    local fastqfile=${2}
    local id=${fastqfile%.fastq}

    bwa mem mask.fasta <((unzip -p ${zipfile} ${fastqfile})) | \
        samtools view -hb - | \
        samtools sort - > ${id}.bam

    samtools index ${id}.bam
    bcftools mpileup \
        --targets-file targets.txt \
        --fasta-ref mask.fasta \
        ${id}.bam | \
        bcftools call \
        --ploidy 1 -m -Ob | \
        bcftools view | \
        sed 's/1:.*$/1/g' | \
        grep -v "^##" | awk '{print $1,$2,$4,$5,$10}' | sed 's/\.bam$//g' > ${id}.call
    
    rm ${id}.{bam,bam.bai}
}
export -f map_and_call


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
            --fasta1) assert_argument "$1" "$opt"; FASTA1="$1"; shift;;
            --strain2) assert_argument "$1" "$opt"; STRAIN2="$1"; shift;;
            --fasta2) assert_argument "$1" "$opt"; FASTA2="$1"; shift;;
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
    --fasta1            fasta file for strain1(required)
    --strain2           second strain in cross (required)
    --fasta2            fasta file for strain2 (required)
    --chr               chromosome being mapped (required)
    --fastqs            zip file containing ONLY relevant sample fastqs for cross (required)
EndOfHelp
exit 0
fi


## Check for required args
argexit='false'
if [ "${CHR}" == '' ]; then echo "ERROR: --chr is required"; argexit='true'
else echo "--chromosome is ${CHR}"; fi

if [ "${STRAIN1}" == '' ]; then echo "ERROR: --strain1 is required"; argexit='true'
else echo "--strain1 is ${STRAIN1}"; fi

if [ "${FASTA1}" == '' ]; then echo "ERROR: --fasta1 is required"; argexit='true'
else echo "--fasta1 is ${FASTA1}"; fi

if [ "${STRAIN2}" == '' ]; then echo "ERROR: --strain2 is required"; argexit='true'
else echo "--strain2 is ${STRAIN2}"; fi

if [ "${FASTA2}" == '' ]; then echo "ERROR: --fasta2 is required"; argexit='true'
else echo "--fasta2 is ${FASTA2}"; fi

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

# Exit if final output already exists
if [ -f ${GITDIR}/data/output/${STRAIN1}_${STRAIN2}_${CHR}.call.tar.gz ]; then
    echo "final output already exists!"
    echo "${GITDIR}/data/output/${STRAIN1}_${STRAIN2}_${CHR}.call.tar.gz"
    echo "Exiting"
    exit 0
fi
# RUN 



## Ensure paths are not relative
export FASTA1=$(realpath ${FASTA1})
export FASTA2=$(realpath ${FASTA2})
export FASTQS=$(realpath ${FASTQS})
export GITDIR="${GITDIR:=${PWD}}"
export GITDIR=$(realpath ${GITDIR})
echo "gitdir is ${GITDIR}"


## Copy required fastas to temp working directory
echo "copying ${FASTA1} to ${SCRATCH} as strain1.fasta"
copy_fasta ${FASTA1} ${SCRATCH}/strain1.fasta
echo "copying ${FASTA2} to ${SCRATCH} as strain2.fasta"
copy_fasta ${FASTA2} ${SCRATCH}/strain2.fasta
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
bwa index strain1.fasta



## Perform chr-to-chr alignment
echo "Aligning chromosome-to-chromosome"
nucmer strain1.fasta strain2.fasta && \
dnadiff -d out.delta


## Ensure alignment map is only from perspective of first strain
awk -v s1=${STRAIN1} '$11~s1' out.snps > out.fixed.snps
rm out.snps
mv out.fixed.snps out.snps


## Prepare masked reference fasta
echo "Filtering snp output from mummer"
python3 ${GITDIR}/src/filter_snps.py out.snps 20 > mask.filter.bed
echo "Getting SNP targets for bcftools"
python3 ${GITDIR}/src/print_target_snps.py out.snps 20 > targets.txt
echo "Generating masked reference fasta for bwa"
python3 ${GITDIR}/src/mask_fasta.py \
    strain1.fasta \
    mask.filter.bed \
    --out mask.fasta
echo "Indexing masked fasta"
bwa index mask.fasta


## Get list of files within zip
files=($(unzip -l ${FASTQS} | tail -n +4 | head -n -2 | awk '{print $4}'))

## Iterate over all files within fastqs.zip
echo "Iterating over ${#files[@]} fastq files within ${FASTQS}"
parallel -j 8 map_and_call {} ::: ${FASTQS} ::: ${files[@]}


## combine and export 
echo "Archiving output to directory ${GITDIR}/data"
echo "as file name ${STRAIN1}_${STRAIN2}_${CHR}.call.tar.gz"
mkdir -p ${GITDIR}/data/output
ls *.call | tar -czv --files-from - -f ${GITDIR}/data/output/${STRAIN1}_${STRAIN2}_${CHR}.call.tar.gz

cd
exit 0