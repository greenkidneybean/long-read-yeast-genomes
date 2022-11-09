#!/usr/bin/env bash
#SBATCH --ntasks-per-node=8
#SBATCH --nodes=1
#SBATCH --mem=48G
#SBATCH --time=3:59:59
#SBATCH --gres=lscratch:100
#SBATCH --partition=quick,norm
#SBATCH --output %j.slurm.out
#SBATCH --error %j.slurm.out


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
EndOfHelp
exit 0
fi


## Check for required args
argexit='false'
if [ "${STRAIN1}" == '' ]; then echo "ERROR: --strain1 is required"; argexit='true'
else echo "--strain1 is ${STRAIN1}"; fi

if [ "${STRAIN2}" == '' ]; then echo "ERROR: --strain2 is required"; argexit='true'
else echo "--strain2 is ${STRAIN2}"; fi

if [ "${CHR}" == '' ]; then echo "ERROR: --chr is required"; argexit='true'
else echo "--chromosome is ${CHR}"; fi

if [ "${FASTA}" == '' ]; then echo "ERROR: --fasta is required"; argexit='true'
else echo "--fasta is ${FASTA}"; fi

if [ ${argexit} == 'true' ]; then echo 'exiting due to missing arguments'; exit 1; fi


# RUN 
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
export GITDIR="${GITDIR:=${PWD}}"
export GITDIR=$(realpath ${GITDIR})
echo "gitdir is ${GITDIR}"


# Exit if final output already exists
if [ -f ${GITDIR}/data/output/${STRAIN1}_${STRAIN2}_${CHR}.ld.txt ]; then
    if [ -f ${GITDIR}/data/output/${STRAIN1}_${FILESTEM}_${STRAIN2}_${CHR}.qassoc.txt ]; then
        echo "final outputs already exists!"
        echo "${GITDIR}/data/output/${STRAIN1}_${FILESTEM}_${STRAIN2}_${CHR}.qassoc.txt"
        echo "${GITDIR}/data/output/${STRAIN1}_${STRAIN2}_${CHR}.ld.txt"
        echo "Exiting"
        exit 0
    fi
fi


# RUN 
cd ${SCRATCH}



module load R/3.6.3
module load plink/1.9

# build genotype matrix
## unarchive calls 

if [ -f ${GITDIR}/data/output/${STRAIN1}_${STRAIN2}_${CHR}.call.tar.gz ]; then
    tar -zxvf ${GITDIR}/data/output/${STRAIN1}_${STRAIN2}_${CHR}.call.tar.gz
elif [ -f ${GITDIR}/data/output/${STRAIN2}_${STRAIN1}_${CHR}.call.tar.gz ]; then
    tar -zxvf ${GITDIR}/data/output/${STRAIN2}_${STRAIN1}_${CHR}.call.tar.gz
else
    echo "cannot find genotype calls archive!"
    exit 1
fi


## iterate over calls in R to build genotype matrix
Rscript ${GITDIR}/src/build_genotype_matrix_bins.R
python3 ${GITDIR}/src/write_ped.py genotypes.mat > genotypes.ped
python3 ${GITDIR}/src/write_map.py genotypes.mat 1 > genotypes.map

## Do LD mapping
plink --file genotypes --map3 --missing-genotype N
plink --bfile plink --r2 inter-chr --ld-window-r2 0

mv plink.ld ${GITDIR}/data/output/${STRAIN1}_${STRAIN2}_${CHR}.ld.txt

# build phenotype file, writes to current dir as phenotypes.txt
Rscript ${GITDIR}/src/print_phenotypes.R ${GITDIR}/data/output/${STRAIN1}_${FILESTEM}_${STRAIN2}.counts.txt

plink --file genotypes --pheno phenotypes.txt --missing-genotype N --assoc

mv plink.qassoc ${GITDIR}/data/output/${STRAIN1}_${FILESTEM}_${STRAIN2}_${CHR}.qassoc.txt

<<comment
library(data.table)
library(ggplot2)
library(viridis)
library(ggthemes)

args <- commandArgs(trailingOnly=TRUE)
strain1 <- args[1]
tfgene <- args[2]
strain2 <- args[3]
chromosome <- args[4]
datadir <- '/data/SBGE/cory/pacbio/'

assoc.filename <- paste0(datadir, strain1, '_', tfgene, '_', strain2, '_', chromosome, '.qassoc.txt')
ld.filename <- paste0(datadir, strain1, '_', strain2, '_', chromosome, '.ld.txt')

dat.assoc <- fread(assoc.filename)
dat.ld <- fread(ld.filename)
ggplot(dat.assoc, aes(x=BP, y=1, fill=-1*log10(P))) + geom_tile() + scale_fill_viridis() + theme_few() + theme(plot.background = element_rect(fill = "black"), panel.background = element_blank())

ggplot(dat.assoc, aes(x=BP, y=-1*log10(P))) + geom_point(shape=21, alpha=0.5) + theme_few()
ggplot(dat.assoc, aes(x=BP, y=-1*log10(P))) + geom_bar(stat='identity',position='dodge') + theme_few() +

ggplot(dat.ld, aes(x=BP_A, y=BP_B, fill=R2)) + geom_tile() +
theme(panel.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()
        ) +
scale_fill_viridis() +
labs(x=chromosome, y=chromosome)
comment
