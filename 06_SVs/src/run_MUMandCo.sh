#!/usr/bin/env bash
#SBATCH --ntasks-per-node=4
#SBATCH --nodes=1
#SBATCH --mem=100G
#SBATCH --time=2:00:00
#SBATCH --gres=lscratch:120
#SBATCH --partition=quick,norm


## run as --array with job number ranging from [1,16]
module load mummer/4.0.0rc1
module load samtools/1.16.1
module load blast/2.13.0+

OUTDIR="${PWD}/reference_alignment/"
mkdir -p ${OUTDIR}

# Define reference and assembly full path sequences
ref="${PWD}/S288C-reference.fasta"
assemblies=($(cd assemblies/ && ls -d $PWD/*.fasta))
N=${SLURM_ARRAY_TASK_ID}

TMPDIR=/lscratch/${SLURM_JOB_ID}
cd ${TMPDIR}
git clone https://github.com/SAMtoBAM/MUMandCo.git && cd MUMandCo

ln -s ${ref} S288C.fasta
assembly=${assemblies[$N]}
query=$(basename ${assembly%.fasta})
ln -s ${assembly} ${query}.fasta

bash mumandco_v3.8.sh \
    -r ${ref} \
    -q ${query}.fasta \
    -g 12500000 \
    -o ${query} \
    -t 4 \
    -b

mv ${query}.SVs_all.tsv ${OUTDIR}/${query}.SVs.tsv
mv ${query}_alignments/${query}_ref.delta_filter.coords ${OUTDIR}/${query}.coords

cd
exit
