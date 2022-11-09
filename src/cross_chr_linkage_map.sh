#!/usr/bin/env bash
#SBATCH --ntasks=3
#SBATCH --nodes=1
#SBATCH --mem=24G
#SBATCH --time=1:00:00
#SBATCH --gres=lscratch:200
#SBATCH --partition=quick,norm
#SBATCH --output %j.slurm.out
#SBATCH --error %j.slurm.out


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


module load R/3.6.3
module load plink/1.9

# build genotype matrix
## unarchive calls 

# do CHR1

tar -zxvf ${DATADIR}/${STRAIN1}_${STRAIN2}_${CHR1}_${CHR2}_translocation.call.tar.gz
# remove problematic individual without calls
[if -f SRR9330831_G1_04.call ]; then rm SRR9330831_G1_04.call; fi
Rscript ${GITDIR}/src/build_genotype_matrix_bins.R
python3 ~/pacbio-yeast-genomes/src/write_ped.py genotypes.mat > ${CHR1}.ped
python3 ~/pacbio-yeast-genomes/src/write_map.py genotypes.mat 1 > ${CHR1}.map
rm genotypes.mat
ls *.call | xargs rm

# do CHR2
tar -zxvf ${DATADIR}/${STRAIN1}_${STRAIN2}_${CHR2}_${CHR1}_translocation.call.tar.gz
# remove problematic individual without calls
[if -f SRR9330831_G1_04.call ]; then rm SRR9330831_G1_04.call; fi
Rscript ${GITDIR}/src/build_genotype_matrix_bins.R
python3 ~/pacbio-yeast-genomes/src/write_ped.py genotypes.mat > ${CHR2}.ped
python3 ~/pacbio-yeast-genomes/src/write_map.py genotypes.mat 2 > ${CHR2}.map
rm genotypes.mat
ls *.call | xargs rm

# merge
~/pacbio-yeast-genomes/mergePed.py ${CHR1}.ped ${CHR2}.ped  > merged.ped
cat ${CHR1}.map ${CHR2}.map > merged.map

# calculate LD
plink --file merged --map3 --missing-genotype N
plink --bfile plink --r2 inter-chr --ld-window-r2 0

# offload final output
mv plink.ld ${DATADIR}/${STRAIN1}_${STRAIN2}_${CHR1}_${CHR2}_translocation.ld.txt


