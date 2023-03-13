#!/bin/bash
#$ -o $JOB_ID_$JOB_NAME.stdout
#$ -e $JOB_ID_$JOB_NAME.error
#$ -pe make-dedicated 16
#$ -l mem_free=4G
#$ -l h_vmem=4G
#$ -w e
#$ -cwd

source ~/.bash_profile 

#sample, source name
PT=$1
DR=$PT
SZ=12000000
FQ=`find *.fastq.gz`

echo "[start canu-2.2 assembly]"
NOW=$(date +"%Y%m%d_%H%M%S")
echo "[time] $NOW"

#add commands

module load all/Java/1.8.0_60
module load base/perl/5.20.2
module load GCC/8.2.0-2.31.1

export PERL5LIB=/home/morganp/HMPgenomes/pipeline/software/canu-2.2/lib/site_perl/canu:$PERL5LIB

#run with default coverage
#change SGE options param name
/home/morganp/HMPgenomes/pipeline/software/canu-2.2/bin/canu -p $PT -d $DR genomeSize=$SZ java="/opt/sw/software/Java/1.8.0_60/bin/java" gridEngineResourceOption="-pe make-dedicated THREADS -l mem_free=MEMORY -l h_vmem=MEMORY" gridOptions="-V -S /bin/bash" -pacbio-hifi $FQ

module unload all/Java/1.8.0_60
module unload base/perl/5.20.2 
module unload GCC/8.2.0-2.31.1

echo "[finish]"
NOW=$(date +"%Y%m%d_%H%M%S")
echo "[time] $NOW"
