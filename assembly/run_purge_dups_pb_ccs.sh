#!/bin/bash
#$ -o $JOB_ID_$JOB_NAME.stdout
#$ -e $JOB_ID_$JOB_NAME.error
#$ -l mem_free=10G
#$ -l h_vmem=10G
#$ -w e
#$ -cwd

source ~/.bash_profile
export PATH=/home/morganp/Analysis/software/minimap2-2.17_x64-linux/:/home/morganp/Analysis/software/purge_dups-v1.0.1/bin:$PATH

FQZ=`find *.fastq.gz`
if [ -f $FQZ ]
then
  echo "[input fastq.gz] $FQZ"
else 
  echo "[error] missing input fastq.gz"
  exit
fi

ASM=`find *.fasta`
if [ -f $ASM ]
then
  echo "[assembly] $ASM"
else
  echo "[error] missing assembly fasta"
  exit
fi

module load compiler/gcc/4.8.5
module load libs/zlib/1.2.8

echo "[start] step 0"
NOW=$(date +"%Y%m%d_%H%M%S")
echo "[time] $NOW"

echo "[cmd] minimap2 -x asm20 $ASM $FQZ | gzip -c - > map.paf.gz"
minimap2 -x asm20 $ASM $FQZ | gzip -c - > map.paf.gz

NOW=$(date +"%Y%m%d_%H%M%S")
echo "[time] $NOW"
echo "[cmd] pbcstat map.paf.gz"
pbcstat map.paf.gz

NOW=$(date +"%Y%m%d_%H%M%S")
echo "[time] $NOW"
echo "[cmd] calcuts PB.stat > cutoffs 2>calcuts.log"
calcuts PB.stat > cutoffs 2>calcuts.log

echo "[start] step 1"
NOW=$(date +"%Y%m%d_%H%M%S")
echo "[time] $NOW"

echo "[cmd] split_fa $ASM > $ASM.split"
split_fa $ASM > $ASM.split

NOW=$(date +"%Y%m%d_%H%M%S")
echo "[time] $NOW"
echo "[cmd] minimap2 -x asm5 -DP $ASM.split $ASM.split | gzip -c - > $ASM.split.self.paf.gz"
minimap2 -x asm5 -DP $ASM.split $ASM.split | gzip -c - > $ASM.split.self.paf.gz

echo "[start] step 2"
NOW=$(date +"%Y%m%d_%H%M%S")
echo "[time] $NOW"

echo "[cmd] purge_dups -2 -T cutoffs -c PB.base.cov $ASM.split.self.paf.gz > dups.bed 2> purge_dups.log"
purge_dups -2 -T cutoffs -c PB.base.cov $ASM.split.self.paf.gz > dups.bed 2> purge_dups.log

echo "[start] step 3"
NOW=$(date +"%Y%m%d_%H%M%S")
echo "[time] $NOW"

echo "[cmd] get_seqs -e dups.bed $ASM"
get_seqs -e dups.bed $ASM

echo "[stop]"
NOW=$(date +"%Y%m%d_%H%M%S")
echo "[time] $NOW"


module unload compiler/gcc/4.8.5
module unload /libs/zlib/1.2.8

