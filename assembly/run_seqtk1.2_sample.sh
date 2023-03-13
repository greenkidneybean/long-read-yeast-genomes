#!/bin/bash
#$ -o $JOB_ID_$JOB_NAME.stdout
#$ -e $JOB_ID_$JOB_NAME.error
#$ -l mem_free=10G
#$ -l h_vmem=10G
#$ -w e
#$ -cwd

PRG=/home/morganp/Analysis/software/seqtk/seqtk

NUM=$1

IN=`find *.fastq.gz`
OUT=`basename $IN .fastq.gz`
OUT+=".$NUM.fastq.gz"

echo "[start] seqtk"
NOW=$(date +"%Y%m%d_%H%M%S")
echo "[time] $NOW"

echo "[input] $IN"
echo "[number] $NUM"
echo "[output] $OUT"

$PRG sample -s1000 $IN $NUM | gzip > $OUT

echo "[done]"
NOW=$(date +"%Y%m%d_%H%M%S")
echo "[time] $NOW"

