#!/bin/bash
#$ -o $JOB_ID_$JOB_NAME.stdout
#$ -e $JOB_ID_$JOB_NAME.error
#$ -l mem_free=20G
#$ -l h_vmem=20G
#$ -w e
#$ -cwd

QR=$1
DB=/cluster/ifs/users/morganp/PacBio/Assembly/Sadhu_yeast/S288c_blastdb/Sc_S288c

NOW=$(date +"%Y%m%d_%H%M%S")
OUT="blastn."
OUT+="$NOW"
OUT+=".out"

blastn -db $DB -query $QR -out $OUT -evalue 1e-100 -num_descriptions 5 -num_alignments 5 


