#!/bin/sh

module load blast

blastp -query orfs/ref_orfs.fa -db blast_db/ref/ref_orfs.fa -outfmt 6 -max_hsps=1 -out blast_out/ref/ref-ref_blast.tsv
blastp -query orfs/ref_orfs.fa -db blast_db/MSY24/MSY24_orfs.fa -outfmt 6 -max_hsps=1 -out blast_out/ref/ref-MSY24_blast.tsv
blastp -query orfs/ref_orfs.fa -db blast_db/MSY25/MSY25_orfs.fa -outfmt 6 -max_hsps=1 -out blast_out/ref/ref-MSY25_blast.tsv

blastp -query orfs/MSY24_orfs.fa -db blast_db/ref/ref_orfs.fa -outfmt 6 -max_hsps=1 -out blast_out/MSY24/MSY24-ref_blast.tsv
blastp -query orfs/MSY24_orfs.fa -db blast_db/MSY24/MSY24_orfs.fa -outfmt 6 -max_hsps=1 -out blast_out/MSY24/MSY24-MSY24_blast.tsv
blastp -query orfs/MSY24_orfs.fa -db blast_db/MSY25/MSY25_orfs.fa -outfmt 6 -max_hsps=1 -out blast_out/MSY24/MSY24-MSY25_blast.tsv

blastp -query orfs/MSY25_orfs.fa -db blast_db/ref/ref_orfs.fa -outfmt 6 -max_hsps=1 -out blast_out/MSY25/MSY25-ref_blast.tsv
blastp -query orfs/MSY25_orfs.fa -db blast_db/MSY24/MSY24_orfs.fa -outfmt 6 -max_hsps=1 -out blast_out/MSY25/MSY25-MSY24_blast.tsv
blastp -query orfs/MSY25_orfs.fa -db blast_db/MSY25/MSY25_orfs.fa -outfmt 6 -max_hsps=1 -out blast_out/MSY25/MSY25-MSY25_blast.tsv
