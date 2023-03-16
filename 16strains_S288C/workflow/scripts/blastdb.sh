#!/bin/sh

module load blast

makeblastdb -in blast_db/ref/ref_orfs.fa -dbtype prot
makeblastdb -in blast_db/MSY24/MSY24_orfs.fa -dbtype prot
makeblastdb -in blast_db/MSY25/MSY25_orfs.fa -dbtype prot
