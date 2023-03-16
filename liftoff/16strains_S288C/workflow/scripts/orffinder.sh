#!/bin/sh
module load ORFfinder

ORFfinder -outfmt=0 -in 16_genomes/MSY24.fa -out orfs/MSY24_orfs.fa
ORFfinder -outfmt=0 -in 16_genomes/MSY25.fa -out orfs/MSY25_orfs.fa
ORFfinder -outfmt=0 -in ref/s288c.fa -out orfs/ref_orfs.fa
