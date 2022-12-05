#!/usr/bin/env bash

source /data/SBGE/conda/etc/profile.d/conda.sh

conda create --name assembly biopython
conda config --add channels defaults
conda config --add channels bioconda
conda config --add channels conda-forge

conda activate assembly

conda install liftoff

liftoff pacbio_genomes/MSY24_adjij_m3_assembly.fasta pacbio_genomes/reference/SacCer3.fna -g pacbio_genomes/reference/sacCer3_sgd_short.gff -o pacbio_genomes/output/MSY24_liftoff2.gff -u pacbio_genomes/output/MSY24_unmapped_features2.txt -p 4 -infer_transcripts -copies -s .3 -sc .3

liftoff pacbio_genomes/MSY25_adjik_m3_assembly.fasta pacbio_genomes/reference/SacCer3.fna -g pacbio_genomes/reference/sacCer3_sgd_short.gff -o pacbio_genomes/output/MSY25_liftoff2.gff -u pacbio_genomes/output/MSY25_unmapped_features2.txt -p 4 -infer_transcripts -copies -s .3 -sc .3

liftoff pacbio_genomes/MSY26_adjil_m3_assembly.fasta pacbio_genomes/reference/SacCer3.fna -g pacbio_genomes/reference/sacCer3_sgd_short.gff -o pacbio_genomes/output/MSY26_liftoff2.gff -u pacbio_genomes/output/MSY26_unmapped_features2.txt -p 4 -infer_transcripts -copies -s .3 -sc .3

liftoff pacbio_genomes/MSY27_adjim_m3_assembly.fasta pacbio_genomes/reference/SacCer3.fna -g pacbio_genomes/reference/sacCer3_sgd_short.gff -o pacbio_genomes/output/MSY27_liftoff2.gff -u pacbio_genomes/output/MSY27_unmapped_features2.txt -p 4 -infer_transcripts -copies -s .3 -sc .3

liftoff pacbio_genomes/MSY28_adjin_m3_assembly.fasta pacbio_genomes/reference/SacCer3.fna -g pacbio_genomes/reference/sacCer3_sgd_short.gff -o pacbio_genomes/output/MSY28_liftoff2.gff -u pacbio_genomes/output/MSY28_unmapped_features2.txt -p 4 -infer_transcripts -copies -s .3 -sc .3

liftoff pacbio_genomes/MSY29_adjio_m3_assembly.fasta pacbio_genomes/reference/SacCer3.fna -g pacbio_genomes/reference/sacCer3_sgd_short.gff -o pacbio_genomes/output/MSY29_liftoff2.gff -u pacbio_genomes/output/MSY29_unmapped_features2.txt -p 4 -infer_transcripts -copies -s .3 -sc .3

liftoff pacbio_genomes/MSY30_adjip_m3_assembly.fasta pacbio_genomes/reference/SacCer3.fna -g pacbio_genomes/reference/sacCer3_sgd_short.gff -o pacbio_genomes/output/MSY30_liftoff2.gff -u pacbio_genomes/output/MSY30_unmapped_features2.txt -p 4 -infer_transcripts -copies -s .3 -sc .3

liftoff pacbio_genomes/MSY31_adjiq_m3_assembly.fasta pacbio_genomes/reference/SacCer3.fna -g pacbio_genomes/reference/sacCer3_sgd_short.gff -o pacbio_genomes/output/MSY31_liftoff2.gff -u pacbio_genomes/output/MSY31_unmapped_features2.txt -p 4 -infer_transcripts -copies -s .3 -sc .3

liftoff pacbio_genomes/MSY32_adjir_m3_assembly.fasta pacbio_genomes/reference/SacCer3.fna -g pacbio_genomes/reference/sacCer3_sgd_short.gff -o pacbio_genomes/output/MSY32_liftoff2.gff -u pacbio_genomes/output/MSY32_unmapped_features2.txt -p 4 -infer_transcripts -copies -s .3 -sc .3

liftoff pacbio_genomes/MSY33_adjis_m3_assembly.fasta pacbio_genomes/reference/SacCer3.fna -g pacbio_genomes/reference/sacCer3_sgd_short.gff -o pacbio_genomes/output/MSY33_liftoff2.gff -u pacbio_genomes/output/MSY33_unmapped_features2.txt -p 4 -infer_transcripts -copies -s .3 -sc .3

liftoff pacbio_genomes/MSY34_adjit_m3_assembly.fasta pacbio_genomes/reference/SacCer3.fna -g pacbio_genomes/reference/sacCer3_sgd_short.gff -o pacbio_genomes/output/MSY34_liftoff2.gff -u pacbio_genomes/output/MSY34_unmapped_features2.txt -p 4 -infer_transcripts -copies -s .3 -sc .3

liftoff pacbio_genomes/MSY35_adjiu_m3_assembly.fasta pacbio_genomes/reference/SacCer3.fna -g pacbio_genomes/reference/sacCer3_sgd_short.gff -o pacbio_genomes/output/MSY35_liftoff2.gff -u pacbio_genomes/output/MSY35_unmapped_features2.txt -p 4 -infer_transcripts -copies -s .3 -sc .3

liftoff pacbio_genomes/MSY36_adjiv_m3_assembly.fasta pacbio_genomes/reference/SacCer3.fna -g pacbio_genomes/reference/sacCer3_sgd_short.gff -o pacbio_genomes/output/MSY36_liftoff2.gff -u pacbio_genomes/output/MSY36_unmapped_features2.txt -p 4 -infer_transcripts -copies -s .3 -sc .3

liftoff pacbio_genomes/MSY37_adjiw_m3_assembly.fasta pacbio_genomes/reference/SacCer3.fna -g pacbio_genomes/reference/sacCer3_sgd_short.gff -o pacbio_genomes/output/MSY37_liftoff2.gff -u pacbio_genomes/output/MSY37_unmapped_features2.txt -p 4 -infer_transcripts -copies -s .3 -sc .3

liftoff pacbio_genomes/MSY38_adjix_m3_assembly.fasta pacbio_genomes/reference/SacCer3.fna -g pacbio_genomes/reference/sacCer3_sgd_short.gff -o pacbio_genomes/output/MSY38_liftoff2.gff -u pacbio_genomes/output/MSY38_unmapped_features2.txt -p 4 -infer_transcripts -copies -s .3 -sc .3

liftoff pacbio_genomes/MSY39_adjiy_m3_assembly.fasta pacbio_genomes/reference/SacCer3.fna -g pacbio_genomes/reference/sacCer3_sgd_short.gff -o pacbio_genomes/output/MSY39_liftoff2.gff -u pacbio_genomes/output/MSY39_unmapped_features2.txt -p 4 -infer_transcripts -copies -s .3 -sc .3


liftoff pacbio_genomes/MSY24_adjij_m3_assembly.fasta pacbio_genomes/reference/SacCer3.fna -g pacbio_genomes/reference/sacCer3_sgd_short_extra_features_as_genes.gff -o pacbio_genomes/gff_files/221019/MSY24extra.gff -u pacbio_genomes/gff_files/unmapped_feature_files/221019msy24tns.txt -p 4 -infer_transcripts -copies -s .3 -sc .3

liftoff pacbio_genomes/MSY25_adjik_m3_assembly.fasta pacbio_genomes/reference/SacCer3.fna -g pacbio_genomes/reference/sacCer3_sgd_short_extra_features_as_genes.gff -o pacbio_genomes/gff_files/221019/MSY25extra.gff -u pacbio_genomes/gff_files/unmapped_feature_files/221019msy25tns.txt -p 4 -infer_transcripts -copies -s .3 -sc .3

liftoff pacbio_genomes/MSY26_adjil_m3_assembly.fasta pacbio_genomes/reference/SacCer3.fna -g pacbio_genomes/reference/sacCer3_sgd_short_extra_features_as_genes.gff -o pacbio_genomes/gff_files/221019/MSY26extra.gff -u pacbio_genomes/gff_files/unmapped_feature_files/221019msy26tns.txt -p 4 -infer_transcripts -copies -s .3 -sc .3

liftoff pacbio_genomes/MSY27_adjim_m3_assembly.fasta pacbio_genomes/reference/SacCer3.fna -g pacbio_genomes/reference/sacCer3_sgd_short_extra_features_as_genes.gff -o pacbio_genomes/gff_files/221019/MSY27extra.gff -u pacbio_genomes/gff_files/unmapped_feature_files/221019msy27tns.txt -p 4 -infer_transcripts -copies -s .3 -sc .3

liftoff pacbio_genomes/MSY28_adjin_m3_assembly.fasta pacbio_genomes/reference/SacCer3.fna -g pacbio_genomes/reference/sacCer3_sgd_short_extra_features_as_genes.gff -o pacbio_genomes/gff_files/221019/MSY28extra.gff -u pacbio_genomes/gff_files/unmapped_feature_files/221019msy28tns.txt -p 4 -infer_transcripts -copies -s .3 -sc .3

liftoff pacbio_genomes/MSY29_adjio_m3_assembly.fasta pacbio_genomes/reference/SacCer3.fna -g pacbio_genomes/reference/sacCer3_sgd_short_extra_features_as_genes.gff -o pacbio_genomes/gff_files/221019/MSY29extra.gff -u pacbio_genomes/gff_files/unmapped_feature_files/221019msy29tns.txt -p 4 -infer_transcripts -copies -s .3 -sc .3

liftoff pacbio_genomes/MSY30_adjip_m3_assembly.fasta pacbio_genomes/reference/SacCer3.fna -g pacbio_genomes/reference/sacCer3_sgd_short_extra_features_as_genes.gff -o pacbio_genomes/gff_files/221019/MSY30extra.gff -u pacbio_genomes/gff_files/unmapped_feature_files/221019msy30tns.txt -p 4 -infer_transcripts -copies -s .3 -sc .3

liftoff pacbio_genomes/MSY31_adjiq_m3_assembly.fasta pacbio_genomes/reference/SacCer3.fna -g pacbio_genomes/reference/sacCer3_sgd_short_extra_features_as_genes.gff -o pacbio_genomes/gff_files/221019/MSY31extra.gff -u pacbio_genomes/gff_files/unmapped_feature_files/221019msy31tns.txt -p 4 -infer_transcripts -copies -s .3 -sc .3

liftoff pacbio_genomes/MSY32_adjir_m3_assembly.fasta pacbio_genomes/reference/SacCer3.fna -g pacbio_genomes/reference/sacCer3_sgd_short_extra_features_as_genes.gff -o pacbio_genomes/gff_files/221019/MSY32extra.gff -u pacbio_genomes/gff_files/unmapped_feature_files/221019msy32tns.txt -p 4 -infer_transcripts -copies -s .3 -sc .3

liftoff pacbio_genomes/MSY33_adjis_m3_assembly.fasta pacbio_genomes/reference/SacCer3.fna -g pacbio_genomes/reference/sacCer3_sgd_short_extra_features_as_genes.gff -o pacbio_genomes/gff_files/221019/MSY33extra.gff -u pacbio_genomes/gff_files/unmapped_feature_files/221019msy33tns.txt -p 4 -infer_transcripts -copies -s .3 -sc .3

liftoff pacbio_genomes/MSY34_adjit_m3_assembly.fasta pacbio_genomes/reference/SacCer3.fna -g pacbio_genomes/reference/sacCer3_sgd_short_extra_features_as_genes.gff -o pacbio_genomes/gff_files/221019/MSY34extra.gff -u pacbio_genomes/gff_files/unmapped_feature_files/221019msy34tns.txt -p 4 -infer_transcripts -copies -s .3 -sc .3

liftoff pacbio_genomes/MSY35_adjiu_m3_assembly.fasta pacbio_genomes/reference/SacCer3.fna -g pacbio_genomes/reference/sacCer3_sgd_short_extra_features_as_genes.gff -o pacbio_genomes/gff_files/221019/MSY35extra.gff -u pacbio_genomes/gff_files/unmapped_feature_files/221019msy35tns.txt -p 4 -infer_transcripts -copies -s .3 -sc .3

liftoff pacbio_genomes/MSY36_adjiv_m3_assembly.fasta pacbio_genomes/reference/SacCer3.fna -g pacbio_genomes/reference/sacCer3_sgd_short_extra_features_as_genes.gff -o pacbio_genomes/gff_files/221019/MSY36extra.gff -u pacbio_genomes/gff_files/unmapped_feature_files/221019msy36tns.txt -p 4 -infer_transcripts -copies -s .3 -sc .3

liftoff pacbio_genomes/MSY37_adjiw_m3_assembly.fasta pacbio_genomes/reference/SacCer3.fna -g pacbio_genomes/reference/sacCer3_sgd_short_extra_features_as_genes.gff -o pacbio_genomes/gff_files/221019/MSY37extra.gff -u pacbio_genomes/gff_files/unmapped_feature_files/221019msy37tns.txt -p 4 -infer_transcripts -copies -s .3 -sc .3

liftoff pacbio_genomes/MSY38_adjix_m3_assembly.fasta pacbio_genomes/reference/SacCer3.fna -g pacbio_genomes/reference/sacCer3_sgd_short_extra_features_as_genes.gff -o pacbio_genomes/gff_files/221019/MSY38extra.gff -u pacbio_genomes/gff_files/unmapped_feature_files/221019msy38tns.txt -p 4 -infer_transcripts -copies -s .3 -sc .3

liftoff pacbio_genomes/MSY39_adjiy_m3_assembly.fasta pacbio_genomes/reference/SacCer3.fna -g pacbio_genomes/reference/sacCer3_sgd_short_extra_features_as_genes.gff -o pacbio_genomes/gff_files/221019/MSY39extra.gff -u pacbio_genomes/gff_files/unmapped_feature_files/221019msy39tns.txt -p 4 -infer_transcripts -copies -s .3 -sc .3
