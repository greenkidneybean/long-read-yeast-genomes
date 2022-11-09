# Translocation between chrX and chrXV

## Call genotypes for both chromosomes
```bash
sbatch ../src/cross_chr_genotype.sh CBS2888 CLIB219 chrX chrXV 
sbatch ../src/cross_chr_genotype.sh CBS2888 CLIB219 chrXV chrX 
sbatch ../src/cross_chr_genotype.sh CBS2888 YJM981 chrX chrXV  
sbatch ../src/cross_chr_genotype.sh CBS2888 YJM981 chrXV chrX  
```

## Linkage mapping
```bash
sbatch ../src/cross_chr_linkage_map.sh CBS2888 CLIB219 chrX chrXV 
sbatch ../src/cross_chr_linkage_map.sh CBS2888 YJM981 chrX chrXV  
```

## Plotting
```bash
module load R/3.6.3
Rscript ../src/plot_translocation.R ../data/output/CBS2888_CLIB219_chrX_chrXV_translocation.ld.txt
Rscript ../src/plot_translocation.R ../data/output/CBS2888_YJM981_chrX_chrXV_translocation.ld.txt
```