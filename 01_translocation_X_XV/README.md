# Translocation between chrX and chrXV

## Call genotypes for both chromosomes
```bash
PROJDIR='/home/wellerca/pacbio-yeast-genomes/01_translocation_X_XV'
sbatch ${PROJDIR}/map_translocation.slurm CBS2888 CLIB219 chrX chrXV 
sbatch ${PROJDIR}/map_translocation.slurm CBS2888 CLIB219 chrXV chrX 
sbatch ${PROJDIR}/map_translocation.slurm CBS2888 YJM981 chrX chrXV  
sbatch ${PROJDIR}/map_translocation.slurm CBS2888 YJM981 chrXV chrX  
```

## Linkage mapping
```bash
sbatch ${PROJDIR}/linkage_map.slurm CBS2888 CLIB219 chrX chrXV 
sbatch ${PROJDIR}/linkage_map.slurm CBS2888 YJM981 chrX chrXV  
```

## Plotting
```bash
module load R/3.6.3
Rscript plot_translocation.R
```