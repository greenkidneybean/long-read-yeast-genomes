Screen genomes of 16 strains against the [1,011 strain pangenome](http://1002genomes.u-strasbg.fr/files/)

Quickstart for Biowulf HPC:
- download `orfs_{SAMPLE}.fasta` to `data/16_genomes/orf_fastas`
- Run Snakemake to generate BLASTP and BLASTN results
  `module load snakemake`
  `snakemake -j 50 --latency-wait 120 --cluster="sbatch"`
- Use Jupyter notebooks to parse the BLAST results
  ```
  conda create --name 16strains_env --file env/environment.yaml

  # create specific environment for Biowulf HPC
  # conda create --name 16strains_env --file env/spec-file.txt

  conda activate 16strains_env
  conda install ipykernel
  ipython kernel install --user --name=16strains_env
  conda activate base
  jupyter lab .
  ```

16strains_pangenome  
├── data  
│   └── 16_genomes  
│       └── orf_fastas  
│           ├── allORFs_16genomes.fasta  
│           ├── orfs_MSY24.fasta  
│           ├── orfs_MSY25.fasta  
│           ├── orfs_MSY26.fasta  
│           ├── orfs_MSY27.fasta  
│           ├── orfs_MSY28.fasta  
│           ├── orfs_MSY29.fasta  
│           ├── orfs_MSY30.fasta  
│           ├── orfs_MSY31.fasta  
│           ├── orfs_MSY32.fasta  
│           ├── orfs_MSY33.fasta  
│           ├── orfs_MSY34.fasta  
│           ├── orfs_MSY35.fasta  
│           ├── orfs_MSY36.fasta  
│           ├── orfs_MSY37.fasta  
│           ├── orfs_MSY38.fasta  
│           └── orfs_MSY39.fasta  
├── env  
│   ├── environment.yml  
│   └── spec-file.txt  
├── img  
│   ├── dag.svg  
│   └── rule_graph.svg  
├── readme.md  
└── workflow  
    ├── config.yaml  
    ├── nb  
    │   ├── 0_16genomes_check-inframe.ipynb  
    │   ├── 0_pangenome_check-inframe.ipynb  
    │   ├── 1_blastn_16genomes_top-percent-whole.ipynb  
    │   ├── 2a_blastn_16strainORFs-to-pangenome_top-percent-whole.ipynb  
    │   ├── 2b_blastn_pangenome_top-percent-whole.ipynb  
    │   ├── 2c_blastn_pangenome-pangenome.ipynb  
    │   └── 3_plot_pangenome-figure.ipynb  
    ├── scripts  
    │   └── translate_fasta.py  
    └── Snakefile  
