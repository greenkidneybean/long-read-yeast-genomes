Screen genomes of 16 strains against the [S288C reference genome](http://sgd-archive.yeastgenome.org/sequence/S288C_reference)

Quickstart for Biowulf HPC:
- download `{SAMPLE}.fasta` to `data/16_genomes/`
- Run Snakemake to generate ORFs and BLASTP results
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

16strains_S288C  
├── data  
│   └── 16_genomes  
│       ├── MSY24.fasta  
│       ├── MSY25.fasta  
│       ├── MSY26.fasta  
│       ├── MSY27.fasta  
│       ├── MSY28.fasta  
│       ├── MSY29.fasta  
│       ├── MSY30.fasta  
│       ├── MSY31.fasta  
│       ├── MSY32.fasta  
│       ├── MSY33.fasta  
│       ├── MSY34.fasta  
│       ├── MSY35.fasta  
│       ├── MSY36.fasta  
│       ├── MSY37.fasta  
│       ├── MSY38.fasta  
│       └── MSY39.fasta  
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
    │   ├── 2a_ref-blast-hits-table.ipynb  
    │   ├── 2b_s288c-orf_blast-hits-table.ipynb  
    │   ├── 3a_unique-from-ref.ipynb  
    │   ├── 4_blast-gene-hits-concat.ipynb  
    │   └── extract_unique_orfs.ipynb  
    ├── scripts  
    │   ├── blastdb.sh  
    │   ├── blast.sh  
    │   ├── clean_fastas.sh  
    │   ├── genomes.txt  
    │   ├── orffinder.sh  
    │   └── unique_orfs.py  
    └── snakefile  
