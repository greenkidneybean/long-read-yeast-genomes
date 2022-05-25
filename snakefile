import pandas as pd
#samples = pd.read_table("samples.tsv").set_index("sample", drop=False)
configfile: "config.yaml"
STRAINS = config['strains']


# SAMPLES = list(samples['sample'])   # SAMPLES: list of samples of given BATCH that will be processed
# BATCH = config['run']['batch']      # BATCH: identifies the group of samples we're concerned with, e.g. date of sequencing run
# GENES = list(config['amplicons'].keys())

# format wells as two digits, e.g. leading 0 if 1-9
wells = ['0' + str(x) if len(str(x))==1 else str(x) for x in range(1,97)]
plates = ['G' + str(x) for x in range(1,6)] + ['R' + str(x) for x in range(1,6)] 


# rule targets identifies the list of final output files to generate
rule targets:
    """Defines the set of files desired at end of pipeline"""
    input: "data/pacbio/blasts/.done"
    #input: expand("data/blasts/{cross}_{plate}_{well}.blast.gz", cross='3003', plate=plates, well=wells)
    #input: expand("data/imputed/{cross}_{plate}_{well}.txt", cross="3003", plate=plates, well=wells)


rule make_blast_db:
    input: 
        expand("data/pacbio/pacbio_{strain}.fasta", strain = STRAINS)
    output: 
        "data/pacbio/pacbio_db.fasta"
    threads: 1
    resources:
        mem_mb = 1024*2,
        runtime_min = 5,
    priority: 1
    shell:
        """
        cat {input}> {output}
        ./blast/bin/makeblastdb \
            -dbtype nucl \
            -in {output}
        """


# rule convert_sam_to_fasta:
#     input: "data/3003.sam.zip"
#     output: temp("data/fasta/{cross}_{plate}_{well}.fasta")
#     threads: 1
#     resources:
#         mem_mb = 1024*2,
#         runtime_min = 20,
#     container: "library://wellerca/pseudodiploidy/mapping:latest"
#     priority: 2
#     shell:
#         """
#         samtools fasta \
#         -1 data/fasta/{wildcards.cross}_{wildcards.plate}_{wildcards.well}.fasta \
#         -2 data/fasta/{wildcards.cross}_{wildcards.plate}_{wildcards.well}.fasta \
#         -0 data/fasta/{wildcards.cross}_{wildcards.plate}_{wildcards.well}.fasta \
#         -s data/fasta/{wildcards.cross}_{wildcards.plate}_{wildcards.well}.fasta \
#         <(cat data/header-{wildcards.cross}.txt <(unzip -p  {input} {wildcards.cross}_{wildcards.plate}_{wildcards.well}.sam))
#         """

rule run_blastn:
    input: blast_db = "data/pacbio/pacbio_db.fasta"
    output: touch("data/blasts/.done")
    threads: 8
    resources:
        mem_mb = 1024*8,
        runtime_min = 120,
    params:
        pct_identity = lambda wildcards: config["pct_identity"],
        e_value_cutoff = lambda wildcards: config["e_value_cutoff"],
        cross = lambda wildcards: config["cross"]
    priority: 3
    container: "library://wellerca/pseudodiploidy/mapping:latest"
    shell:
        """
        src/convert_sam_to_fasta.sh  data/{params.cross}.sam.zip {params.pct_identity} {threads} {params.e_value_cutoff} {input.blast_db}
        """

rule impute_haplotypes:
    input: "data/blasts/.done"
    output: touch("data/imputed/.done")
    threads: 1
    resources:
        mem_mb = 1024*2,
        runtime_min = 5,
    params:
        parent1 = STRAINS[0]
        parent2 = STRAINS[1]
    priority: 4
    shell:
        """
        module load R/3.6.3
        Rscript src/linkage_mapping.R {params.parent1} {params.parent2}
        """
