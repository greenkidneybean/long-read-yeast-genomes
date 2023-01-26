# Structural Variant Analysis for 16 Yeast Isolates
Note: ensure 16 reference assembly `fasta` files are in the [`assemblies/`](/assemblies/) directory.

```
273614.fasta
BY.fasta
CBS2888.fasta
CLIB219.fasta
CLIB413.fasta
I14.fasta
M22.fasta
PW5.fasta
RM.fasta
Y10.fasta
YJM145.fasta
YJM454.fasta
YJM978.fasta
YJM981.fasta
YPS1009.fasta
YPS163.fasta
```

## Retrieve reference sequence
The script [`download_reference.sh`](src/download_reference.sh) retrieves S288C reference genome
(if it isn't already downloaded) and saves it as `S288C-reference.fasta`. The script also generates
the file [`chromosome_ids.txt`](chromosome_ids.txt) for translating the NCBI codes to simpler 'chr#'
style identifiers.
```bash
src/download_reference.sh
```

## Aligning long-read assemblies to reference genome
The script [`run_MUMandCo.sh`](src/run_MUMandCo.sh) generates a genome-genome alignment for one
of the sixteen assemblies in the [`assemblies/`](/assemblies/) directory. Run all of them as a
job array with `--array=0-15` to submit 16 jobs at once. Each individual job in the array will
check the [`assemblies/`](/assemblies/) directory and then generate an alignment for that Nth
assembly. Alignments are saved to the [`reference_alignment/`](/reference_alignment/) directory.
```bash
sbatch --array=0-15 src/run_MUMandCo.sh
```

## Remove excess lines from TSVs
The coordinate aligned output from MUMandCo includes some extra (non-data) lines at the top.
Trim these out by running [fix_coords.sh](src/fix_coords.sh)
```bash
src/fix_coords.sh
```


### Generate alignment plot
```bash
module load R/3.6.3
Rscript src/plot_MUMandCo_alignment.R
```

![](/reference_alignment/all-aligned-S288C.png)