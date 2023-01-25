# Structural Variant Analysis for 16 Yeast Isolates
Note: ensure 16 reference assembly `fasta` files are in the [`assemblies/`](assemblies/) directory.
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

## Structural Variants vis-a-vis S288C Reference
### Aligning long-read assemblies to reference genome
The script [`align_reference.sh`](src/align_reference.sh) retrieves S288C reference genome (if it
isn't already downloaded) and saves it as `S288C-reference.fasta`. The script then generates `.paf`
alignments using `minimap2` for 16 genome assemblies within [`assemblies`](assemblies). Output is
saved in the [`reference_alignment`](reference_alignment) directory. 

## Generating ID table of S288C chromosome IDs
Used for translating accession IDs to chr<#>
```bash
grep '>' S288C-reference.fasta | \
sed -r 's/^>.*NC_/NC_/g' | \
sed -r 's/\| .*genomic\] \[/ /g' | \
sed -r 's/chromosome=/chr/g' | \
sed -r 's/locat.*$/mitochondrion/g' | \
sed 's/]$//g' > chromosome_ids.txt
```
---


### Generate plot with [`plot_reference_alignments.R`](src/plot_reference_alignments.R)
```bash
module load R/3.6.3
Rscript src/plot_reference_alignments.R
```

![](/06_SVs/reference_alignment/all.png)
---

## To align long-read genome assemblies to each other
The script [`align_pairwise.sh`](src/align_pairwise.sh) generates `paf` using `minimap2` for all
pairwise combinations of `fasta` within the [`assemblies`](assemblies) directory (ignoring alignment
to itself). Alignments are done reciprocally (A to B and B to A) for `N*(N-1)` total alignments, 
output to the [`pairwise_alignments`](pairwise_alignments) directory.
---


