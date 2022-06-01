# README

# Notes


Cross 273614 x YJM981 hosted on [NCBI](https://www.ncbi.nlm.nih.gov/sra/SRX6097474%5baccn%5d).

View run `SRR9330809` then `reads`

Reads listed like:
```
1. SRR9330809.1 SRS4998012
name: 1, member: 3003_G1_44
```

Where `3003_G1_44` describes the cross (`3003`) and segregant (`G1_44`)


## Download `sra-toolkit` 
```
wget https://ftp-trace.ncbi.nlm.nih.gov/sra/sdk/3.0.0/sratoolkit.3.0.0-ubuntu64.tar.gz

tar -zxvf sratoolkit.3.0.0-ubuntu64.tar.gz
rm sratoolkit.3.0.0-ubuntu64.tar.gz
mv sratoolkit.3.0.0-ubuntu64 sra
sra/bin/vdb-config --restore-defaults
```

## Retrieve reads
Use sra-toolkit to retrieve the desired reads in compressed format
```bash
acc='SRR9330809'
sra/bin/sam-dump  --gzip ${acc} > ${acc}.sam.gz
```

## Sorting samfile
The goal is to split a large sam file into multiple smaller `sam` files, where each one
represents a single sequenced segregant. This information is contained within the 12th
field of the sam file (`RG`, read group). The field is formatted as `RG:Z:3003_R1_20`
where `3003` denotes this specific cross, `R1` designates the plate number, and `20`
denotes the specific well in that plate (an individual segregant). To do this, I parse
the file with `awk` and print any given line to a file named for the segregant. The
following code uses a ramdisk on the biowulf HPC (because we're generating ~10,000 unique
files) to iterate line-by-line of `SRR9330809.sam.gz`:

1. parse the read group field (`$12`)
2. split the field by `:` into an array `a`
3. extract the third value in array `a`, e.g. `3003_R1_20`
4. print the entire line to a file named for its read group, e.g. `3003_R1_20.sam`
5. once all lines have been split into separate files, combine them into a single `.zip` archive.


```bash
# Split large sam file contents into separate files while on biowulf ramdisk
sinteractive --mem=200G                             # request interactive node with 200G memory
TMPDIR="/dev/shm/${USER}/${SLURM_JOB_ID}"           # name tmp directory on /dev/shm/
mkdir -p ${TMPDIR} && cd ${TMPDIR}                  # cd into tmp directory

samfile="/data/wellerca/projects/pacbio-yeast-genomes/SRR9330809.sam.gz"

zcat ${samfile} | \
awk '{split($12,a,":"); s=a[3]".sam"; print > s}'   # for each line in samfile, print to '3003_G5_96.sam', e.g.

ls 3003*.sam | zip -@ 3003.sam.zip                  # create zip archive from all newly generated files

cp 3003.sam.zip /data/wellerca/projects/            # copy off the ramdisk to permenant storage

cd && rm -rf ${TMPDIR}                              # clean up tmp directory
```

## Reconstruct header
```bash
# Extract full header
zgrep -F "@" SRR9330809.sam.gz > header-full.txt

# Include all @HD lines
grep -F "@HD" header-full.txt > header-3003.txt

# Include all @SQ lines
grep -F "@SQ" header-full.txt >> header-3003.txt

# Include all @RG lines that include 3003
grep -F "ID:3003" header-full.txt >> header-3003.txt

# Include all @PG lines that include 3003
grep -F "@PG" header-full.txt | grep -F 3003 >> header-3003.txt
```
And these header files were backed up to:

`nihbox:/cloud/pacbio-yeast-genomes/header-3003.txt`

`nihbox:/cloud/pacbio-yeast-genomes/header-full.txt`

## Regenerate fastq file and BLAST
```bash
# get blast binaries from NCBI
bash src/get-blast.sh

# retrieve header data
rclone copy nihbox:/cloud/pacbio-yeast-genomes/header-3003.txt data/
```

## Make blast database
```bash
# make blast database for telomeric sequences
./blast/bin/makeblastdb \
    -dbtype nucl \
    -in ./data/telomeric-SUC.fasta
```

## Run blast searches
```bash
# If on biowulf HPC:
sbatch --array=1-10 src/run-blast.slurm \
    3003 \
    1e-10 \
    ${PWD} \
    telomeric-SUC.fasta

# If on local machine:
for i in $(seq 1 10); do
    bash src/run-blast.slurm \
        3003 \
        1e-10 \
        ${PWD} \
        telomeric-SUC.fasta \
        $i
done &
```

## Blast vs entire chromosome
```blast
# If on local machine:
for i in $(seq 1 10); do
    bash src/run-blast.slurm \
        3003 \
        1e-10 \
        ${PWD} \
        MSY37-36-chrIV.fasta \
        $i \
        98
done &
```

## Retrieve Joshua Bloom parental genotype data
```bash
# Retrieve from cloud storage
rclone copy onedrive:/cloud/pacbio-yeast-genomes/joshua-bloom-data/  data/

# Manually download from Box:
https://www.dropbox.com/sh/jqm7a11zz9laytd/AABaE0EfQxLH6ounPhJ7yYWya?dl=0
# joshua-bloom-data/parents_w_svar_sorted.vcf.gz
# joshua-bloom-data/parents_w_svar_sorted.vcf.gz.tbi
```

## Retrieve S288c reference
```bash
rclone copy nihbox:/cloud/S288C/S288C_reference_sequence_R64-3-1_20210421.fsa.gz data/
gunzip data/S288C_reference_sequence_R64-3-1_20210421.fsa.gz
python3 src/extractContig.py data/S288C_reference_sequence_R64-3-1_20210421.fsa chromosome=IV > data/S288C-chrIV.fasta

```

## Retrieve parental pacbio contigs for chrIV
```bash
# rclone copy nihbox:/cloud/pacbio-yeast-genomes/MSY24_adjij_m3_assembly.fasta
# rclone copy nihbox:/cloud/pacbio-yeast-genomes/MSY25_adjik_m4_assembly.fasta
# rclone copy nihbox:/cloud/pacbio-yeast-genomes/MSY26_adjil_m3_assembly.fasta
# rclone copy nihbox:/cloud/pacbio-yeast-genomes/MSY27_adjim_m1_assembly.fasta
# rclone copy nihbox:/cloud/pacbio-yeast-genomes/MSY28_adjin_m2_assembly.fasta
# rclone copy nihbox:/cloud/pacbio-yeast-genomes/MSY29_adjio_m1_assembly.fasta
# rclone copy nihbox:/cloud/pacbio-yeast-genomes/MSY30_adjip_m1_assembly.fasta
# rclone copy nihbox:/cloud/pacbio-yeast-genomes/MSY31_adjiq_m1_assembly.fasta
# rclone copy nihbox:/cloud/pacbio-yeast-genomes/MSY32_adjir_m1_assembly.fasta
# rclone copy nihbox:/cloud/pacbio-yeast-genomes/MSY33_adjis_m1_assembly.fasta
# rclone copy nihbox:/clouqd/pacbio-yeast-genomes/MSY34_adjit_m1_assembly.fasta
# rclone copy nihbox:/cloud/pacbio-yeast-genomes/MSY35_adjiu_m1_assembly.fasta
rclone copy nihbox:/cloud/pacbio-yeast-genomes/MSY36_adjiv_m1_assembly.fasta data/ 
rclone copy nihbox:/cloud/pacbio-yeast-genomes/MSY37_adjiw_m2_assembly.fasta data/ 
# rclone copy nihbox:/cloud/pacbio-yeast-genomes/MSY38_adjix_m1_assembly.fasta
# rclone copy nihbox:/cloud/pacbio-yeast-genomes/MSY39_adjiy_m1_assembly.fasta
```

## Extract chrIV contigs
```bash
python3 src/extractContig.py data/MSY36_adjiv_m1_assembly.fasta chrIV > data/MSY36-chrIV.fasta
python3 src/extractContig.py data/MSY37_adjiw_m2_assembly.fasta chrIV > data/MSY37-chrIV.fasta
```

## Run minimap2

`singularity pull src/minimap.sif library://wellerca/pacbio/minimap`


```bash
singularity exec --bind $(readlink $PWD) src/minimap.sif \
    minimap2 \
    -x splice \
    -g 7000 \
    data/pacbio/pacbio_273614.fasta \
    data/pacbio/pacbio_YJM981.fasta | \
    cut -f 1-12 > minimap.out

src/mask_minimap_fastas.py data/pacbio/pacbio_YJM981.fasta data/pacbio/YJM981.bed --out data/pacbio/YJM981.mask.fasta
src/mask_minimap_fastas.py data/pacbio/pacbio_273614.fasta data/pacbio/273614.bed --out data/pacbio/273614.mask.fasta


singularity exec --bind $(readlink $PWD) src/minimap.sif \
    minimap2 \
    -x splice \
    -g 7000 \
    data/pacbio/YJM981.mask.fasta \
    data/pacbio/273614.mask.fasta | \
    cut -f 1-12


singularity exec --bind $(readlink $PWD) src/minimap.sif \
    minimap2 \
    -x asm10 \
    -g 10000 \
    data/pacbio/pacbio_273614.fasta \
    data/pacbio/pacbio_YJM981.fasta \
    > 273614-YJM981.paf

singularity exec --bind $(readlink $PWD) src/minimap.sif \
    paftools.js \
    call \
    -L 40000 \
    273614-YJM981.paf \
    > 273614-YJM981-call.txt

module load mummer
nucmer data/pacbio/273614_chrI.fa data/pacbio/YJM981_chrI.fa
show-aligns out.delta 273614_chrI YJM981_chrI
#src/minimap.sif paftools.js view -l 0 aln.sam > aln.maf
```


function splitFasta() {
    inFasta=${1}
    awk -F "(^>|\t| )" '{if($0 ~ /^>/) {s=$2".fasta";  print ">"$2 > s} else print > s}' ${inFasta}
}

export -f splitFasta

module load mummer



(cd data/pacbio; splitFasta pacbio_YJM981.fasta; splitFasta pacbio_273614.fasta)
(cd data/pacbio; nucmer 273614_chrIV.fasta YJM981_chrIV.fasta; dnadiff -d out.delta;) 

src/filter_snps.py data/pacbio/out.snps 20 > data/pacbio/273614.mask.filter.bed
src/print_target_snps.py data/pacbio/out.snps 20 > data/pacbio/273614.mask.targets.txt


# make bed file
src/mask_fasta.py data/pacbio/273614_chrIV.fasta data/pacbio/273614.filter.bed --out data/pacbio/273614.mask.fasta


# map
```bash
module load bwa
module load samtools
bwa index data/pacbio/273614.mask.fasta

src/bwa_map.sh 3003_G1_01.fastq 273614.mask
src/bwa_map.sh 3003_G1_02.fastq 273614.mask



bcftools mpileup \
    --targets-file data/pacbio/273614.targets.txt \
    --fasta-ref data/pacbio/273614.mask.fasta \
    data/bam/3003_G1_01_273614.mask.bam | \
    bcftools call \
    --ploidy 1 -m -Ob | \
    bcftools view | \
    sed 's/1:.*$/1/g' | \
    grep -v "^##" > test.vcf

samtools depth data/bam/3003_G1_01_273614.mask.bam
```

```R

library(data.table)
library(ggplot2)

sample <- '3003_G1_02'

dat <- fread(paste0('data/vcf/', sample, '.vcf'), skip="#CHROM")
setnames(dat, c('CHROM','POS',' ID','REF','ALT','QUAL','FILTER','INFO','FORMAT','GT'))
windowSize <- 2000
dat[, bin := cut(POS, breaks=seq(0,windowSize+max(dat$POS), windowSize))]
dat[, bin := as.numeric(factor(bin))]
ggplot(data=dat[, list(POS, fractionAlt = sum(GT==1)/.N), by=bin][fractionAlt > 0.9 | fractionAlt < 0.1], aes(x=POS, y=fractionAlt)) + geom_point()

library(foreach)

window_width <- 5

dat[, GT := as.numeric(GT)]

o <- foreach(i=1:(nrow(dat)-1), .combine='rbind') %do% {
        dat.sub <- dat[i:(i+window_width)]
        data.table("meanGT" = mean(dat.sub[, GT], na.rm=T),
                    "center" = median(dat.sub[,POS], na.rm=T)
                    )
    }

ggplot(o, aes(x=center, y=meanGT)) + geom_point()

dat.ag <- dat[, list(.N, fractionAlt = sum(GT==1)/.N), by=bin]
dat.ag[, POS := windowSize*bin]
ggplot(dat.ag[fractionAlt %in% c(0,1) & N > 1], aes(x=POS, y=fractionAlt)) + geom_point()

```
