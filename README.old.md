# README

## Retrieve Assemblies

# create and set up / link data folder
/path/to/repo/data
            |---assemblies
            |---shortreads
            |---output
            |---samfiles

# Notes

Bioproject hosted on NCBI [here](https://www.ncbi.nlm.nih.gov/sra?linkname=bioproject_sra_all&from_uid=549760)

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
acc=SRR9330809
sra/bin/sam-dump  --gzip ${acc} > ${acc}.sam.gz
```

## Sorting samfile
The goal is to split a large sam file into multiple smaller `sam` files, where each one
represents a single sequenced segregant. This information is contained within the 12th
field of the sam file (`RG`, read group). The field is formatted as `RG:Z:3003_R1_20`
where `3003` denotes this specific cross, `R1` designates the plate number, and `20`
denotes the specific well in that plate (an individual segregant). To do this, I parse
the file with `awk` and print any given line to a file named for the segregant. The
following code uses a ramdisk on the biowulf HPC (because were generating ~10,000 unique
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
awk {split($12,a,":"); s=a[3]".sam"; print > s}   # for each line in samfile, print to 3003_G5_96.sam, e.g.

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
    awk -F "(^>|\t| )" {if($0 ~ /^>/) {s=$2".fasta";  print ">"$2 > s} else print > s} ${inFasta}
}

export -f splitFasta

module load mummer


# split genome into independent chromosome fasta files
(cd data/pacbio; splitFasta pacbio_YJM981.fasta; splitFasta pacbio_273614.fasta)

# create 1-to-1 genome mapping
(cd data/pacbio; nucmer 273614_chrIV.fasta YJM981_chrIV.fasta; dnadiff -d out.delta;) 

# filter SNPs to good regions
src/filter_snps.py data/pacbio/out.snps 20 > data/pacbio/273614.mask.filter.bed

# make target SNP file for variant calling (building VCF)
src/print_target_snps.py data/pacbio/out.snps 20 > data/pacbio/273614.mask.targets.txt


# make bed file
src/mask_fasta.py data/pacbio/273614_chrIV.fasta data/pacbio/273614.mask.filter.bed --out data/pacbio/273614.mask.fasta


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
    sed s/1:.*$/1/g | \
    grep -v "^##" > test.vcf

samtools depth data/bam/3003_G1_01_273614.mask.bam
```

```R

library(data.table)
library(ggplot2)

sample <- 3003_G1_02

dat <- fread(paste0(data/vcf/, sample, .vcf), skip="#CHROM")
setnames(dat, c(CHROM,POS, ID,REF,ALT,QUAL,FILTER,INFO,FORMAT,GT))
windowSize <- 2000
dat[, bin := cut(POS, breaks=seq(0,windowSize+max(dat$POS), windowSize))]
dat[, bin := as.numeric(factor(bin))]
ggplot(data=dat[, list(POS, fractionAlt = sum(GT==1)/.N), by=bin][fractionAlt > 0.95 | fractionAlt < 0.05], aes(x=POS, y=fractionAlt)) + geom_point()

library(foreach)

window_width <- 5

dat[, GT := as.numeric(GT)]

o <- foreach(i=1:(nrow(dat)-1), .combine=rbind) %do% {
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


# Make genotypes
using `src/calculateLD.R` to write `data/vcf/chrIV/genotypes.mat` then `src/write_ped.py` to convert
to ped format for `PLINK`

```
python src/write_ped.py data/vcf/chrIV/genotypes.mat > data/vcf/chrIV/genotypes.ped
python src/write_map.py data/vcf/chrIV/genotypes.mat > data/vcf/chrIV/genotypes.map

plink --file data/vcf/chrIV/genotypes --map3 --missing-genotype N
plink --bfile plink --r2 inter-chr --ld-window-r2 0

library(data.table)
library(ggplot2)
library(viridis)
dat <- fread(plink.ld)
ggplot(dat, aes(x=BP_A, y=BP_B, fill=R2)) + geom_tile() + theme(plot.background = element_rect(fill = "black")) +
theme(panel.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()
        ) +
scale_fill_viridis()
```

### UPDATED 10-19-2022

function get_samfile_list() {
    local RUN=${1}
    local samfiles=$(zip -sf ${RUN}.sam.zip | tail -n +2 | head -n -1)
    echo ${samfiles}
}

# Print only reads from fastq file
zcat $file | tail -n +2 | awk NR%4==1 > onlyreads.txt

function build_fastq() {
    local RUN=${1}
    local ZIP_FILENAME=${RUN}.sam.zip
    local SAM_FILENAME=${2}
    local FILESTEM=$(echo ${SAM_FILENAME%.sam})
    local CROSS=$(echo $FILESTEM | cut -d _ -f 1)
    local PLATE=$(echo $FILESTEM | cut -d _ -f 2)
    local WELL=$(echo $FILESTEM | cut -d _ -f 3)
    samtools fastq \
    -1 ${RUN}_${PLATE}_${WELL}.fastq \
    -2 ${RUN}_${PLATE}_${WELL}.fastq \
    -0 ${RUN}_${PLATE}_${WELL}.fastq \
    -s ${RUN}_${PLATE}_${WELL}.fastq \
    <(cat ${RUN}.header.txt <(unzip -p  ${ZIP_FILENAME} ${CROSS}_${PLATE}_${WELL}.sam))
}

export -f build_fastq

function build_header() {
    local RUN=${1}
    local GZ_FILE=${RUN}.sam.gz
    local ZIP_FILE=${RUN}.sam.zip
    # get cross ID
    local CROSS=$(zip -sf ${ZIP_FILE} | head -n 2 | tail -n 1 | xargs echo -n | cut -d "_" -f 1)
    # Extract full header
    zcat ${GZ_FILE} | head -n 30000 | grep "^@" > ${RUN}.header.full.txt

    # Include all @HD lines
    grep "^@HD" ${RUN}.header.full.txt > ${RUN}.header.txt

    # Include all @SQ lines
    grep "^@SQ" ${RUN}.header.full.txt >> ${RUN}.header.txt

    # Include all @RG lines that include $RUN
    grep "^@RG" ${RUN}.header.full.txt | grep "ID:$CROSS"  >> ${RUN}.header.txt

    # Include all @PG lines that include $RUN
    grep "^@PG" ${RUN}.header.full.txt | grep "$CROSS" >> ${RUN}.header.txt

    # remove initial tmp header
    rm ${RUN}.header.full.txt
}

export -f build_header

## Reconstruct header

RUN=SRR9330831
RUN=SRR9330810

samfile_list=($(get_samfile_list $RUN))

for samfile in ${samfile_list[@]}; do
    build_fastq ${RUN} ${samfile}
done &

ls /data/SBGE/cory/${RUN}*.fastq | zip -@ -j ${RUN}.fastq.zip &

sinteractive --mem=8G --gres=lscratch:200

cd /lscratch/${SLURM_JOB_ID}
module load samtools

RUN=SRR9330831
RUN=SRR9330810

samfile_list=($(get_samfile_list $RUN))

for samfile in ${samfile_list[@]}; do
    build_fastq ${RUN} ${samfile}
done &

Now have `SRR9330831.fastq.zip` and `SRR9330810.fastq.zip`




Make blast DBs

for file in *.fasta; do makeblastdb -dbtype nucl -in ${file}; done

# Submit Jobs
```bash
sbatch ~/pacbio-yeast-genomes/src/map_cross.slurm SRR9330831 CBS2888 CLIB219 chrVII
sbatch ~/pacbio-yeast-genomes/src/map_cross.slurm SRR9330831 CBS2888 CLIB219 chrX
sbatch ~/pacbio-yeast-genomes/src/map_cross.slurm SRR9330831 CBS2888 CLIB219 chrXV

sbatch ~/pacbio-yeast-genomes/src/map_cross.slurm SRR9330810 RM YPS163 chrVII
```


awk $11~"CBS2888_pacbio_chrX$" out.snps > out.fixed.snps

# Plot
```R
library(data.table)
library(ggplot2)

files <- list.files(pattern="*.call")
dat <- fread(files[1])
dat[, "#CHROM" := NULL]
setkey(dat, POS)

for(file in files[2:length(files)]) {
    dat.tmp <- fread(file, select=c(2,5))
    dat.tmp <- dat.tmp[!duplicated(POS)]
    setkey(dat.tmp, POS)
    dat <- merge(dat, dat.tmp, all=TRUE)
}

dat.long <- melt(dat, measure.vars=colnames(dat)[4:length(colnames(dat))], variable.name="sample",value.name="haplotype")[!is.na(haplotype)]

ggplot(dat.long, aes(x=POS, y=sample, color=haplotype, fill=haplotype)) + geom_tile()

calculate_ld <- function(ab, a, b) {
    (ab - a*b)^2 / (a*(1-a)*b*(1-b))
}

dat.long[ALT != .]

get_ld <- function(pos1, pos2, DT) {
    DT.tmp <- DT[POS %in% c(pos1, pos2)]
    bothalleles <- DT.tmp[, .N, by=sample][N==2][,sample]
    DT.tmp <- DT.tmp[sample %in% bothalleles]
    N <- length(bothalleles)
    ab_N <- nrow(DT.tmp[, .N, by=list(sample, haplotype)][haplotype==0 & N==2])
    ab <- ab_N/N
    a <- nrow(DT[POS==pos1 & haplotype == 0])/N
    b <- nrow(DT[POS==pos2 & haplotype == 0])/N
    r2 <- (ab - a*b)^2 / (a*(1-a)*b*(1-b))
    return(data.table(pos1, pos2, r2))
}

positions <- sort(unique(dat.long$POS))[1:100]

out <- foreach(idx1=1:(length(positions)-1), .combine=rbind) %do% {
    foreach(idx2=idx1+1 : length(positions), .combine=rbind) %do% {
        pos_1 <- positions[idx1]
        pos_2 <- positions[idx2]
        get_ld(pos_1, pos_2, dat.long)
    }
}

ggplot(out[!is.na(pos2)], aes(x=pos1, y=pos2, fill=r2))+ geom_tile()

```

python ~/pacbio-yeast-genomes/src/write_ped.py genotypes.mat > genotypes.ped
python ~/pacbio-yeast-genomes/src/write_map.py genotypes.mat > genotypes.map


plink --file genotypes --map3 --missing-genotype N
plink --bfile plink --r2 inter-chr --ld-window-r2 0

library(data.table)
library(ggplot2)
library(viridis)
dat <- fread(plink.ld)
ggplot(dat, aes(x=BP_A, y=BP_B, fill=R2)) + geom_tile() + theme(plot.background = element_rect(fill = "black")) +
theme(panel.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()
        ) +
scale_fill_viridis()

```python
import sys

file1 = sys.argv[1]
file2 = sys.argv[2]

with open(file1, r) as infile:
    file1_text = [x.strip() for x in infile.readlines()]

with open(file2, r) as infile:
    file2_text = [x.strip() for x in infile.readlines()]

for i in range(len(file1_text)):
    firstGenos = file1_text[i]
    secondGenos = \t.join(file2_text[i].split()[6:])
    output = firstGenos + \t + secondGenos
    print(output)
```

~/pacbio-yeast-genomes/mergePed.py chrX.ped chr15.ped  > merged.ped
cat chrX.map chr15.map > merged.map

plink --file merged --map3 --missing-genotype N
plink --bfile plink --r2 inter-chr --ld-window-r2 0

# same chrom
library(data.table)
library(ggplot2)
library(viridis)
dat <- fread(plink.ld)
ggplot(dat[, aes(x=BP_A, y=BP_B, fill=R2)) + geom_tile() + theme(plot.background = element_rect(fill = "black")) +
theme(panel.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()
        ) +
scale_fill_viridis()

# cross chrom
library(data.table)
library(ggplot2)
library(viridis)
dat <- fread(plink.ld)

g1 <- ggplot(dat[CHR_A =="23" & CHR_B == "23"], aes(x=BP_A, y=BP_B, fill=R2)) + 
geom_rect(data=NULL,aes(xmin=-Inf,xmax=Inf,ymin=-Inf,ymax=Inf), fill="black") +
geom_tile() +
theme(panel.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()
        ) +
scale_fill_viridis() +
labs(x="chrX", y="chrX")


g2 <- ggplot(dat[CHR_A =="15" & CHR_B == "15"], aes(x=BP_A, y=BP_B, fill=R2)) + 
geom_rect(data=NULL,aes(xmin=-Inf,xmax=Inf,ymin=-Inf,ymax=Inf), fill="black") +
geom_tile() +
theme(panel.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()
        ) +
scale_fill_viridis() +
labs(x="chrXV", y="chrXV")


g3 <- ggplot(dat[CHR_A =="15" & CHR_B == "23"], aes(x=BP_A, y=BP_B, fill=R2)) + 
geom_rect(data=NULL,aes(xmin=-Inf,xmax=Inf,ymin=-Inf,ymax=Inf), fill="black") +
geom_tile() +
theme(panel.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()
        ) +
scale_fill_viridis() +
labs(x="chrXV", y="chrX")

ggsave(g1, file="chrX_chrX.png")
ggsave(g2, file="chrXV_chrXV.png")
ggsave(g3, file="chrXV_chrX.png")

# MALR


python3 ~/pacbio-yeast-genomes/src/write_ped.py genotypes.mat > genotypes.ped
python3 ~/pacbio-yeast-genomes/src/write_ped.py MALR.mat | sed s/SRR9330831_//g > MALR.ped

python3 ~/pacbio-yeast-genomes/src/write_map.py genotypes.mat > genotypes.map
python3 ~/pacbio-yeast-genomes/src/write_map.py MALR.mat > MALR.map


python3 ~/pacbio-yeast-genomes/mergePed.py genotypes.ped MALR.ped > merged.ped
cat genotypes.map MALR.map > merged.map

plink --file merged --map3 --missing-genotype N
plink --bfile plink --r2 inter-chr --ld-window-r2 0


#

Get SRA link for wget

# M22 x CLIB219
run="SRR9330832"
curl -v https://locate.ncbi.nlm.nih.gov/sdl/2/retrieve?acc=${run}
wget -O ${run}.sra https://sra-pub-run-odp.s3.amazonaws.com/sra/${run}/${run}

# BY x RM
run="SRR9330837"
wget -O ${run}.sra  https://sra-pub-run-odp.s3.amazonaws.com/sra/${run}/${run}

# YJM145 x CLIB413
run="SRR9330813"
wget -O ${run}.sra  https://sra-pub-run-odp.s3.amazonaws.com/sra/${run}/${run}

# CLIB413 x YJM978
run="SRR9330815"
wget -O ${run}.sra  https://sra-pub-run-odp.s3.amazonaws.com/sra/${run}/${run}

# PW5 x 273641
run="SRR9330830"
wget -O ${run}.sra  https://sra-pub-run-odp.s3.amazonaws.com/sra/${run}/${run}

module load sratoolkit

sam-dump



sbatch ~/pacbio-yeast-genomes/src/split_sam.slurm SRR9330832

```bash

cat <<EOT > runs.txt
YJM981_CBS2888 SRX6097475 SRR9330808 3004
CBS2888_CLIB219 SRX6097452 SRR9330831 3043
CLIB219_M22 SRX6097451 SRR9330832 3028
M22_BY SRX6097472 SRR9330811 375
BY_RM SRX6097446 SRR9330837 A
RM_YPS163 SRX6097473 SRR9330810 376
YPS163_YJM145 SRX6097447 SRR9330836 B
YJM145_CLIB413 SRX6097470 SRR9330813 377
CLIB413_YJM978 SRX6097468 SRR9330815 393
YJM978_YJM454 SRX6097471 SRR9330812 381
YJM454_YPS1009 SRX6097450 SRR9330833 3008
YPS1009_I14 SRX6097469 SRR9330814 2999
I14_Y10 SRX6097466 SRR9330817 3000
Y10_PW5 SRX6097467 SRR9330816 3001
PW5_273641 SRX6097453 SRR9330830 3049
273641_YJM981 SRX6097474 SRR9330809 3003
EOT

while read cross srx run cramid; do
    sbatch ~/pacbio-yeast-genomes/src/split_sam.slurm ${run} ${cramid}
done < runs.txt

```

## Analyses
# 1 Linkage between chrX and chrXV
mostly done
# 2 Linkage of MALR genes to chrVII
mostly done

# 3 Linkage of MALR genes on chrXI
# BY x RM 01 - 11              SRR9330837  chrXI
sbatch ~/pacbio-yeast-genomes/src/map_cross.slurm SRR9330837 BY RM chrXI

# RM x YPS163           SRR9330810  chrXI
sbatch ~/pacbio-yeast-genomes/src/map_cross.slurm SRR9330810 RM YPS163 chrXI

    # CLIB413 x YJM145      SRR9330813  chrXI
sbatch ~/pacbio-yeast-genomes/src/map_cross.slurm SRR9330813 CLIB413 YJM145 chrXI

    # CLIB413 x YJM978      SRR9330815  chrXI
sbatch ~/pacbio-yeast-genomes/src/map_cross.slurm SRR9330815 CLIB413 YJM978 chrXI


# 4 Linkage of SGE1 copies to chrVII and chrXI
26x27   RM x YPS163         SRR9330810  chrXI
sbatch ~/pacbio-yeast-genomes/src/map_cross.slurm SRR9330810 RM YPS163 chrXI

36x35   273614 x PW5        SRR9330830  chrXI
sbatch ~/pacbio-yeast-genomes/src/map_cross.slurm SRR9330830 273614 PW5 chrXI


36x37   273614 x YJM981     SRR9330809  chrXI
sbatch ~/pacbio-yeast-genomes/src/map_cross.slurm SRR9330809 273614 YJM981 chrXI

39x38   CBS2888 x CLIB219   SRR9330831  chrVII
sbatch ~/pacbio-yeast-genomes/src/map_cross.slurm SRR9330831 CBS2888 CLIB219 chrVII

39x24   CLIB219 x M22       SRR9330832  chrVII
sbatch ~/pacbio-yeast-genomes/src/map_cross.slurm SRR9330832 CLIB219 M22 chrVII

# 5 Linkage of SUC genes to chrI, chrIII, and chrIX.


31x30   YJM978 x YJM454     SRR9330812  chrIX
sbatch ~/pacbio-yeast-genomes/src/map_cross.slurm SRR9330812 YJM454 YJM978 chrIX

31x32   YJM454 x YPS1009    SRR9330833  chrIX
sbatch ~/pacbio-yeast-genomes/src/map_cross.slurm SRR9330833 YJM454 YPS1009 chrIX

34x33   I14 x Y10           SRR9330817  chrI + chrIII
34x35   Y10 x PW5           SRR9330816  chrI + chrIII

## SUC genes


sbatch ~/pacbio-yeast-genomes/src/map_cross.slurm SRR9330812 YJM454 YJM978 chrIX
sbatch ~/pacbio-yeast-genomes/src/map_cross.slurm SRR9330833 YJM454 YPS1009 chrIX
sbatch ~/pacbio-yeast-genomes/src/map_cross.slurm SRR9330817 Y10 I14 chrI
sbatch ~/pacbio-yeast-genomes/src/map_cross.slurm SRR9330817 Y10 I14 chrIII
sbatch ~/pacbio-yeast-genomes/src/map_cross.slurm SRR9330816 Y10 PW5 chrI
sbatch ~/pacbio-yeast-genomes/src/map_cross.slurm SRR9330816 Y10 PW5 chrIII


YJM981_CBS2888 SRX6097475 SRR9330808 3004
CBS2888_CLIB219 SRX6097452 SRR9330831 3043
CLIB219_M22 SRX6097451 SRR9330832 3028
M22_BY SRX6097472 SRR9330811 375
BY_RM SRX6097446 SRR9330837 A   ### BAD
RM_YPS163 SRX6097473 SRR9330810 376
YPS163_YJM145 SRX6097447 SRR9330836 B ### BAD
YJM145_CLIB413 SRX6097470 SRR9330813 377
CLIB413_YJM978 SRX6097468 SRR9330815 393
YJM978_YJM454 SRX6097471 SRR9330812 381
YJM454_YPS1009 SRX6097450 SRR9330833 3008
YPS1009_I14 SRX6097469 SRR9330814 2999
I14_Y10 SRX6097466 SRR9330817 3000
Y10_PW5 SRX6097467 SRR9330816 3001
PW5_273641 SRX6097453 SRR9330830 3049
273641_YJM981 SRX6097474 SRR9330809 3003


```bash
cat <<EOT > MS-strains.txt
MSY24	M22
MSY25	BY
MSY26	RM
MSY27	YPS163
MSY28	YJM145
MSY29	CLIB413
MSY30	YJM978
MSY31	YJM454
MSY32	YPS1009
MSY33	I14
MSY34	Y10
MSY35	PW5
MSY36	273614
MSY37	YJM981
MSY38	CBS2888
MSY39	CLIB219
EOT
```

while read MSID realID; do
    sed "s/${MSID}/${realID}/g" ${MSID}.fasta > ${realID}.fasta
done

function splitFasta() {
    inFasta=${1}
    awk -F "(^>|\t| )" '{if($0 ~ /^>/) {s=$2".fasta";  print ">"$2 > s} else print > s}' ${inFasta}
}

export -f splitFasta

for file in *.fasta; do
    splitFasta ${file}
done

ls *.fasta | zip -@ pacbio-assemblies.zip      


# SUC genes
>31-chrIX-SUC_het30_het32
31x30   YJM978  YJM454  SRR9330812
31x32   YJM454  YPS1009 SRR9330833

>34-chrI_AND_chrIII-SUC_het33_het35
34x33 Y10 x I14     SRR9330817
34x35  Y10 x PW5    SRR9330816

