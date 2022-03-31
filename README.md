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


