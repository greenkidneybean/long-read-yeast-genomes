# SUC mapping


## Defining diagnostic FASTA sequences
```bash
mkdir -p ./fastas
cat <<EOT > ./fastas/chrIX-SUC.fasta
>chrIX-SUC
TCCATGTCTTTGGTCCGCAAGTTTTCTTTGAACACTGAATATCAAGCTAATCCAGAGACTGAATTGATCAATTTGAAAGC
tgaaccaatactgaacattagtaatgccggtccatggttacactttgccagcaactccactttgacgaaagccaattcgt
tcagcgtagatttgagtaattctaccggcacat
TAGAGTTTGAGTTGGTTTACGCTGTTAACACCACACAATCTGTTTCTAAAAGTGTCTTTTCAGACTTATCTCTTTGGTTT
EOT

cat <<EOT > ./fastas/chrI-chrIII-SUC.fasta
>chrI-chrIII-SUC
ACGACCTGAAGTCCTGGAAGCTAGAATCTGCATTTGCCAATGAAGGTTTTTTAGGTTACCAATATGAGTGTCCCGGTTTG
a
TCGAAGTCCCAACTGAGCAAGATCCTTCTAAATCGCATTGGGTCATGTTTATTTCTATCAACCCAGGTGCACCTGCTGGC
NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN
NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN
TATCAAGCTAATCCAGAGACTGAATTGATCAATTTGAAAGCTGAACCAATACTGAACATTAGTAATGCCGGTCCATGGTT
a
CACTTTGCCAGCAACTCCACTTTGACGAAAGCCAATTCGTTCAGCGTAGATTTGAGTAATTCTACCGGCACATTAGAGTT
NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN
NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN
GAAGATCCTGAAGAATATTTAAGAATGGGTTTTGAAGCCAGTGCTTCTTCCTTCTTTTTGGACCGTGGTAACTCTAAGGT
c
AAGTTTGTCAAGGAGAACCCATATTTCACAAACAGAATGTCTGTCAACAACCAACCATTCAAGTCTGAGAACGACCTAAG
EOT
```

## Call Genotypes
```bash
source ../src/fun.sh
call_genos 'YJM454' 'YJM978'    'chrIX'     'SRR9330812'
call_genos 'YJM454' 'YPS1009'   'chrIX'     'SRR9330833'
call_genos 'Y10'    'I14'       'chrI'      'SRR9330817'
call_genos 'Y10'    'PW5'       'chrI'      'SRR9330816'
call_genos 'Y10'    'I14'       'chrIII'    'SRR9330817'
call_genos 'Y10'    'PW5'       'chrIII'    'SRR9330816'
```

## Call phenotypes
```bash
source ../src/fun.sh
call_phenos 'YJM454'    'YJM978'    'SRR9330812'    'chrIX-SUC.fasta'
call_phenos 'YJM454'    'YPS1009'   'SRR9330833'    'chrIX-SUC.fasta'
call_phenos 'Y10'       'I14'       'SRR9330817'    'chrI-chrIII-SUC.fasta'
call_phenos 'Y10'       'PW5'       'SRR9330816'    'chrI-chrIII-SUC.fasta'
```

## Association testing
```bash
source ../src/fun.sh
qtl_map 'YJM454'    'YJM978'    'chrIX-SUC.fasta'       'chrIX'
qtl_map 'YJM454'    'YPS1009'   'chrIX-SUC.fasta'       'chrIX'
qtl_map 'Y10'       'I14'       'chrI-chrIII-SUC.fasta' 'chrI'
qtl_map 'Y10'       'PW5'       'chrI-chrIII-SUC.fasta' 'chrI'
qtl_map 'Y10'       'I14'       'chrI-chrIII-SUC.fasta' 'chrIII'
qtl_map 'Y10'       'PW5'       'chrI-chrIII-SUC.fasta' 'chrIII'
```

## Plotting
```bash
source ../src/fun.sh
plot_association 'YJM454'    'YJM978'    'chrIX-SUC.fasta'          'chrIX'
plot_association 'YJM454'    'YPS1009'   'chrIX-SUC.fasta'          'chrIX'
plot_association 'Y10'       'I14'       'chrI-chrIII-SUC.fasta'    'chrI'
plot_association 'Y10'       'PW5'       'chrI-chrIII-SUC.fasta'    'chrI'
plot_association 'Y10'       'I14'       'chrI-chrIII-SUC.fasta'    'chrIII'
plot_association 'Y10'       'PW5'       'chrI-chrIII-SUC.fasta'    'chrIII'
```

## Upload to onedrive
```bash
module load rclone
rclone copy  --max-depth 1 --include "*.pdf" ./plots nihonedrive:/Data/pacbio_genomes/05_SUC
```
