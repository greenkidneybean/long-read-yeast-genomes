# SGE-1 mapping

| Strain1   | Strain2   | SRA Run    | Chr     | Fasta                |
|:----------|:----------|:-----------|:--------|:---------------------|
|  RM       | YPS163    | SRR9330810 | chrXI   | chrXI-SGE1.fasta     |
|  273614   | PW5       | SRR9330830 | chrXI   | chrXI_659119.fasta   |
|  273614   | YJM981    | SRR9330809 | chrXI   | chrXI_659119.fasta   |
|  CLIB219  | CBS2888   | SRR9330831 | chrVII  | chrVII_1068298.fasta |
|  CLIB219  | M22       | SRR9330832 | chrVII  | chrVII_1068298.fasta |

## Defining diagnostic FASTA sequences
```bash
mkdir -p ./fastas
cat <<EOT > ./fastas/chrXI-SGE1.fasta
>chrXI-SGE1
ATTCGGGTATCATTATTTGCTTTTTTACCGTGGGCCCTATCTTATTGTTACTATTTTGTGCTTACGACTTTCATTTTCTG
tCATTATcGGGGCTTCAt
TATGACAACAAGCGGATCAAACCGTTACTGACATGGAATATTGCCTCAAATTGTGGCATATTTACAAGCTCCATAACAGG
EOT

cat <<EOT > ./fastas/chrXI-659119.fasta
>chrXI-659119
TCAATTCTCCCTGGAATAGCTTTTGGTAGTATTTTCCAAGCAACGTTATTAAGCTCCCAGGTGCAGATAACATCAGACGA
g
CCAGACTTTCAAAACAAGTTTATTGAAGTCACAGCTTTCAACTCGTTCGCCAAATCCTTGGGCTTTGCGTTTGGAGGGAA
EOT

cat <<EOT > ./fastas/chrVII-1068298.fasta
>chrVII-1068298
AAAAAGCCTACATTAGCGAGTATACATCTTTGGGAACTATCAATTCCAGCTATGATTGCAACTATGGCCATAGCATATCT
gAATTCAAAATATGGgATt
ATCAAACCGGCAATTGTTTTTGGTGTGCTTTGTGGGATTGTTGGATCTGGTTTATTTACGCTAATCAATGGAGAACTCTC
EOT
```

## Call Genotypes
```bash
source ../src/fun.sh
call_genos 'RM'         'YPS163'    'chrXI'     'SRR9330810'
call_genos '273614'     'PW5'       'chrXI'     'SRR9330830'
call_genos '273614'     'YJM981'    'chrXI'     'SRR9330809'
call_genos 'CLIB219'    'CBS2888'   'chrVII'    'SRR9330831'
call_genos 'CLIB219'    'M22'       'chrVII'    'SRR9330832'
```


## Call phenotypes
```bash
source ../src/fun.sh
call_phenos 'RM'        'YPS163'    'SRR9330810'    'chrXI-SGE1.fasta'
call_phenos '273614'    'PW5'       'SRR9330830'    'chrXI-659119.fasta'
call_phenos '273614'    'YJM981'    'SRR9330809'    'chrXI-659119.fasta'
call_phenos 'CLIB219'   'CBS2888'   'SRR9330831'    'chrVII-1068298.fasta'
call_phenos 'CLIB219'   'M22'       'SRR9330832'    'chrVII-1068298.fasta'
```

## Association testing
```bash
source ../src/fun.sh
qtl_map 'RM'        'YPS163'    'chrXI-SGE1.fasta'      'chrXI'
qtl_map '273614'    'PW5'       'chrXI-659119.fasta'    'chrXI'
qtl_map '273614'    'YJM981'    'chrXI-659119.fasta'    'chrXI'
qtl_map 'CLIB219'   'CBS2888'   'chrVII-1068298.fasta'  'chrVII'
qtl_map 'CLIB219'   'M22'       'chrVII-1068298.fasta'  'chrVII'
```

## Plotting
```bash
source ../src/fun.sh
plot_association 'RM'        'YPS163'    'chrXI-SGE1.fasta'      'chrXI'
plot_association '273614'    'PW5'       'chrXI-659119.fasta'    'chrXI'
plot_association '273614'    'YJM981'    'chrXI-659119.fasta'    'chrXI'
plot_association 'CLIB219'   'CBS2888'   'chrVII-1068298.fasta'  'chrVII'
plot_association 'CLIB219'   'M22'       'chrVII-1068298.fasta'  'chrVII'
```

## Upload to onedrive
```bash
module load rclone
rclone copy --max-depth 1 --include "*.pdf" ./plots nihonedrive:/Data/pacbio_genomes/04_SGE1/
```
