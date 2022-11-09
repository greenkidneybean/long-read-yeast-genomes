# Mapping of chrVII MALR genes

## Call genotypes
```bash
source ../src/fun.sh
call_genos  'YPS163'      'YJM145'      'chrVII'    'SRR9330836'
call_genos  'M22'         'BY'          'chrVII'    'SRR9330811'
call_genos  'BY'          'RM'          'chrVII'    'SRR9330837'
call_genos  'M22'         'CLIB219'     'chrVII'    'SRR9330832'
call_genos  'RM'          'YPS163'      'chrVII'    'SRR9330810'
call_genos  'YJM145'      'CLIB413'     'chrVII'    'SRR9330813'
call_genos  'YJM978'      'CLIB413'     'chrVII'    'SRR9330815'
call_genos  'YJM978'      'YJM454'      'chrVII'    'SRR9330812'
call_genos  'YPS1009'     'YJM454'      'chrVII'    'SRR9330833'
call_genos  'YPS1009'     'I14'         'chrVII'    'SRR9330814'
call_genos  'I14'         'Y10'         'chrVII'    'SRR9330817'
call_genos  'PW5'         'Y10'         'chrVII'    'SRR9330816'
call_genos  'PW5'         '273614'      'chrVII'    'SRR9330830'
call_genos  'YJM981'      'CBS2888'     'chrVII'    'SRR9330808'
call_genos  'CBS2888'     'YJM981'      'chrVII'    'SRR9330808'
call_genos  'CLIB219'     'CBS2888'     'chrVII'    'SRR9330831'
```


## Call phenotypes
```bash
source ../src/fun.sh
call_phenos 'M22'         'BY'          'M22-tf2.fasta'     'SRR9330811'
call_phenos 'M22'         'CLIB219'     'M22-tf2.fasta'     'SRR9330832'
call_phenos 'RM'          'YPS163'      'RM-tf2.fasta'      'SRR9330810'
call_phenos 'YPS163'      'YJM145'      'YPS163-tf1.fasta'  'SRR9330836'
call_phenos 'YPS163'      'RM'          'YPS163-tf2.fasta'  'SRR9330810'
call_phenos 'YPS163'      'YJM145'      'YPS163-tf2.fasta'  'SRR9330836'
call_phenos 'YPS163'      'RM'          'YPS163-tf3.fasta'  'SRR9330810'
call_phenos 'YPS163'      'YJM145'      'YPS163-tf3.fasta'  'SRR9330836'
call_phenos 'YPS163'      'RM'          'YPS163-tf4.fasta'  'SRR9330810'
call_phenos 'YJM145'      'CLIB413'     'YJM145-tf.fasta'   'SRR9330813'
call_phenos 'CLIB413'     'YJM145'      'CLIB413-tf.fasta'  'SRR9330813'
call_phenos 'YJM978'      'CLIB413'     'YJM978-tf1.fasta'  'SRR9330815'
call_phenos 'YJM978'      'YJM454'      'YJM978-tf1.fasta'  'SRR9330812'
call_phenos 'YJM978'      'YJM454'      'YJM978-tf2.fasta'  'SRR9330812'
call_phenos 'YJM978'      'CLIB413'     'YJM978-tf3.fasta'  'SRR9330815'
call_phenos 'YJM978'      'YJM454'      'YJM978-tf3.fasta'  'SRR9330812'
call_phenos 'YJM978'      'CLIB413'     'YJM978-tf4.fasta'  'SRR9330815'
call_phenos 'YJM454'      'YPS1009'     'YJM454-tf.fasta'   'SRR9330833'
call_phenos 'YPS1009'     'YJM454'      'YPS1009-tf1.fasta' 'SRR9330833'
call_phenos 'YPS1009'     'I14'         'YPS1009-tf1.fasta' 'SRR9330814'
call_phenos 'YPS1009'     'YJM454'      'YPS1009-tf2.fasta' 'SRR9330833'
call_phenos 'YPS1009'     'I14'         'YPS1009-tf2.fasta' 'SRR9330814'
call_phenos 'I14'         'YPS1009'     'I14-tf1.fasta'     'SRR9330814'
call_phenos 'I14'         'Y10'         'I14-tf1.fasta'     'SRR9330817'
call_phenos 'PW5'         'Y10'         'PW5-tf1.fasta'     'SRR9330816'
call_phenos 'PW5'         '273614'      'PW5-tf1.fasta'     'SRR9330830'
call_phenos 'PW5'         'Y10'         'PW5-tf2.fasta'     'SRR9330816'
call_phenos 'PW5'         '273614'      'PW5-tf2.fasta'     'SRR9330830'
call_phenos 'PW5'         'Y10'         'PW5-tf3.fasta'     'SRR9330816'
call_phenos 'PW5'         '273614'      'PW5-tf3.fasta'     'SRR9330830'
call_phenos '273614'      'PW5'         '273614-tf.fasta'   'SRR9330830'
call_phenos 'YJM981'      'CBS2888'     'YJM981-tf.fasta'   'SRR9330808'
call_phenos 'CBS2888'     'YJM981'      'CBS2888-tf.fasta'  'SRR9330808'
call_phenos 'CLIB219'     'CBS2888'     'CLIB219-tf1.fasta' 'SRR9330831'
call_phenos 'CLIB219'     'M22'         'CLIB219-tf2.fasta' 'SRR9330832'
call_phenos 'CLIB219'     'CBS2888'     'CLIB219-tf3.fasta' 'SRR9330831'
call_phenos 'CLIB219'     'M22'         'CLIB219-tf3.fasta' 'SRR9330832'
call_phenos 'CLIB219'     'CBS2888'     'CLIB219-tf4.fasta' 'SRR9330831'
call_phenos 'CLIB219'     'M22'         'CLIB219-tf4.fasta' 'SRR9330832'
```


## Association testing
```bash
source ../src/fun.sh
qtl_map 'M22'         'BY'          'M22-tf2.fasta'     'chrVII'
qtl_map 'M22'         'CLIB219'     'M22-tf2.fasta'     'chrVII'
qtl_map 'RM'          'YPS163'      'RM-tf2.fasta'      'chrVII'
qtl_map 'YPS163'      'YJM145'      'YPS163-tf1.fasta'  'chrVII'
qtl_map 'YPS163'      'RM'          'YPS163-tf2.fasta'  'chrVII'
qtl_map 'YPS163'      'YJM145'      'YPS163-tf2.fasta'  'chrVII'
qtl_map 'YPS163'      'RM'          'YPS163-tf3.fasta'  'chrVII'
qtl_map 'YPS163'      'YJM145'      'YPS163-tf3.fasta'  'chrVII'
qtl_map 'YPS163'      'RM'          'YPS163-tf4.fasta'  'chrVII'
qtl_map 'YJM145'      'CLIB413'     'YJM145-tf.fasta'   'chrVII'
qtl_map 'CLIB413'     'YJM145'      'CLIB413-tf.fasta'  'chrVII'
qtl_map 'YJM978'      'CLIB413'     'YJM978-tf1.fasta'  'chrVII'
qtl_map 'YJM978'      'YJM454'      'YJM978-tf1.fasta'  'chrVII'
qtl_map 'YJM978'      'YJM454'      'YJM978-tf2.fasta'  'chrVII'
qtl_map 'YJM978'      'CLIB413'     'YJM978-tf3.fasta'  'chrVII'
qtl_map 'YJM978'      'YJM454'      'YJM978-tf3.fasta'  'chrVII'
qtl_map 'YJM978'      'CLIB413'     'YJM978-tf4.fasta'  'chrVII'
qtl_map 'YJM454'      'YPS1009'     'YJM454-tf.fasta'   'chrVII'
qtl_map 'YPS1009'     'YJM454'      'YPS1009-tf1.fasta' 'chrVII'
qtl_map 'YPS1009'     'I14'         'YPS1009-tf1.fasta' 'chrVII'
qtl_map 'YPS1009'     'YJM454'      'YPS1009-tf2.fasta' 'chrVII'
qtl_map 'YPS1009'     'I14'         'YPS1009-tf2.fasta' 'chrVII'
qtl_map 'I14'         'YPS1009'     'I14-tf1.fasta'     'chrVII'
qtl_map 'I14'         'Y10'         'I14-tf1.fasta'     'chrVII'
qtl_map 'PW5'         'Y10'         'PW5-tf1.fasta'     'chrVII'
qtl_map 'PW5'         '273614'      'PW5-tf1.fasta'     'chrVII'
qtl_map 'PW5'         'Y10'         'PW5-tf2.fasta'     'chrVII'
qtl_map 'PW5'         '273614'      'PW5-tf2.fasta'     'chrVII'
qtl_map 'PW5'         'Y10'         'PW5-tf3.fasta'     'chrVII'
qtl_map 'PW5'         '273614'      'PW5-tf3.fasta'     'chrVII'
qtl_map '273614'      'PW5'         '273614-tf.fasta'   'chrVII'
qtl_map 'YJM981'      'CBS2888'     'YJM981-tf.fasta'   'chrVII'
qtl_map 'CBS2888'     'YJM981'      'CBS2888-tf.fasta'  'chrVII'
qtl_map 'CLIB219'     'CBS2888'     'CLIB219-tf1.fasta' 'chrVII'
qtl_map 'CLIB219'     'M22'         'CLIB219-tf2.fasta' 'chrVII'
qtl_map 'CLIB219'     'CBS2888'     'CLIB219-tf3.fasta' 'chrVII'
qtl_map 'CLIB219'     'M22'         'CLIB219-tf3.fasta' 'chrVII'
qtl_map 'CLIB219'     'CBS2888'     'CLIB219-tf4.fasta' 'chrVII'
qtl_map 'CLIB219'     'M22'         'CLIB219-tf4.fasta' 'chrVII'
```

## Plotting
```bash
source ../src/fun.sh
plot_association 'M22'         'BY'          'M22-tf2.fasta'     'chrVII'
plot_association 'M22'         'CLIB219'     'M22-tf2.fasta'     'chrVII'
plot_association 'RM'          'YPS163'      'RM-tf2.fasta'      'chrVII'
plot_association 'YPS163'      'YJM145'      'YPS163-tf1.fasta'  'chrVII'
plot_association 'YPS163'      'RM'          'YPS163-tf2.fasta'  'chrVII'
plot_association 'YPS163'      'YJM145'      'YPS163-tf2.fasta'  'chrVII'
plot_association 'YPS163'      'RM'          'YPS163-tf3.fasta'  'chrVII'
plot_association 'YPS163'      'YJM145'      'YPS163-tf3.fasta'  'chrVII'
plot_association 'YPS163'      'RM'          'YPS163-tf4.fasta'  'chrVII'
plot_association 'YJM145'      'CLIB413'     'YJM145-tf.fasta'   'chrVII'
plot_association 'CLIB413'     'YJM145'      'CLIB413-tf.fasta'  'chrVII'
plot_association 'YJM978'      'CLIB413'     'YJM978-tf1.fasta'  'chrVII'
plot_association 'YJM978'      'YJM454'      'YJM978-tf1.fasta'  'chrVII'
plot_association 'YJM978'      'YJM454'      'YJM978-tf2.fasta'  'chrVII'
plot_association 'YJM978'      'CLIB413'     'YJM978-tf3.fasta'  'chrVII'
plot_association 'YJM978'      'YJM454'      'YJM978-tf3.fasta'  'chrVII'
plot_association 'YJM978'      'CLIB413'     'YJM978-tf4.fasta'  'chrVII'
plot_association 'YJM454'      'YPS1009'     'YJM454-tf.fasta'   'chrVII'
plot_association 'YPS1009'     'YJM454'      'YPS1009-tf1.fasta' 'chrVII'
plot_association 'YPS1009'     'I14'         'YPS1009-tf1.fasta' 'chrVII'
plot_association 'YPS1009'     'YJM454'      'YPS1009-tf2.fasta' 'chrVII'
plot_association 'YPS1009'     'I14'         'YPS1009-tf2.fasta' 'chrVII'
plot_association 'I14'         'YPS1009'     'I14-tf1.fasta'     'chrVII'
plot_association 'I14'         'Y10'         'I14-tf1.fasta'     'chrVII'
plot_association 'PW5'         'Y10'         'PW5-tf1.fasta'     'chrVII'
plot_association 'PW5'         '273614'      'PW5-tf1.fasta'     'chrVII'
plot_association 'PW5'         'Y10'         'PW5-tf2.fasta'     'chrVII'
plot_association 'PW5'         '273614'      'PW5-tf2.fasta'     'chrVII'
plot_association 'PW5'         'Y10'         'PW5-tf3.fasta'     'chrVII'
plot_association 'PW5'         '273614'      'PW5-tf3.fasta'     'chrVII'
plot_association '273614'      'PW5'         '273614-tf.fasta'   'chrVII'
plot_association 'YJM981'      'CBS2888'     'YJM981-tf.fasta'   'chrVII'
plot_association 'CBS2888'     'YJM981'      'CBS2888-tf.fasta'  'chrVII'
plot_association 'CLIB219'     'CBS2888'     'CLIB219-tf1.fasta' 'chrVII'
plot_association 'CLIB219'     'M22'         'CLIB219-tf2.fasta' 'chrVII'
plot_association 'CLIB219'     'CBS2888'     'CLIB219-tf3.fasta' 'chrVII'
plot_association 'CLIB219'     'M22'         'CLIB219-tf3.fasta' 'chrVII'
plot_association 'CLIB219'     'CBS2888'     'CLIB219-tf4.fasta' 'chrVII'
plot_association 'CLIB219'     'M22'         'CLIB219-tf4.fasta' 'chrVII'

```


## Upload to onedrive
```bash
module load rclone
rclone copy  --dry-run --max-depth 1 --include "*chrVII.pdf" ./plots nihonedrive:/Data/pacbio_genomes/02_MALR_chrVII/
```


