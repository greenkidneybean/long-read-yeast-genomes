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