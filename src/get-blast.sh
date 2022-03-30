#!/usr/bin/env bash

wget -O blast.tar.gz wget https://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/2.13.0/ncbi-blast-2.13.0+-x64-linux.tar.gz
tar -zxvf blast.tar.gz
mv ncbi-blast-2.13.0+ blast