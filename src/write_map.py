#!/usr/bin/env python
 
import sys
infile = sys.argv[1]

with open(infile, 'r', encoding='utf-8') as infile:
    for line in infile:
        if line.startswith('\t'):
            bins = [x.split('_')[1] for x in line.strip().split()]


chromosome = '15'
id = 0

for bin in bins:
    id += 1
    idstr = str(id)
    print('\t'.join([chromosome, idstr, bin]))