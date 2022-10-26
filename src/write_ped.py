#!/usr/bin/env python

import sys
infile = sys.argv[1]

with open(infile, 'r', encoding='utf-8') as infile:
    for line in infile:
        if line.startswith('\t'):
            continue
        if not 'bin' in line:
            #line = line.replace('0','R')
            #line = line.replace('1','A')
            #line = line.replace('NA','N')
            text = line.strip().split()
            identity = text[0]
            ids = '\t'.join([identity] * 2 + ['0'] + ['0'])
            sex = '1'
            phenotype = '-9'
            genotypes = '\t'.join([f"{x}\t{x}" for x in text[1:]])
            genotypes = genotypes.replace('0','R')
            genotypes = genotypes.replace('1','A')
            genotypes = genotypes.replace('NA','N')
            print('\t'.join([ids, sex, phenotype, genotypes]))
