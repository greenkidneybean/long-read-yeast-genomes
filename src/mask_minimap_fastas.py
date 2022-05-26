#!/usr/bin/env python

import sys

fasta_filename = sys.argv[1]
start = int(sys.argv[2])
end = int(sys.argv[3])

def import_fasta(f):
    seqs = {}
    with open(f, 'r', encoding='utf-8') as infile:
        entries = infile.read().split('>')[1:]
        for i in entries:
            header = i[1:1000].split('\n')[0].strip()
            nucleotides = ''.join(i.split('\n')[1:])
            seqs[header] = nucleotides

fasta = import_fasta(fasta_filename)

def mask_fasta(seq, start, end):
    '''replaces nucleotides outside the (start, end) range with N'''
masked_fasta = 'N'*start + 'N'*end
