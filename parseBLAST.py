#!/usr/bin/env python3

import sys

#filename = sys.argv[1]
filename = 'aln.maf'
with open(filename, 'r') as infile:
    text = infile.read().strip().split("//")[:-1]

alignments = {}


for alignment in text:
    info, read1, gaps, read2 = alignment.strip().split('\n')
    splitInfo = info.split('\t')
    seq = splitInfo[9]
    matchLen = len(seq)
    query = splitInfo[0].strip('>')
    if query not in alignments or matchLen > alignments[query][0]:
        alignments[query] = [matchLen, seq, read1, read2]

alignments['tig00000087j113rc'][1]
alignments['tig00000001']

764288