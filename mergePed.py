#!/usr/bin/env python

import sys

file1 = sys.argv[1]
file2 = sys.argv[2]

with open(file1, 'r') as infile:
    file1_text = [x.strip() for x in infile.readlines()]

with open(file2, 'r') as infile:
    file2_text = [x.strip() for x in infile.readlines()]

for i in range(len(file1_text)):
    firstGenos = file1_text[i]
    secondGenos = '\t'.join(file2_text[i].split()[6:])
    output = firstGenos + '\t' + secondGenos
    print(output)

