#!/usr/bin/env python
'''only includes snps from <snps_filename> that are 
at least <dist_threshold> bp away from alignment mismatch'''

import sys

snps_filename = sys.argv[1]
dist_threshold = int(sys.argv[2])

with open(snps_filename, 'r', encoding='utf-8') as infile:
    for line in infile:
        split_line = line.strip().split()
        ref_pos, ref, alt, alt_pos, distance = split_line[:5]
        strand = int(split_line[9])
        ref_name = split_line[10]
        #alt_name = split_line[11]
        if strand != -1:
            distance = int(distance)
            if distance >= dist_threshold:
                print(ref_name, ref_pos, ref, alt, sep="\t")
