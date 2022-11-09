#!/usr/bin/env python3
'''only includes snps from <snps_filename> that are 
at least <dist_threshold> bp away from alignment mismatch'''

import sys

snps_filename = sys.argv[1]
dist_threshold = int(sys.argv[2])

snps = []

with open(snps_filename, 'r', encoding='utf-8') as infile:
    for line in infile:
        split_line = line.strip().split()
        ref_pos, ref, alt, alt_pos, distance = split_line[:5]
        strand = int(split_line[9])
        ref_name = split_line[10]
        #alt_name = split_line[11]
        if strand != -1:
            ref_pos = int(ref_pos) - 1
            alt_pos = int(alt_pos) - 1
            distance = int(distance)
            if distance >= dist_threshold:
                start = ref_pos - distance + 2
                end = ref_pos + distance
                snps.append((ref_name, start, end))

prev_ref_name = ''
prev_start = 0
prev_end = 0

for snp in snps:
    ref_name, start, end = snp
    if prev_ref_name == '':
        prev_ref_name, prev_start, prev_end = snp
        continue
    if ref_name != prev_ref_name and prev_ref_name != '':
        print('this should not happen')
        exit(1)
    if start > prev_end:    # non-overlapping with previous
        print(prev_ref_name, prev_start, prev_end, sep='\t')
        prev_ref_name, prev_start, prev_end = snp
    elif start <= prev_end: # overlapping with previous
        prev_end = end      # update prev_end to current end
print(prev_ref_name, prev_start, prev_end, sep='\t')
