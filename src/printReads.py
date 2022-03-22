#!/usr/bin/env python3
# script currently unused
import argparse
import os

parser = argparse.ArgumentParser()

parser.add_argument('samfile', type=str,
                    help='File stem for gene 1. Ensure proper naming of <GENE1_FILESTEM>.fasta')
parser.add_argument("--debug", 
                        help="Assist with debugging by increasing output verbosity",
                        action="store_true")

args = parser.parse_args()

assert os.path.exists(args.samfile), "Error: samfile %s does not exist" % (args.samfile)

with open(args.samfile, 'r') as infile:
    for line in infile:
        read, score, RG = line.strip().split()
        RG = '_'.join(RG.split(":")[2].split("_")[1:])
        read = "@%s\n%s\n+\n%s" % (RG,read,score)
        print(read)