# conda activate alignparse-environment

import argparse
import pandas as pd
from Bio import SeqIO
import re
import os, os.path

# parse arguments
parser = argparse.ArgumentParser()
parser.add_argument("-orfs", help="ORFfinder .fasta output file of all ORFs", type=str)
parser.add_argument("-blast", help="Blast .tsv output file", type=str)
parser.add_argument("-out", help="Output file.fasta", type=str)
args = parser.parse_args()

def safe_open_w(path):
    ''' Open "path" for writing, creating any parent directories as needed.
    '''
    os.makedirs(os.path.dirname(path), exist_ok=True)
    return open(path, 'w')

# Python check if file exists argparse

cols=headers=['query','target','pident','length','mismatch','gapopen','qstart','qend','tstart','tend','evalue','bitscore']

# need to get fasta headers to a list

# set of raw orfs gene names
raw_orfs = []
with open(args.orfs, "r") as f:
    for record in SeqIO.parse(f, "fasta"):
        raw_orfs.append(record.description.split(' ')[0])
raw_orfs = set(raw_orfs)

# set of orfs present in reference blast
df=pd.read_csv(args.blast,sep='\t', names=cols)
blast_orfs = set(df['query'].tolist())

# get unique orfs with no blast hit in reference
unique_orfs = raw_orfs - blast_orfs

fasta_sequences = SeqIO.parse(open(args.orfs),'fasta')

with safe_open_w(args.out) as f:
    for seq in fasta_sequences:
        if seq.description.split(' ')[0] in unique_orfs:
            SeqIO.write([seq], f, "fasta")
