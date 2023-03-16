# need to make sure the file exists

import argparse
from Bio.Seq import Seq
from Bio import SeqIO

parser = argparse.ArgumentParser()
parser.add_argument("input", help='input fasta file to be translated')
parser.add_argument("output", help='output fasta file to be created')
args = parser.parse_args()

input_file = "example.fasta"
output_file = open(args.output,'w')

for record in SeqIO.parse(args.input, "fasta"):
    output_file.write(">" + record.id + "\n")
    output_file.write(str(Seq(record.seq).translate()) + "\n")

output_file.close()
