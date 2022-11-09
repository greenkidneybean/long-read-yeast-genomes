#!/usr/bin/env python3
'''combines masking bed file + fasta file and returns masked fasta file'''
import sys
import argparse
import logging
import gzip
import fileinput
import time

class ExitOnExceptionHandler(logging.StreamHandler):
    '''logger class to exit on error or critical'''
    def emit(self, record):
        super().emit(record)
        if record.levelno in (logging.ERROR, logging.CRITICAL):
            raise SystemExit(-1)

logger = logging.getLogger(None)
logging.basicConfig(handlers=[ExitOnExceptionHandler()], level=logging.DEBUG)

def import_bed(f):
    '''imports bed file and returns dictionary of seq: (start, end)'''
    bed = {}
    with open(f, 'r', encoding='utf-8') as infile:
        for line in infile:
            seq, start, end = line.strip().split()
            start = int(start)
            end = int(end)
            if seq not in bed:
                bed[seq] = []
            bed[seq].append((start, end))
    return bed


def import_fasta(f):
    '''loads a fasta file and returns dictionary of header : seq'''
    seqs = {}
    try:
        with fileinput.hook_compressed(f, 'r') as infile:
            try:
                text = infile.read()
            except gzip.BadGzipFile as error:
                logger.error('File ends in .gz but cannot be decoded: %s', error)
                sys.exit(1)
            if f.endswith('.gz'):
                try:
                    text = text.decode()
                except(UnicodeDecodeError, AttributeError, gzip.BadGzipFile) as error:
                    logger.error('Cannot uncompress from gzip, even though \
                                  file ends in .gz: %s', error)
    except FileNotFoundError as error:
        logger.error('File %s does not exist, check spelling and location', f)
    entries = text.split('>')[1:]
    for i in entries:
        header = i[:1000].split('\n')[0].strip()
        nucleotides = ''.join(i.split('\n')[1:])
        seqs[header] = nucleotides
    return seqs

def get_mask_ranges(bed_ranges, seq):
    '''returns the negation of the included bed_ranges,
    i.e. everything that lies outside of those ranges.'''
    mask_ranges = [(1, bed_ranges[0][0] - 1)]
    mask_ranges +=[(bed_ranges[i][1] + 1, bed_ranges[i+1][0] - 1) for i in range(len(bed_ranges)-1)]
    mask_ranges += [(bed_ranges[-1][1], len(seq))]
    return mask_ranges


def wrap_fasta(seq, width=60):
    '''inserts newline every <width> characters'''
    return '\n'.join([seq[i:(i+width)] for i in range(0, len(seq), width)])


def mask_fasta(nucleotides, bed_retain_ranges):
    '''replaces nucleotides outside the (start, end) range with N.
    <nucleotides> is a pure string of nucleotides (no newlines),
    <bed> is a list of tuples, each item being (start, end) positions.'''
    mask_ranges = get_mask_ranges(bed_retain_ranges, nucleotides)
    masked_seq = ''
    for mask_range, retain_range in zip(mask_ranges, bed_retain_ranges):
        masked_seq += 'N' * (1 + mask_range[1] - mask_range[0])     # add masked nucleotides
        masked_seq += nucleotides[slice(retain_range[0]-1, retain_range[1])]
    masked_seq += 'N' * (1 + mask_ranges[-1][1] - mask_ranges[-1][0])
    return masked_seq


def main():
    '''main'''
    parser = argparse.ArgumentParser()
    if '--debug' not in sys.argv:
        parser.add_argument('fasta_filename',
                            type=str,
                            help='Filename for input fasta file. Can be gzipped.')
        parser.add_argument('bed_filename',
                            type=str,
                            help='Filename for masking the input fasta sequences')
        parser.add_argument('--out',
                            type=str,
                            help='Filename for writing output')

        args = parser.parse_args()
        fasta_filename = args.fasta_filename
        bed_filename = args.bed_filename
        out_filename = args.out

    start_time = time.time()
    # Initialize logging

    fastas = import_fasta(fasta_filename)
    bed = import_bed(bed_filename)

    for i in list(fastas):
        if i in list(bed):
            masked_sequence = wrap_fasta(mask_fasta(fastas[i], bed[i]))
            fastas[i] = masked_sequence

    with open(out_filename, 'w', encoding='utf-8') as outfile:
        for i in list(fastas):
            print('>' + i, file=outfile)
            print(fastas[i], file=outfile)

    run_time = (time.time() - start_time)
    logger.info("Total runtime: %s", run_time)
    sys.exit()


if __name__ == '__main__':
    main()
