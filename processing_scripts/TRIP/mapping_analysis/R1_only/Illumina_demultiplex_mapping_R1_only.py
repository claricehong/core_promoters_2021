#!/usr/bin/python

import sys
import re
import argparse
import os
import gzip

parser = argparse.ArgumentParser(
    description = 'Demultiplexes Illumina sequencing fastq files')
parser.add_argument(
    'mpBC', help = 'file containing multiplexing barcodes')
parser.add_argument(
    'fastq_R1', help = 'input fastq file R1')
parser.add_argument(
    'fastq_R2', help = 'input fastq file R2')
parser.add_argument(
    '-o', '--output', help = 'output file basename')
args = parser.parse_args()

def process_fastq(fastq_file):

    current_record = {}

    for name, seq, crap, quality in zip(*[iter(fastq_file)]*4):
        current_record['name'] = name.strip('\n')
        current_record['seq'] = seq.strip('\n')
        current_record['quality'] = quality.strip('\n')

        yield current_record

def get_fastq_id(fastq_name):

    return fastq_name.split(' ')[0]

# set output files basename
# defaults to empty string if no input given
if args.output is not None:
    basename = args.output + '_'
else:
    basename = ''

lookup = open(basename + 'lookup.txt', 'w') 
R1_discarded = open(basename + 'discarded_R1', 'w')
R2_discarded = open(basename + 'discarded_R2', 'w')

R1_bc_file = {}
all_bcs = []
d_matching_bcs = {}

with open(args.mpBC, 'r') as bc_file:
    header = bc_file.readline()
    for entry in bc_file:
        entry = entry.strip().split('\t')
        bc_name = entry[0]
        bcs = entry[1:]
        all_bcs.extend(bcs)
        main_bc = entry[1]
        for each_bc in bcs:
            d_matching_bcs[each_bc] = main_bc

        current_base = basename + bc_name 
        R1_bc_file[main_bc] = open(current_base + '_R1', 'w')
        lookup.write('{current}_R1\t{current}\n'.
            format(current=current_base))

MPBCs = '(' + '|'.join(all_bcs) + ')'
mpBC_re = re.compile(MPBCs + r'([ATCGN]+$)')

# open fastq files
if str(args.fastq_R1).endswith('.gz'):
    R1_open_file = gzip.open(args.fastq_R1, 'rt')
else:
    R1_open_file = open(args.fastq_R1, 'r')

if str(args.fastq_R2).endswith('.gz'):
    R2_open_file = gzip.open(args.fastq_R2, 'rt')
else:
    R2_open_file = open(args.fastq_R2, 'r')

# read fastq files
R1_reader = process_fastq(R1_open_file)
R2_reader = process_fastq(R2_open_file)

def check_random_bc(seq):
    if seq.startswith('TGATC'):
        return seq[5:]
    else:
        return seq[:16]

for R1_record, R2_record in zip(R1_reader, R2_reader):
    
    if R1_record['name'].split(' ')[0] != R2_record['name'].split(' ')[0]:
        print(R1_record['name'])
        print(R2_record['name'])
        break

    bc_match = mpBC_re.search(R1_record['seq'])
    if bc_match != None:

        current_mpbc = d_matching_bcs[bc_match.group(1)]
        if R2_record['seq'].startswith('AAC'):
            insulator_bc = R2_record['seq'][20:32]
            if insulator_bc == 'GGTTAGCGGTAA':
                random_bc = R2_record['seq'][53:69]
            else:
                random_bc = check_random_bc(R2_record['seq'][53:74])
        else:
            insulator_bc = R2_record['seq'][17:29]
            if insulator_bc == 'GGTTAGCGGTAA':
                random_bc = R2_record['seq'][53:69]
            else:
                random_bc = check_random_bc(R2_record['seq'][50:71])

        R1_name = R1_record['name'].split(' ')[0] + ':' + insulator_bc + ':' + random_bc
        R1_seq = bc_match.group(2)[31:]
        R1_quality = R1_record['quality'][-len(R1_seq):]

        if len(R1_seq) > 30:
            R1_bc_file[current_mpbc].write('{name}\n{seq}\n+\n{quality}\n'.format(
                name = R1_name, seq = R1_seq, quality = R1_quality))
    
    else:
        R1_discarded.write(R1_record['seq'] + '\n')

for open_file in R1_bc_file.values():
    open_file.close()

lookup.close()
R1_discarded.close()
