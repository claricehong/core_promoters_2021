#!/usr/bin/python

import sys
import argparse
import os
import gzip
import random

parser = argparse.ArgumentParser(
    description = 'Change the barcodes with N in them to a random ATCG because the free barcodes software does not accept them and find LP barcode and adds LP name to read name')
parser.add_argument(
    'fastq_R1', help = 'fastq file read 1')
parser.add_argument(
    'fastq_R2', help = 'fastq file read 2')
parser.add_argument(
    'gBC', help = 'file containing genomic barcodes')
parser.add_argument(
    '-o', '--output', help = 'output file basename')
args = parser.parse_args()

def process_fastq(fastq_file):

    current_record = {}

    for name, seq, crap, quality in zip(*[iter(fastq_file)]*4):
        current_record['name'] = name.strip('\n')
        current_record['seq'] = seq.strip('\n')
        current_record['crap'] = crap.strip('\n')
        current_record['quality'] = quality.strip('\n')

        yield current_record

def rev_comp(seq):

    final = ''
    d = {'A':'T', 'T':'A', 'C':'G', 'G':'C', 'N':'N'}

    for base in reversed(seq):
        final += d[base]

    return final

# set output files basename
# defaults to empty string is basename not given
if args.output is not None:
    basename = args.output
else:
    basename = ''

# read LP barcodes
gBCs = {}

with open(args.gBC, 'r') as f:
    header = f.readline()
    for line in f:
        if not line.startswith('*'):
            line = line.strip('\n').split('\t')
            gBCs[line[1]] = line[0]

if str(args.fastq_R1).endswith('.gz'):
    R1_open_file = gzip.open(args.fastq_R1, 'rt')
else:
    R1_open_file = open(args.fastq_R1, 'r')

if str(args.fastq_R2).endswith('.gz'):
    R2_open_file = gzip.open(args.fastq_R2, 'rt')
else:
    R2_open_file = open(args.fastq_R2, 'r')

R1_reader = process_fastq(R1_open_file)
R2_reader = process_fastq(R2_open_file)
output_file = open(basename + '.fastq', 'w')
discarded = open(basename + '_discarded', 'w')

count = 0
total_count = 0

random.seed(23948)

for read1, read2 in zip(R1_reader, R2_reader):
    total_count += 1
    check_R2_re = read2['seq'][:45].split('GATCC')
    check_R1_re = read1['seq'][:44].split('CTAGA')

    if len(check_R2_re) == 2 and len(check_R1_re) == 2:
        if check_R1_re[1] != '':
            current_gBC = rev_comp(check_R2_re[1][:12])
            current_cBC = check_R1_re[1][:12]
            current_name = read1['name'].split(' ')
            if current_gBC in gBCs.keys():
                current_gBC_name = gBCs[current_gBC]
                current_name_with_gBC = current_name[0] + ':LP' + current_gBC_name + ' ' + current_name[1]
                if 'N' in current_cBC:
                    chosen_base = random.choice(['A', 'T', 'C', 'G'])
                    new_seq = check_R1_re[1].replace('N', chosen_base)
                    output_file.write('{name}\n{seq}\n{crap}\n{quality}\n'.format(
                        name = current_name_with_gBC, seq = new_seq, crap = read1['crap'][:len(check_R1_re[1])], quality = read1['quality'][:len(check_R1_re[1])]))
                else:
                    output_file.write('{name}\n{seq}\n{crap}\n{quality}\n'.format(
                        name = current_name_with_gBC, seq = check_R1_re[1], crap = read1['crap'][:len(check_R1_re[1])], quality = read1['quality'][:len(check_R1_re[1])]))
            else:
                count += 1
                discarded.write(read1['seq'] + '\n')
    else:
        count += 1
        discarded.write(read1['seq'] + '\n')

print('proportion of abandoned reads: ' + '{:.2f}'.format(count/total_count))
