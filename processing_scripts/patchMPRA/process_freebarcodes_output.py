#!/usr/bin/python

import sys
import argparse
import os
from collections import defaultdict

parser = argparse.ArgumentParser(
    description = 'Process freebarcodes output and output counts per barcode pair')
parser.add_argument(
    'freebarcodes_output', help = 'decoded freebarcodes file')
parser.add_argument(
    '-o', '--output', help = 'output file basename')
args = parser.parse_args()

basename = os.path.basename(args.freebarcodes_output).split('.')[0]

d_counts = defaultdict(int)

with open(args.freebarcodes_output, 'r') as f:
    for line in f:
        line = line.strip('\n').split('\t')
        gBC = line[0].split(':')[-1]
        pBC = line[1]
        d_counts[(pBC, gBC)] += 1

with open(basename + '_counts.txt', 'w') as f:
    f.write('pBC\tgBC\tcount\n')
    for bcs, current_count in d_counts.items():
        f.write('{pBC}\t{gBC}\t{count}\n'.format(pBC = bcs[0], gBC = bcs[1], count = current_count))