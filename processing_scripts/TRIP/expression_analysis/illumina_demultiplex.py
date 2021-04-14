#!/usr/bin/python

import sys
import re
import argparse
import os
from collections import defaultdict
import gzip

parser = argparse.ArgumentParser(description = 'Demultiplexes Illumina sequencing fastq files')
parser.add_argument('fastq', help = 'input fastq file (gzip or not gzip is fine)')
parser.add_argument('MPBC', help = 'file containing multiplexing barcodes')
parser.add_argument('middle', help = 'middle sequence between multiplexing barcode and labelling barcode')
parser.add_argument('BClength', type = int, help = 'length of labelling BC')
parser.add_argument('-o', '--output', help = 'output file basename')
args = parser.parse_args()
			
# set output files basename
# defaults to empty string if no input given
if args.output is not None:
	basename = args.output + '_'
else:
	basename = os.path.basename(args.fastq)

# get multiplexing barcodes from file as dictionary
d_MPBC_BCname = {}
MPBCs = []
outfiles = []
lookup = open(basename + "lookup.txt", "w")

with open(args.MPBC, 'r') as f:
	header = f.readline()
	for idx, line in enumerate(f):
		line = line.strip().split('\t')
		MPBCs.append(line[1])
		d_MPBC_BCname[line[1]] = idx
		outfiles.append(open(basename + line[0], 'w'))
		lookup.write(basename + line[0] + '\n')

discarded = open(basename + 'discarded', 'w')

#generate regex to search sequence
MPBCs = '(' + '|'.join(MPBCs) + ')'
middle = '(' + args.middle + ')'
p = re.compile(MPBCs + middle + r'([ATCGN]+$)')	
q = re.compile(r'([ATCGN]+$)')

#checks if the line is a sequence file 
#checks if it matches the pattern defined by the regex above 
#writes BC to appropriate file if it is

if str(args.fastq).endswith('.gz'):
	with gzip.open(args.fastq, 'rt') as f:
		for line in f:
			line = line.strip()
			m = p.search(line)
			if m != None:
				current_mpbc_name = d_MPBC_BCname[m.group(1)]
				bc = m.group(3)[:args.BClength]
				outfiles[current_mpbc_name].write(bc + '\n')
			else:
				n = q.fullmatch(line)
				if n != None:
					discarded.write(line + '\n')
else:
	with open(args.fastq, 'r') as f:
		for line in f:
			line = line.strip()
			m = p.search(line)
			if m != None:
				current_mpbc_name = d_MPBC_BCname[m.group(1)]
				bc = m.group(3)[:args.BClength]
				outfiles[current_mpbc_name].write(bc + '\n')

for i in outfiles:
	i.close()

discarded.close()

							
