#!/usr/bin/python3

import sys
import argparse
import os
from collections import defaultdict

parser = argparse.ArgumentParser(description = 'Gets counts for each BC pair from demultiplexed output file')
parser.add_argument('readsfile', help = 'demultiplexed read file')
parser.add_argument('TRIP_BCs', help = 'insulator names and BCs file')
args = parser.parse_args()

basename = os.path.basename(args.readsfile)

d = defaultdict(int)

#generate lists of the BCs that we're using in this experiment
iBCs = {}

#find the gBCs that we're using based on the insulators that manually selected in the file with an asterisk
with open(args.TRIP_BCs, 'r') as f: 
    header = f.readline()
    for line in f:
        line = line.strip().split('\t')
        if not line[0].startswith('#'):
            iBCs[line[1]] = line[0]

#pull the tBCs out of the fastq file 
with open(args.readsfile) as f:
    for line in f:
        line = line.strip()
        if line.startswith('GATCA'):
            tBC = line[5:21]
            iBC = line[-12:]
            if iBC in iBCs.keys():
                d[(tBC, iBCs[iBC])] += 1
        else:
            tBC = line[:16]
            iBC = line[-12:]
            if iBC in iBCs.keys():
                d[(tBC, iBCs[iBC])] += 1


#write a list of the cBC, then the gBC, then the count of that particular combination
with open(basename + '_counts', 'w') as f:
    f.write('tBC\tiBC\tcount\n')
    for key, value in d.items():
        f.write('{tBC}\t{iBC}\t{counts}\n'.format(tBC = key[0], iBC = key[1], counts = value))
