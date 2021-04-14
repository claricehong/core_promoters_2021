#!/usr/bin/python3

import sys
import argparse
import os
import re
import gzip
from collections import defaultdict

parser = argparse.ArgumentParser(description = 'Gets insulator/TRIP barcode counts from fastq files')
parser.add_argument('fastq', help = 'input fastq file (either gzip or uncompressed)')
parser.add_argument('primer', help = 'primer sequence used for PCR')
parser.add_argument('TRIP_BCs', help = 'insulator names and barcodes')
parser.add_argument('output', help = 'output file basename')
args = parser.parse_args()

def hamdist(str1, str2): # From http://code.activestate.com/recipes/499304-hamming-distance/
    """Count the number of differences between equal length strings str1 and str2"""
    diffs = 0
    for ch1, ch2 in zip(str1, str2):
        if ch1 != ch2:
            diffs += 1
    return diffs

# some of the libraries were cloned slightly differently so the barcodes are in a slightly different position
# this function checks and returns the correct barcode
def check_random_bc(seq):
    """Finds barcode at correct position"""
    if seq.startswith('GATCA'):
        return seq[5:]
    else:
        return seq[:16]

# set output files basename
basename = os.path.basename(args.output)

iBCs = {}

# find the insulator BCs that we're using based on the insulators that are manually selected in the file without a # sign
with open(args.TRIP_BCs, 'r') as f: 
    header = f.readline()
    for line in f:
        line = line.strip().split('\t')
        if not line[0].startswith('#'):
            iBCs[line[1]] = line[0]

# generate regex to search sequence
# use multiple regexes to ensure that I get as many reads as possible
primer = '(' + args.primer + ')'
regex1 = re.compile(primer + r'([ATCGN]+$)')
seq_before_iBC = '(GCCCTT)'
regex2 = re.compile(seq_before_iBC + r'([ATCGN]+$)')

# regex to check whether the line is a sequence
seq_check = re.compile(r'([ATCGN]+$)')

# open fastq file
if str(args.fastq).endswith('.gz'):
    fastq_file = gzip.open(args.fastq, 'rt')
else:
    fastq_file = open(args.fastq, 'r')

# open file to write unmatched reads to
discarded = open(basename + '_discarded', 'w')

# initialise dictionary to store barcode counts
d_counts = defaultdict(int)
discard_count = 0
total_count = 0

# read and process fastq file
for line in fastq_file:
    total_count += 1
    line = line.strip()
    one = regex1.search(line)
    two = regex2.search(line)
    
    # find the trip (random) barcodes and insulator barcodes (iBC)
    # add barcodes to count dictionary if the iBC matches an insulator 
    if one != None:
        seq = one.group(2)
        iBC = seq[42:54]
        tBC = check_random_bc(seq[:21])
        if iBC in iBCs.keys():
            d_counts[(tBC, iBCs[iBC])] += 1
        else:
            discard_count += 1

    elif two != None:
        iBC = two.group(2)[:12]
        end_coordinate = 75 - 15 - len(two.group(0))
        seq = line[end_coordinate - 21 : end_coordinate]
        tBC = check_random_bc(seq)
        if iBC in iBCs.keys():
            d_counts[(tBC, iBCs[iBC])] += 1
        else:
            discard_count += 1

    # write the reads that didn't match to discard file
    else:
        n = seq_check.fullmatch(line)
        if n != None:
            discarded.write(line + '\n')
            discard_count += 1

# filter barcodes so that each barcode pair is supported by more than 5 reads
d_counts_filtered = {bc:count for (bc, count) in d_counts.items() if count > 5}
sorted_counts = sorted(d_counts_filtered.items(), key = lambda x: x[1], reverse = True)

# initialise dictionary to store collapsed barcodes counts
# barcodes are collapsed such that barcodes that are less than 3 hamming distance apart are considered the same barcode
# this analysis is similar to that described in the original TRIP paper (Akhtar et al.,2013, PMID: 23953119)
d_collapsed_counts = {}
multiple_match_count = 0

for barcodes, count in sorted_counts:
    tBC = barcodes[0]
    iBC = barcodes[1]
    matched_bcs = 0
    for existing_tBC, existing_iBC in d_collapsed_counts.keys():
        if iBC == existing_iBC and hamdist(tBC, existing_tBC) < 3:
            matched_bcs += 1
            
    if matched_bcs == 0:
        d_collapsed_counts[(tBC, iBC)] = count
    elif matched_bcs == 1:
        discard_count += count
    else:
        multiple_match_count += 1
        discard_count += count

print('number of reads discarded: ' + str(discard_count/(total_count/4)))
print('number of barcodes with multiple matches: ' + str(multiple_match_count))

# write barcodes and counts to file
with open(basename + '_counts', 'w') as f:
    f.write('tBC\tiBC\tcount\n')
    for key, value in d_collapsed_counts.items():
        f.write('{tBC}\t{iBC}\t{counts}\n'.format(tBC = key[0], iBC = key[1], counts = value))
