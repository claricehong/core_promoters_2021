#!/usr/bin/python

import sys
import argparse
import pysam
from collections import Counter
from collections import defaultdict

parser = argparse.ArgumentParser(
    description = 'Read bam file and output barcode with annotated locations')
parser.add_argument(
    'bam', help = 'bamfile')
parser.add_argument(
    'insulator_bcs', help = 'insulator BC file')
parser.add_argument(
    'output', help = 'output file name')
args = parser.parse_args()

# create barcodes with one substitution because the mapping is full of errors
def generate_one_substitution(seq):

    out = [seq]

    for position, nt in enumerate(seq):
        nucleotides = ['A', 'T', 'C', 'G']
        nucleotides.remove(nt)
        for other_nt in nucleotides:
            new = seq[:position] + other_nt + seq[position+1:]
            out.append(new)

    return out

def rev_comp(seq):
    
    new_seq = ''
    conversion = {'A': 'T', 'T': 'A', 'G':'C', 'C':'G'}
    
    for i in reversed(seq):
        new_seq += conversion[i]
        
    return new_seq

def check_Ns(seq):

    valid = ['A', 'T', 'C', 'G']

    for i in seq:
        if i not in valid:
            return False
    
    return True

def get_strand(read):
    
    if read.is_reverse:
        return '-'
    
    return '+'


d_insulator_bcs = {}

with open(args.insulator_bcs, 'r') as f:
    header = f.readline()
    for line in f:
        name, bc = line.strip('\n').split('\t')
        related_bcs = generate_one_substitution(bc)
        for each_bc in related_bcs:
            d_insulator_bcs[each_bc] = name

bamfile = pysam.AlignmentFile(args.bam, "rb")

count_mapped_barcodes = 0
count_matched_barcodes = 0
all_locations = defaultdict(list)

for read in bamfile:
    if not read.is_unmapped:
        if not read.is_secondary:
            if read.is_read1:
                count_mapped_barcodes += 1
                insulator_bc = read.query_name.split(':')[-2]
                trip_bc = read.query_name.split(':')[-1]
                chromosome = read.reference_name
                location = str(read.get_reference_positions()[0])
                if insulator_bc in d_insulator_bcs.keys() and check_Ns(trip_bc) == True:
                    count_matched_barcodes += 1
                    insulator = d_insulator_bcs[insulator_bc]
                    strand = get_strand(read)
                    all_locations[(insulator, rev_comp(trip_bc))].append((chromosome, location, strand))

all_locations = {bcs: locations for bcs, locations in all_locations.items() if len(locations) > 2}

print('Number of mapped barcodes: ' + str(count_mapped_barcodes))
print('Number of matched barcodes: ' + str(count_matched_barcodes))

barcode_locations = {}

for bcs, locations in all_locations.items():
    count_items = Counter(locations)
    top = count_items.most_common(2)
    total = sum(count_items.values())
    if top[0][1]/total > 0.9:
        barcode_locations[bcs] = [top[0][0], top[0][1]]
    elif (top[0][1] + top[1][1])/total > 0.9:
        loc1 = top[0][0]
        loc2 = top[1][0]
        if loc1[0] == loc2[0]:
            if abs(int(loc1[1]) - int(loc2[1])) < 1000:
                barcode_locations[bcs] = [top[0][0], top[0][1] + top[1][1]]

all_tBCs = Counter([i[1] for i in barcode_locations.keys()])
repeat_tBCs = [k for k, v in all_tBCs.items() if v >= 2]
repeat_barcode_locations = defaultdict(list)

for bc, location in barcode_locations.items():
    if bc[1] in repeat_tBCs:
        repeat_barcode_locations[bc[1]].append((bc[0], location[0], location[1]))

keys_to_drop = []

for tBC, info in repeat_barcode_locations.items():
    counts = [i[2] for i in info]
    locations = [i[1] for i in info]
    iBCs = [i[0] for i in info]
    if len(set(locations)) == 1:
        if counts[0]/sum(counts) > 0.9:
            for iBC in iBCs[1:]:
                keys_to_drop.append((iBC, tBC))

for key in keys_to_drop:
    if key in barcode_locations.keys():
        del barcode_locations[key]

with open(args.output, 'w') as f:
    for bcs, location in barcode_locations.items():
        f.write('chr{chrom}\t{loc}\t{loc}\t{iBC}\t{tBC}\t{count}\t{strand}\n'.format(
            chrom = location[0][0], loc = location[0][1], iBC = bcs[0], tBC = bcs[1], count = location[1], strand = location[0][2]))

print('Number of unique barcodes :' + str(len(all_locations)))
print('Number of barcodes that match only to one location: ' + str(len(barcode_locations)))

bamfile.close()
