#!/usr/bin/python

import sys
import argparse
import pysam
import collections
import itertools

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
        name, bc = line.strip().split('\t')
        related_bcs = generate_one_substitution(bc)
        for each_bc in related_bcs:
            d_insulator_bcs[each_bc] = name

bamfile = pysam.AlignmentFile(args.bam, "rb")

count_mapped_barcodes = 0
count_matched_barcodes = 0
all_locations = collections.defaultdict(list)

for read in bamfile:
    if not read.is_unmapped:
        if not read.is_secondary:
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
barcode_locations2 = {}

def check_percentage(location, total, top_count):

    if location[1]/total > 0.1 and location[1]*5 > top_count:
        return 'over'
    else:
        return 'below'

def check_all_percentages(all_locations, total, top_count):

    over_results = []

    for location in all_locations:
        if check_percentage(location, total, top_count) == 'over':
            over_results.append(location[0])

    if len(over_results) == 0:
        return 'good'
    else:
        return over_results

def compare_locations(loc1, loc2):

    if loc1[0] == loc2[0] and loc1[2] == loc2[2] and abs(int(loc1[1]) - int(loc2[1])) < 1000:
        return 'close'
    else:
        return 'far'


for bcs, locations in all_locations.items():
    count_items = collections.Counter(locations)
    top = count_items.most_common()
    total = sum(count_items.values())

    collated_locs = {}
    collated_locs[top[0][0]] = top[0][1]

    for loc, count in top:
        added = 'no'
        existing_collated_locs = [k for k, v in collated_locs.items()]
        for existing_loc in existing_collated_locs:
            if compare_locations(existing_loc, loc) == 'close':
                collated_locs[existing_loc] += count
                added = 'yes'
        if added == 'no':
            # elif compare_locations(existing_loc, loc) == 'far':
            collated_locs[loc] = count

    collated_locs_sorted = sorted(collated_locs.items(), key = lambda x: x[1], reverse = True)
    top_count = collated_locs_sorted[0][1]

    if top_count/total > 0.7:
        other_locs = collated_locs_sorted[1:]
        if check_all_percentages(other_locs, total, top_count) == 'good':
            barcode_locations[bcs] = collated_locs_sorted[0]
    #     else:
    #         print(collated_locs_sorted)
    # else:
    #     print(collated_locs_sorted)
    # elif total > 100:
    #     if top[0][0][0] == top[1][0][0]:
    #         print(top)
    #         print(collated_locs_sorted)

all_tBCs = collections.Counter([i[1] for i in barcode_locations.keys()])
repeat_tBCs = [k for k, v in all_tBCs.items() if v >= 2]
repeat_barcode_locations = collections.defaultdict(list)

for bc, location in barcode_locations.items():
    if bc[1] in repeat_tBCs:
        repeat_barcode_locations[bc[1]].append((bc[0], location[0], location[1]))

keys_to_drop = []

for tBC, info in repeat_barcode_locations.items():
    counts = [i[2] for i in info]
    locations = [i[1] for i in info]
    iBCs = [i[0] for i in info]
    if len(set(locations)) == 1:
        if counts[0]/sum(counts) > 0.8:
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
