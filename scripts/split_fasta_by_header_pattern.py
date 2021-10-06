#!/usr/bin/env python

import argparse
import re
import sys
from Bio import SeqIO

parser = argparse.ArgumentParser(description="Search for pattern in fasta header and return sequences whose header includes (or does not include) the match")

parser.add_argument('--fasta', required=True, help='input fasta file', action='store')
parser.add_argument('--pattern', required=True, help='pattern to match in header', action='store')
parser.add_argument('--complement', required=False, default=False, help='Select only sequences that do not match the pattern', action='store_true')

args=parser.parse_args()

# read in sequences
records = list(SeqIO.parse(args.fasta, "fasta"))

# create empty lists
found_sequences = []
missing_sequences = []

# If --complement is provided, return sequences that do not match the search string. Otherwise return matches.

for record in records:

        if re.search(args.pattern, record.description):
                found_sequences.append(record)
        else:
                missing_sequences.append(record)

# If search pattern is not found in any header, print error message only
# If search pattern is found, print sequences that match (or do not match, if --complement argument is provided) to screen

if args.complement:
        if len(missing_sequences) == len(records):
                print("ERROR: Could not find " + args.pattern + " in any fasta header!", file=sys.stderr)
        else:
                print("Writing out %i sequences that DO NOT match the search string!" % len(missing_sequences), file=sys.stderr)
                SeqIO.write(missing_sequences, sys.stdout, "fasta")
else:
        print("Writing out %i sequences that match the search string!" % len(found_sequences), file=sys.stderr)
        SeqIO.write(found_sequences, sys.stdout, "fasta")

