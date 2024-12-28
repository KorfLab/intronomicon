import argparse
import gzip
import sys

import bakeoff

parser = argparse.ArgumentParser('convert sam to pairwise alignment')
parser.add_argument('ref', help='reference genome')
parser.add_argument('sam', help='sam file')
arg = parser.parse_args()


chrom = {}
for name, seq in bakeoff.readfasta(arg.ref):
	chrom[name] = seq

for query, seq, pos, cigar, bflag in bakeoff.readsam(arg.sam):
	print(cigar)
