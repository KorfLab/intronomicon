import argparse
import gzip
import sys

from ftx import FTX

parser = argparse.ArgumentParser(description='alignment evaluator')
parser.add_argument('files', nargs='+', help='run-aligner.py output files')
parser.add_argument('--table2', required=True, metavar='<tsv>',
	help='2 exon table')
parser.add_argument('--debug', action='store_true')
parser.add_argument('--verbose', action='store_true')
arg = parser.parse_args()
if not arg.table2.endswith('tsv'): sys.exit('hey dummy, use .tsv extension')

# open file pointers to all files
fps = []
progs = []
for file in arg.files:
	if   file.endswith('.ftx.gz'): fp = gzip.open(file, 'rt')
	elif file.endswith('.ftx'):    fp = open(file)
	else: sys.exit(f'unexpected file {file}')
	fps.append(fp)
	progs.append(file[:file.index('.')])

# aggregate alignment data
data = {}
while True:
	lines = []
	for fp in fps: lines.append(fp.readline().rstrip().split())
	if len(lines[0]) == 0: break
	for prog, (ref, ali) in zip(progs, lines):
		if ref not in data: data[ref] = {}
		data[ref][prog] = ali

# gather statistics - incomplete

###############
# Single Exon #
###############


#############
# Two Exons #
#############
oh = {}
for rstr in data:
	ref = FTX.parse(rstr)
	if len(ref.exons) != 2: continue
	emin = min(ref.exon_length(0), ref.exon_length(1))
	for prog, astr in data[rstr].items():
		if astr == 'None': r = 'unaligned'
		else:
			ali = FTX.parse(astr)
			if len(ali.exons) == 1: r = 'unspliced'
			elif len(ali.exons) > 2: r = 'artefact'
			else:
				found = 0
				if ref.exons[0][1] == ali.exons[0][1]: found += 1
				if ref.exons[1][0] == ali.exons[1][0]: found += 1

				if   found == 2: r = 'match'
				elif found == 1: r = 'partial'
				else:            r = 'artefact'
		if emin not in oh: oh[emin] = {}
		if prog not in oh[emin]: oh[emin][prog] = {}
		if r not in oh[emin][prog]: oh[emin][prog][r] = 0
		oh[emin][prog][r] += 1

categories = ('match', 'partial', 'unspliced', 'unaligned', 'artefact')

with open(arg.table2, 'w') as fp:
	print('length\tprogram\t', '\t'.join(categories), file=fp)
	for length in sorted(oh):
		for prog in sorted(progs):
			print(length, prog, sep='\t', end='\t', file=fp)
			for r in categories:
				if r not in oh[length][prog]: print('0', end='\t', file=fp)
				else: print(oh[length][prog][r], end='\t', file=fp)
			print(file=fp)


###############
# Three Exons #
###############

