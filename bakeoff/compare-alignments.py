import argparse
import gzip
import json
import sys

from ftx import FTX

parser = argparse.ArgumentParser(description='alignment evaluator')
parser.add_argument('files', nargs='+', help='run-aligner.py output files')
parser.add_argument('--tablename', required=True,
	help='base name for various tables')
parser.add_argument('--minexon', type=int, default=20,
	help='minimum exon length for table3 [%(default)i]')
parser.add_argument('--debug', action='store_true')
parser.add_argument('--verbose', action='store_true')
arg = parser.parse_args()

# open file pointers to all files
fps = []
progs = []
for file in arg.files:
	if   file.endswith('.ftx.gz'): fp = gzip.open(file, 'rt')
	elif file.endswith('.ftx'):    fp = open(file)
	else: sys.exit(f'unexpected file {file}')
	fps.append(fp)
	progs.append(file[:file.index('.')])

# aggregate alignment data by reference read
data = {}
while True:
	lines = []
	for fp in fps: lines.append(fp.readline().rstrip().split())
	if len(lines[0]) == 0: break
	for prog, (ref, ali) in zip(progs, lines):
		if ref not in data: data[ref] = {}
		data[ref][prog] = ali

###############
# Single Exon #
###############

d1 = {}
for rstr in data:
	ref = FTX.parse(rstr)
	if len(ref.exons) != 1: continue
	for prog, astr in data[rstr].items():
		if astr == 'None': r = 'unaligned'
		else:
			ali = FTX.parse(astr)
			if ref.beg == ali.beg and ref.end == ali.end: r = 'match'
			else: r = 'artefact'
		if prog not in d1: d1[prog] = {}
		if r not in d1[prog]: d1[prog][r] = 0
		d1[prog][r] += 1

categories = ('match', 'unaliged', 'artefact')

with open(f'{arg.tablename}_table1.tsv', 'w') as fp:
	print('program\t', '\t'.join(categories), file=fp)
	for prog in d1:
		print(prog, end='\t', file=fp)
		for cat in categories:
			if cat not in d1[prog]: print(0, end='\t', file=fp)
			else: print(d1[prog][cat], end='\t', file=fp)
		print(file=fp)

#############
# Two Exons #
#############

d2 = {}
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
				score = ref.compare_introns(ali)
				if   score == 2: r = 'match'
				elif score == 1: r = 'partial'
				else:            r = 'artefact'
		if emin not in d2: d2[emin] = {}
		if prog not in d2[emin]: d2[emin][prog] = {}
		if r not in d2[emin][prog]: d2[emin][prog][r] = 0
		d2[emin][prog][r] += 1

categories = ('match', 'partial', 'unspliced', 'unaligned', 'artefact')

with open(f'{arg.tablename}_table2.tsv', 'w') as fp:
	print('length\tprogram\t', '\t'.join(categories), file=fp)
	for length in sorted(d2):
		for prog in sorted(progs):
			print(length, prog, sep='\t', end='\t', file=fp)
			for r in categories:
				if r not in d2[length][prog]: print('0', end='\t', file=fp)
				else: print(d2[length][prog][r], end='\t', file=fp)
			print(file=fp)


###############
# Three Exons #
###############

d3 = {}
for rstr in data:
	ref = FTX.parse(rstr)
	if len(ref.exons) != 3: continue
	if ref.exon_length(0) < arg.minexon: continue
	if ref.exon_length(2) < arg.minexon: continue
	emin = ref.exon_length(1)

	for prog, astr in data[rstr].items():
		if astr == 'None': r = 'unaligned'
		else:
			ali = FTX.parse(astr)
			if len(ali.exons) == 1: r = 'unspliced'
			elif len(ali.exons) > 3: r = 'artefact'
			else:
				score = ref.compare_introns(ali)
				if   score == 4: r = 'match'
				elif score == 3: r = 'partial'
				elif score == 2: r = 'partial'
				elif score == 1: r = 'partial'
				else:            r = 'artefact'

		if emin not in d3: d3[emin] = {}
		if prog not in d3[emin]: d3[emin][prog] = {}
		if r not in d3[emin][prog]: d3[emin][prog][r] = 0
		d3[emin][prog][r] += 1

categories = ('match', 'partial', 'unspliced', 'unaligned', 'artefact')

with open(f'{arg.tablename}_table3.tsv', 'w') as fp:
	print('length\tprogram\t', '\t'.join(categories), file=fp)
	for length in sorted(d3):
		for prog in sorted(progs):
			print(length, prog, sep='\t', end='\t', file=fp)
			for r in categories:
				if r not in d3[length][prog]: print('0', end='\t', file=fp)
				else: print(d3[length][prog][r], end='\t', file=fp)
			print(file=fp)
