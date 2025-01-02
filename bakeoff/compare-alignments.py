import argparse
import gzip
import json
import sys

from ftx import FTX

parser = argparse.ArgumentParser(description='alignment evaluator')
parser.add_argument('files', nargs='+', help='run-aligner.py output files')
parser.add_argument('--basename', required=True,
	help='base name for various output files')
parser.add_argument('--minexon', type=int, default=20,
	help='minimum exon length for table3 [%(default)i]')
parser.add_argument('--experimental', action='store_true',
	help='special processing for the synthetic genome experiment')
parser.add_argument('--debug', action='store_true')
parser.add_argument('--verbose', action='store_true')
arg = parser.parse_args()

# Aggregate alignment data by reference read

fps = []
progs = []
for file in arg.files:
	if   file.endswith('.ftx.gz'): fp = gzip.open(file, 'rt')
	elif file.endswith('.ftx'):    fp = open(file)
	else: sys.exit(f'unexpected file {file}')
	fps.append(fp)
	progs.append(file[:file.index('.')])

data = {}
while True:
	lines = []
	for fp in fps: lines.append(fp.readline().rstrip().split())
	if len(lines[0]) == 0: break
	for prog, (ref, ali) in zip(progs, lines):
		if ref not in data: data[ref] = {}
		data[ref][prog] = ali

dfp = open(f'{arg.basename}_details.txt', 'w') # for detailed output

# Single exon

d1 = {}
for rstr in data:
	ref = FTX.parse(rstr)
	if len(ref.exons) != 1: continue
	print(ref, file=dfp)
	for prog, astr in data[rstr].items():
		if astr == 'None':
			r = 'unaligned'
			ali = 'None'
		else:
			ali = FTX.parse(astr)
			if len(ali.exons) != 1:
				r = 'artefact'
			elif ref.beg <= ali.beg and ref.end >= ali.beg:
				if ref.beg == ali.beg and ref.end == ali.end: r = 'hit'
				else: r = 'near'
			else: r = 'miss'
		print('\t', ali, r, file=dfp)
		if prog not in d1: d1[prog] = {}
		if r not in d1[prog]: d1[prog][r] = 0
		d1[prog][r] += 1

categories = ('hit', 'miss', 'near', 'artefact', 'unaligned')

with open(f'{arg.basename}_table1.tsv', 'w') as fp:
	print('program\t', '\t'.join(categories), file=fp)
	for prog in d1:
		print(prog, end='\t', file=fp)
		for cat in categories:
			if cat not in d1[prog]: print(0, end='\t', file=fp)
			else: print(d1[prog][cat], end='\t', file=fp)
		print(file=fp)

# Two exons
# Experimental: segment by site

d2 = {}
for rstr in data:
	ref = FTX.parse(rstr)
	if len(ref.exons) != 2: continue
	print(ref, file=dfp)
	emin = min(ref.exon_length(0), ref.exon_length(1))
	for prog, astr in data[rstr].items():
		if astr == 'None':
			r = 'unaligned'
			ali = 'None'
		else:
			ali = FTX.parse(astr)
			if len(ali.exons) == 1: r = 'under'
			elif len(ali.exons) > 2: r = 'over'
			else:
				score = ref.compare_introns(ali)
				if   score == 2: r = 'correct'
				elif score == 1: r = 'partial'
				else:            r = 'wrong'
		print('\t', ali, r, file=dfp)
		if emin not in d2: d2[emin] = {}
		if prog not in d2[emin]: d2[emin][prog] = {}
		if r not in d2[emin][prog]: d2[emin][prog][r] = 0
		d2[emin][prog][r] += 1

categories = ('correct', 'partial', 'over', 'under', 'unaligned', 'wrong')

with open(f'{arg.basename}_table2.tsv', 'w') as fp:
	print('length\tprogram\t', '\t'.join(categories), file=fp)
	for length in sorted(d2):
		for prog in sorted(progs):
			print(length, prog, sep='\t', end='\t', file=fp)
			for r in categories:
				if r not in d2[length][prog]: print('0', end='\t', file=fp)
				else: print(d2[length][prog][r], end='\t', file=fp)
			print(file=fp)

# Three Exons
# Experimental: focus on middle exon

d3 = {}
for rstr in data:
	ref = FTX.parse(rstr)
	if len(ref.exons) != 3: continue
	if ref.exon_length(0) < arg.minexon: continue
	if ref.exon_length(2) < arg.minexon: continue

	print(ref, file=dfp)
	emin = ref.exon_length(1)

	for prog, astr in data[rstr].items():
		if astr == 'None':
			r = 'unaligned'
			ali = 'None'
		else:
			ali = FTX.parse(astr)
			if len(ali.exons) == 1: r = 'unspliced'
			else:
				score = ref.compare_introns(ali)
				if   score == 4: r = 'correct'
				elif score == 3: r = 'found3'
				elif score == 2: r = 'found2'
				elif score == 1: r = 'found1'
				else:            r = 'artefact'
		print('\t', ali, r, file=dfp)
		if emin not in d3: d3[emin] = {}
		if prog not in d3[emin]: d3[emin][prog] = {}
		if r not in d3[emin][prog]: d3[emin][prog][r] = 0
		d3[emin][prog][r] += 1

categories = ('correct', 'found1', 'found2', 'found3', 'unaligned', 'artefact')

with open(f'{arg.basename}_table3.tsv', 'w') as fp:
	print('length\tprogram\t', '\t'.join(categories), file=fp)
	for length in sorted(d3):
		for prog in sorted(progs):
			print(length, prog, sep='\t', end='\t', file=fp)
			for r in categories:
				if r not in d3[length][prog]: print('0', end='\t', file=fp)
				else: print(d3[length][prog][r], end='\t', file=fp)
			print(file=fp)

# Four Exons
# This doesn't happen in the synthetic data experiments
# In real data, we are ignoring this for now

dfp.close()
