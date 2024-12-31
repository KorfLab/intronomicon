import argparse
import gzip
import json

from ftx import FTX

parser = argparse.ArgumentParser(description='alignment evaluator')
parser.add_argument('files', nargs='+', help='run-aligner.py output files')
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

# iterate through files in parallel
data = {}
while True:
	lines = []
	for fp in fps: lines.append(fp.readline().rstrip().split())
	if len(lines[0]) == 0: break
	for prog, (ref, ali) in zip(progs, lines):
		if ref not in data: data[ref] = {}
		data[ref][prog] = ali

# gather statistics - incomplete
for rstr in data:
	ref = FTX.parse(rstr)
	print(ref)
	for prog, astr in data[rstr].items():
		if astr == 'None':
			pass
		else:
			ali = FTX.parse(astr)
			ali.info = prog
			print(ali)
	print()
