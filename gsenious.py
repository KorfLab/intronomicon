import argparse
import os
import sys

parser = argparse.ArgumentParser(description='experimental GSE extractor')
parser.add_argument('sheet', help='GSElog-sheet2.tsv')
parser.add_argument('dir', help='directory of GSEs')
parser.add_argument('--build', default='build/extracts',
	help='build directory')
arg = parser.parse_args()

os.system(f'mkdir -p {arg.build}')

data = {}
with open(arg.sheet) as fp:
	header = next(fp)
	for line in fp:
		f = line.split('\t')
		gse, gsm, file, sheet, strat, gene, count, note = f
		if count == '': continue
		if gse not in data: data[gse] = []
		data[gse].append({
			'gsm': gsm,
			'file': file,
			'strat': strat,
			'gene': gene,
			'count': count})

for gse, gsms in data.items():
	for gsm in gsms:
		sid = gsm['gsm']
		g = gsm['gene']
		c = gsm['count']
		infile = '/'.join((arg.dir, gse, gsm['file']))
		outfile = f'{arg.build}/{sid}.txt'
		print(infile, g, c, outfile, flush=True)
		os.system(f'./genecount-extractor {infile} {g} {c} | python3 name-resolver.py - worm_genes.json > {outfile}')
#		sys.exit()
