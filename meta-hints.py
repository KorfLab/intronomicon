import argparse
from difflib import SequenceMatcher
import gzip
import json
import matplotlib.pyplot as plt
import numpy as np
import os
from scipy.cluster.hierarchy import linkage, dendrogram
import sqlite3
import sys

def distance(a, b):
	if a == b: return 0

	tags = set()
	for tag in a: tags.add(tag)
	for tag in b: tags.add(tag)

	d = 0
	for tag in tags:
		if tag not in a or tag not in b:
			d += 1
			continue
		x = 1 - SequenceMatcher(None, a[tag], b[tag]).ratio()
		d += x

	return d

def set_tagval(d, s):
	tag, val = s.split(':', maxsplit=1)
	tag = tag.rstrip().lstrip()
	val = val.rstrip().lstrip()
	if tag in d: d[tag] += '; {val}'
	d[tag] = val

def simplify(text):
	fields = ('Sample_title', 'Sample_source_name_ch1', 'tissue',
		'developmental state', 'genotype')
	simple = [text[field] for field in fields if field in text]
	return ' | '.join(simple)


cluster_methods = ('single', 'complete', 'average', 'weighted', 'centroid',
	'median', 'ward')

parser = argparse.ArgumentParser(description='attempt the impossible')
parser.add_argument('db', help='database')
parser.add_argument('out', help='output directory')
parser.add_argument('--method', default='weighted',
	help=f'clustering method [%(default)s] {cluster_methods}')
arg = parser.parse_args()

os.system(f'mkdir -p {arg.out}/png {arg.out}/r1/original')

con = sqlite3.connect(arg.db)
cur = con.cursor()

# gsm stuff
cur.execute(f'SELECT gse_id, gsm_id, gsm_txt FROM run INNER JOIN experiment ON experiment.srx_id = run.srx_id WHERE platform = "ILLUMINA" and rlen >= 100')
d = {}
for gse_id, gsm_id, gsm_txt in cur.fetchall():
	sd = {}
	for thing in gsm_txt.split(' #+# '):
		if thing.startswith('Sample_character'):
			for text in thing[28:].split('#-#'): set_tagval(sd, text)
		else:
			set_tagval(sd, thing)
	if gse_id not in d: d[gse_id] = {}
	d[gse_id][gsm_id] = sd

# gse text
gtxt = {}
cur.execute(f'SELECT gse_id, gse_txt FROM series')
for gse_id, gse_txt in cur.fetchall():
	gtxt[gse_id] = gse_txt

# something
for gse_id, gsmd in d.items():
	gsm_ids = list(gsmd.keys())
	if len(gsm_ids) < 4: continue # skipping small series for now

	# create pairwise distances among GSMs
	condmat = []
	for i in range(len(gsm_ids)):
		for j in range(i+1, len(gsm_ids)):
			condmat.append(distance(gsmd[gsm_ids[i]], gsmd[gsm_ids[j]]))

	# create dendrogram
	plt.figure(figsize=(12,8))
	plt.title(gse_id)
	Z = linkage(condmat, method='weighted')
	dn = dendrogram(Z, labels=gsm_ids, orientation='left')
	plt.savefig(f'{arg.out}/png/{gse_id}.png')
	plt.close()

	# create grouped GSE ids
	group = {}
	for leaf, color in zip(dn['leaves'], dn['leaves_color_list']):
		if color not in group: group[color] = set()
		group[color].add(leaf)

	# create page
	with open(f'{arg.out}/src/{gse_id}.html', 'w') as fp:
		print(f'<html><head><title>{gse_id}</title></head>\n<body>', file=fp)
		print(f'<h1>{gse_id}</h1>', file=fp)
		print(f'<img src="png/{gse_id}.png" style="float:right;width:600px;height:400px;">', file=fp)
		section = gtxt[gse_id].split('#+#')
		for txt in section:
			txt = txt.rstrip().lstrip()
			if txt.startswith('Series_sample_id:'): continue
			print(f'<li>{txt}', file=fp)
		for color in group:
			print(f'<dl><dt>{color}</dt>', file=fp)
			gsms = [gsm_ids[i] for i in group[color]]
			for gsm_id in gsms:
				print(f'<dd><b>{gsm_id}</b> {simplify(d[gse_id][gsm_id])}</dd>', file=fp)
			print('</dl>', file=fp)
		print('</body>\n</html>', file=fp)
