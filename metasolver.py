import argparse
from difflib import SequenceMatcher
import gzip
import json
import numpy as np
from scipy.cluster.hierarchy import linkage, dendrogram
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

parser = argparse.ArgumentParser(description='attempt the impossible')
parser.add_argument('json', help='metadata in json format')
parser.add_argument('--verbose', action='store_true')
arg = parser.parse_args()


if arg.json.endswith('.gz'): fp = gzip.open(arg.json, 'rt')
else:                        fp = open(arg.json)

# cluster_methods = ('complete', 'average', 'weighted', 'centroid', 'median', 'ward')

d = json.load(fp)

dcount = {}
for gse_id, gsmd in d.items():
	gsm_ids = list(gsmd.keys())
	if len(gsm_ids) < 4: continue
	condmat = []
	for i in range(len(gsm_ids)):
		for j in range(i+1, len(gsm_ids)):
			condmat.append(distance(gsmd[gsm_ids[i]], gsmd[gsm_ids[j]]))

	Z = linkage(condmat, method='weighted')
	dn = dendrogram(Z)
	group = {}
	for leaf, color in zip(dn['leaves'], dn['leaves_color_list']):
		if color not in group: group[color] = set()
		group[color].add(leaf)

	print(gse_id)
	for col in group:
		gsms = [gsm_ids[i] for i in group[col]]
		print('\t', gsms)

"""
	# make all the trees
	tree = {}
	for meth in cluster_methods:
		Z = linkage(condmat, method=meth)
		dn = dendrogram(Z)
		group = {}
		for leaf, color in zip(dn['leaves'], dn['leaves_color_list']):
			if color not in group: group[color] = set()
			group[color].add(leaf)
		tree[meth] = group

	# are all the trees the same?
	different = []
	same = []
	for i in range(len(cluster_methods)):
		a = tree[cluster_methods[i]]
		for j in range(len(cluster_methods)):
			b = tree[cluster_methods[j]]
			acol = set(a.keys())
			bcol = set(b.keys())
			if acol != bcol:
				different.append((cluster_methods[i], cluster_methods[j]))
				break
			for col in acol:
				if a[col] != b[col]:
					different.append((cluster_methods[i], cluster_methods[j]))
					break

	if different:
		for m1, m2 in different:
			if m1 not in dcount: dcount[m1] = {}
			if m2 not in dcount[m1]: dcount[m1][m2] = 0
			dcount[m1][m2] += 1

print(json.dumps(dcount, indent=4))

"""
