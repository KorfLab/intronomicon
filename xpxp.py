import json
import math
import sys
import matplotlib.pyplot as plt
import numpy as np
import os
from scipy.cluster.hierarchy import linkage, dendrogram

def dtc(P, Q):
	d = 0
	for k in P: d += abs(P[k] - Q[k])
	return d


def count_type(d):
	ints = 0
	floats = 0
	for gene, val in d.items():
		if math.isclose(0, val % 1): ints += 1
		else: floats += 1
	if   ints == len(d): return int
	elif floats > 0.99 * len(d): return float
	else: sys.exit('what to do')


def read_expression(file):
	t = {}
	with open(file) as fp:
		for line in fp:
			gene, val = line.split()
			t[gene] = float(val)
	return t

# get gene lengths (hard-coded)
gene_len = {}
with open('worm_genes.txt') as fp:
	header = next(fp)
	for line in fp:
		gene, size = line.split()
		gene_len[gene] = int(size)

# get gene synonyms
gene_names = {}
with open('worm_genes.json') as fp:
	gene_names = json.load(fp)


# read all of the gene expression files
all_exp = {}
for file in sys.argv[1:]: all_exp[file] = read_expression(file)

# find genes that are in all files
fave_gene = set()
for gene in all_exp[list(all_exp.keys())[0]]:
	missing = False
	for file in all_exp:
		if gene not in all_exp[file]:
			missing = True
			break
	if missing: continue
	fave_gene.add(gene)

# create fav expression subset
fav_exp = {}
for file in all_exp:
	fav_exp[file] = {}
	for gene, val in all_exp[file].items():
		if gene in fave_gene:
			fav_exp[file][gene] = val

# normalize to gene length if the counts are all ints
for file in fav_exp:
	if count_type(fav_exp[file]) == int:
		for gene in fav_exp[file]:
			fav_exp[file][gene] /= gene_len[gene]
	#	print(fav_exp[file])

# turn expressions into probability distributions
for file in fav_exp:
	total = sum(fav_exp[file].values())
	for gene in fav_exp[file]:
		fav_exp[file][gene] /= total

# create pairwise distances
files = list(fav_exp.keys())
condmat = []
for i in range(len(files)):
	for j in range(i + 1, len(files)):
		condmat.append(dtc(fav_exp[files[i]], fav_exp[files[j]]))

# cluster and plot
labels = []
for file in files:
	f = file.split('/')
	labels.append(f[-1])

plt.figure(figsize=(12,8))
plt.title('RNA-seq Clustering?')
Z = linkage(condmat, method='weighted')
dn = dendrogram(Z, labels=labels, orientation='left')
plt.savefig('foo.png')
plt.close()
