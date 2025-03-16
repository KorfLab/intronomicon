import argparse
import glob
import json
import re
import sys

def tx2gene(uid):
	f = uid.split('.')
	if   len(f) < 2: return ''
	elif len(f) == 2: clone, gene = f
	elif len(f) == 3: clone, gene, iso = f
	else: sys.exit(f'wtf {f}')

	if   gene[-1].isdigit(): return f'{clone}.{gene}'
	elif gene[-1].isalpha(): return f'{clone}.{gene[:-1]}'

parser = argparse.ArgumentParser()
parser.add_argument('file', help='file of gene names, use - for stdin')
parser.add_argument('index', help='index file (json)')
parser.add_argument('--build', help='build index file')
arg = parser.parse_args()

names = {}
if arg.build:
	for filename in glob.glob(f'{arg.build}/*'):
		with open(filename) as fp: gene = json.load(fp)
		ids = set()
		for db in gene:
			if 'display_id' in db: ids.add(db['display_id'])
			if 'primary_id' in db: ids.add(db['primary_id'])
		wbid = None
		for uid in ids:
			if uid.startswith('WBGene'):
				wbid = uid
				break
		if wbid not in names: names[wbid] = []
		for uid in ids:
			if uid != wbid: names[wbid].append(uid)
	with open(arg.index, 'w') as fp:
		print(json.dumps(names, indent=2), file=fp)
else:
	with open(arg.index) as fp: names = json.load(fp)

lookup = {}
for wbid in names:
	for xid in names[wbid]:
		if xid not in lookup: lookup[xid] = []
		lookup[xid].append(wbid)
	lookup[wbid] = [wbid] # just in case

if arg.file == '-': fp = sys.stdin
else: fp = open(arg.file)

missing = set()
multiple = set()
for line in fp:
	line = line.rstrip()
	f = line.split(maxsplit=1)
	if len(f) == 1:
		uid = line
		stuff = None
	else: uid, stuff = f

	if uid not in lookup:
		uid = tx2gene(uid) # try transcript ID instead
		if uid not in lookup:
			missing.add(uid)
			continue

	if len(lookup[uid]) > 1:
		multiple.add(uid)
		continue
	if stuff is None:
		print(lookup[uid])
	else:
		print(lookup[uid][0], stuff, sep='\t')


print('missing:', missing, file=sys.stderr)
print('multiple:', multiple, file=sys.stderr)
