import argparse
import json
import re
import sqlite3
import sys

def commonsize(sizes):
	n = {}
	for size in sizes:
		if size not in n: n[size] = 0
		n[size] += 1
	common = []
	for size, c in n.items():
		if c > 0.05 * len(sizes) or c > 100: common.append(f'{size}:{c}')
	return common

def set_tagval(d, s):
	tag, val = s.split(':', maxsplit=1)
	tag = tag.rstrip().lstrip()
	val = val.rstrip().lstrip()
	if tag in d: d[tag] += '; {val}'
	d[tag] = val

actions = ('stats', 'training', 'metamorph')

parser = argparse.ArgumentParser(description='reads intronomica and stuff')
parser.add_argument('database', help='database name (e.g. something.db)')
parser.add_argument('action', help='the thing you want to do')
arg = parser.parse_args()

if arg.action not in actions: sys.exit(f'error: unknown action, try {actions}')

# connect to db
con = sqlite3.connect(arg.database)
cur = con.cursor()

if arg.action == 'stats':
	tables = ('series', 'experiment', 'run', 'intron')
	print("Counts")
	for table in tables:
		cur.execute(f'SELECT COUNT(*) FROM {table}')
		print(table, cur.fetchall()[0][0])

	print("Platforms")
	cur.execute(f'SELECT rlen, platform FROM run INNER JOIN experiment ON experiment.srx_id = run.srx_id')
	pcount = {}
	plen = {}
	for rlen, plat in cur.fetchall():
		if plat not in pcount: pcount[plat] = 0
		pcount[plat] += 1
		if plat not in plen: plen[plat] = []
		plen[plat].append(rlen)
	for p, n in pcount.items(): print(p, n)

	print("Read lengths")
	for p, a in plen.items(): print(p, min(a), max(a), commonsize(a))

if arg.action == 'training':
	d = {}
	cur.execute(f'SELECT gse_id, rlen, platform FROM run INNER JOIN experiment ON experiment.srx_id = run.srx_id where platform = "ILLUMINA"')
	for gse_id, rlen, plat in cur.fetchall():
		if gse_id not in d: d[gse_id] = []
		d[gse_id].append(rlen)
	for gse_id, sizes in d.items():
		if min(sizes) < 100: continue
		print(gse_id, len(sizes), min(sizes), max(sizes), sep='\t')

if arg.action == 'metamorph':
	cur.execute(f'SELECT gse_id, gsm_id, gsm_txt FROM run INNER JOIN experiment ON experiment.srx_id = run.srx_id where platform = "ILLUMINA" and rlen >= 100')
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

	print(json.dumps(d, indent=2))
