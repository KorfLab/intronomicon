
import argparse
import json
import os
import sqlite3

parser = argparse.ArgumentParser(description='metadata training set maker')
parser.add_argument('database', help='database name (e.g. something.db)')
parser.add_argument('dir', help='output directory')
arg = parser.parse_args()

# connect to db
con = sqlite3.connect(arg.database)
cur = con.cursor()

# get all gse
gse = {}
cur.execute('SELECT * from series')
for gse_id, gse_txt in cur.fetchall():
	gse[gse_id] = {'gse': gse_id, 'text': gse_txt, 'exps': []}

# add experiment
cur.execute('SELECT * from experiment')
for srx_id, gse_id, gsm_id, gsm_txt, plt, mod, par in cur.fetchall():
	gse[gse_id]['exps'].append({'srx': srx_id, 'gsm': gsm_id, 'text': gsm_txt})

os.system(f'mkdir -p {arg.dir}')
for gse_id in gse:
	with open(f'{arg.dir}/{gse_id}.json', 'w') as fp:
		print(json.dumps(gse[gse_id], indent=4), file=fp)
