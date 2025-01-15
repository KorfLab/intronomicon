import argparse
import glob
import json
import requests
import sqlite3
import sys
import time

import sraxml

def get_gsm_text(gsm, arg):
	url = f'https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc={gsm}&targ=self&view=brief&form=text'

	success = False
	for _ in range(arg.retry + 1):
		if arg.verbose: print(f'requesting: {url}', file=sys.stderr)
		response = requests.get(url)
		if response.status_code != 200:
			print(f'ERROR {response.status_code}, will retry in {arg.delay}')
		else:
			success = True
			break
		time.sleep(arg.delay)
	if not success: sys.exit('problems with downloading GEO records')

	lines = response.text.split('!')
	for line in lines[1:]:
		f = line.split()
		k = f[0]
		v = ' '.join(f[2:])

	sys.exit('ok')

	return 'placeholder'

parser = argparse.ArgumentParser(description='create intronomicon database')
parser.add_argument('name', help='database name (including .db suffix)')
parser.add_argument('xml', help='directory of xml files, (e.g. build)')
parser.add_argument('taxid', help='taxid (e.g. 6239 is C. elegans)')
parser.add_argument('--delay', type=float, default=0.5,
	help='delay between requests [%(default).2f]')
parser.add_argument('--retry', type=int, default=3,
	help='number of times to retry after download failure [%(default)i]')
parser.add_argument('--verbose', action='store_true',
	help='print status messages')
parser.add_argument('--test', action='store_true', help='stops after 10 files')
arg = parser.parse_args()

if not arg.name.endswith('.db'): sys.exit('database name must end in .db')

con = sqlite3.connect(arg.name)
cur = con.cursor()

tables = [
	"""CREATE TABLE intron(
		intron_id TEXT PRIMARY KEY,
		chrom TEXT,
		start INTEGER,
		end INTEGER,
		count INTEGER,
		strand INTEGER CHECK (strand IN (1, -1)),
		srx_id TEXT,
		FOREIGN KEY (srx_id) REFERENCES experiment (srx_id))""",
	"""CREATE TABLE experiment(
		srx_id TEXT PRIMARY KEY,
		gsm_id TEXT,
		srx_txt TEXT,
		gsm_txt TEXT,
		platform TEXT,
		model TEXT,
		paired INTEGER CHECK (paired in (0, 1)))""",
	"""CREATE TABLE runs(
		run_id TEXT PRIMARY KEY,
		nts INTEGER,
		seqs INTEGER,
		srx_id TEXT,
		aligned INTEGER CHECK (aligned in (0, 1)),
		labeled INTEGER CHECK (labeled in (0, 1)),
		locked INTEGER CHECK (locked in (0, 1)),
		FOREIGN KEY (srx_id) REFERENCES experiment (srx_id))"""
]

# create tables
for table in tables: cur.execute(table)

# fill tables
n = 0
for filename in glob.glob(f'{arg.xml}/*'):
	if arg.test and n > 10: break

	with open(filename) as fp: data, status = sraxml.read(fp)
	if data is None: continue
	if data['taxid'] != arg.taxid:
		print(data['taxid'], arg.taxid)
		sys.exit(f'taxid mismatch ')
	n += 1

	# experiment table
	texts = [f'{k}: {v}' for k, v in data['info'].items()]
	srx = data['srx_id']
	gsm = data['geo_id'] if data['geo_id'] is not None else ''
	stx = '; '.join(texts).replace('"', '')
	gtx = '' # get_gsm_text(gsm, arg) if gsm else ''
	plt = data['platform']
	mod = data['model']
	par = data['paired']
	rows = '(srx_id, gsm_id, srx_txt, gsm_txt, platform, model, paired)'
	vals = f'("{srx}", "{gsm}", "{stx}", "{gtx}", "{plt}", "{mod}", {par})'
	try:
		cur.execute(f'INSERT OR IGNORE INTO experiment {rows} VALUES {vals}')
	except:
		sys.exit(f'experiment table write error {rows} {vals}')

	# runs table
	for run in data['runs']:
		rid = run['run_id']
		nts = run['nts']
		seqs = run['seqs']
		rows = '(run_id, nts, seqs, srx_id, aligned, labeled, locked)'
		vals = f'("{rid}", {nts}, {seqs}, "{srx}", 0, 0, 0)'
		try:
			cur.execute(f'INSERT OR IGNORE INTO runs {rows} VALUES {vals}')
		except:
			sys.exit(f'runs table write error {rows} {vals}')

con.commit()
con.close()
