import argparse
import glob
import os
import re
import sqlite3
import sys
import time

import sraxml

def soft2text(path, fields):
	with open(path) as fp: text = fp.read()
	lines = text.split('!')
	m = re.search(r'^\^\S+ = (\S+)', lines[0])
	uid = m.group(1)
	d = {}
	for line in lines[1:]:
		f = line.split()
		if len(f) == 0: continue
		f = line.split()
		k = f[0]
		v = ' '.join(f[2:])
		if k not in d: d[k] = v
		else: d[k] += f'; {v}'

	output = [f'{k}: {d[k]}' for k in fields if k in d]
	text = '. '.join(output).replace('"', '')

	return uid, text

## CLI ##

parser = argparse.ArgumentParser(description='intronomicon database builder')
parser.add_argument('name', help='database name (including .db suffix)')
parser.add_argument('dir', help='download directory (e.g. build)')
parser.add_argument('--taxid', default='6239', help='taxid [%(default)s]')
parser.add_argument('--debug', action='store_true')
arg = parser.parse_args()

if not arg.name.endswith('.db'): sys.exit('database name must end in .db')
create_tables = False if os.path.exists(arg.name) else True

## SQL ##

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
		gse_id TEXT,
		gsm_id TEXT,
		gsm_txt TEXT,
		platform TEXT,
		model TEXT,
		paired INTEGER CHECK (paired in (0, 1)),
		FOREIGN KEY (gse_id) REFERENCES series (gse_id))""",
	"""CREATE TABLE run(
		run_id TEXT PRIMARY KEY,
		nts INTEGER,
		seqs INTEGER,
		srx_id TEXT,
		aligned INTEGER CHECK (aligned in (0, 1)),
		labeled INTEGER CHECK (labeled in (0, 1)),
		locked INTEGER CHECK (locked in (0, 1)),
		FOREIGN KEY (srx_id) REFERENCES experiment (srx_id))""",
	"""CREATE TABLE series(
		gse_id TEXT PRIMARY KEY,
		gse_txt TEXT)"""
]

# create tables
if create_tables:
	for table in tables: cur.execute(table)

# gse table
n = 0
for filename in glob.glob(f'{arg.dir}/gse/*'):
	if arg.debug and n >= 10: break
	n += 1
	gse_id, gse_txt = soft2text(filename, ('Series_title', 'Series_summary',
		'Series_overall_design', 'Series_type', 'Series_sample_id'))
	rows = '(gse_id, gse_txt)'
	vals = f'("{gse_id}", "{gse_txt}")'
	try:
		cur.execute(f'INSERT OR IGNORE INTO series {rows} VALUES {vals}')
	except:
		sys.exit(f'series table write error {rows} {vals}')

# experiment and run tables
n = 0
for filename in glob.glob(f'{arg.dir}/sra/*'):
	if arg.debug and n >= 10: break

	# read sra and ensure geo
	with open(filename) as fp: data, status = sraxml.read(fp)
	if data is None: continue
	if data['taxid'] != arg.taxid:
		print(data['taxid'], arg.taxid)
		sys.exit(f'taxid mismatch ')
	if not data['gsm_id']: continue
	gsm_id, gsm_txt = soft2text(f'{arg.dir}/gsm/{data["gsm_id"]}.txt',
		('Sample_source_name_ch1', 'Sample_title' 'Sample_description',
		'Sample_molecule_ch1', 'Sample_library_selection',
		'Sample_library_strategy', 'Sample_characteristics_ch1',
		'Sample_series_id'))
	assert(gsm_id == data['gsm_id'])
	m = re.search(r'Sample_series_id: (GSE\d+)', gsm_txt) # keeping 1 GSE
	gse_id = m.group(1) # note: a sample can appear in more than one series
	n += 1

	# experiment table
	srx = data['srx_id']
	plt = data['platform']
	mod = data['model']
	par = data['paired']
	rows = '(srx_id, gse_id, gsm_id, gsm_txt, platform, model, paired)'
	vals = f'("{srx}", "{gse_id}", "{gsm_id}", "{gsm_txt}", "{plt}", "{mod}", {par})'
	try:
		cur.execute(f'INSERT OR IGNORE INTO experiment {rows} VALUES {vals}')
	except:
		sys.exit(f'experiment table write error {rows} {vals}')

	# run table
	for run in data['runs']:
		rid = run['run_id']
		nts = run['nts']
		seqs = run['seqs']
		rows = '(run_id, nts, seqs, srx_id, aligned, labeled, locked)'
		vals = f'("{rid}", {nts}, {seqs}, "{srx}", 0, 0, 0)'
		try:
			cur.execute(f'INSERT OR IGNORE INTO run {rows} VALUES {vals}')
		except:
			sys.exit(f'run table write error {rows} {vals}')

con.commit()
con.close()
