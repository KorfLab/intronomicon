import argparse
import glob
import json
import os
import re
import sqlite3
import sys
import time

from ncbi_reader import sraxml_read, soft_read

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
		gsm_id TEXT,
		gse_id TEXT,
		platform TEXT,
		model TEXT,
		paired INTEGER CHECK (paired in (0, 1)),
		FOREIGN KEY (gse_id) REFERENCES series (gse_id))""",
	"""CREATE TABLE run(
		srr_id TEXT PRIMARY KEY,
		nts INTEGER,
		rlen INTEGER,
		srx_id TEXT,
		aligned INTEGER CHECK (aligned in (0, 1)),
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
for filename in glob.glob(f'{arg.dir}/gse_soft/*'):
	if arg.debug and n >= 10: break
	n += 1
	thing = soft_read(filename)
	if not thing: continue # e.g. microarray
	skey = [x for x in thing.keys() if x.startswith('SERIES')][0]
	series = thing[skey]
	
	# want 1 of platform taxid and sample taxid
	taxids = set()
	check = ('Series_platform_taxid', 'Series_sample_taxid')
	for tag in check:
		if tag not in series: continue
		for taxid in series[tag]: taxids.add(taxid)
	if len(taxids) > 1: continue # for now
	if arg.taxid not in taxid: sys.exit('wtf')
	gse_id = series['Series_geo_accession'][0]
	gse_txt = series['Series_title'][0].replace('""', "'")
	rows = '(gse_id, gse_txt)'
	vals = f'("{gse_id}", "{gse_txt}")'
	try:
		cur.execute(f'INSERT OR IGNORE INTO series {rows} VALUES {vals}')
		con.commit()
	except:
		sys.exit(f'series table write error {rows} {vals}')

# experiment and run tables
n = 0
for filename in glob.glob(f'{arg.dir}/sra_xml/*'):
	n += 1
	if arg.debug and n >= 10: break

	# read sra and ensure geo
	with open(filename) as fp: data, status = sraxml_read(fp)
	if data is None: continue
	if data['taxid'] != arg.taxid:
		print(data['taxid'], arg.taxid)
		sys.exit(f'taxid mismatch ')
	if not data['gsm_id']: continue
	gsm_file = f'{arg.dir}/gsm_soft/{data["gsm_id"]}.txt'
	if not os.path.exists(gsm_file): continue
	
	thing = soft_read(gsm_file)
	skey = [x for x in thing.keys() if x.startswith('SAMPLE')][0]
	sample = thing[skey]
	gsm_id = sample['Sample_geo_accession'][0]
	gsm_txt = sample['Sample_characteristics_ch1'][0].replace('""', "'")
	gse_id = sample['Sample_series_id'][0] # could be more than one?
	#print(sample['Sample_series_id']) yes, but dead before here
	
	# experiment table
	srx = data['srx_id']
	plt = data['platform']
	mod = data['model']
	par = data['paired']
	rows = '(srx_id, gsm_id, gse_id, platform, model, paired)'
	vals = f'("{srx}", "{gsm_id}", "{gse_id}", "{plt}", "{mod}", {par})'
	try:
		cur.execute(f'INSERT OR IGNORE INTO experiment {rows} VALUES {vals}')
		con.commit()
	except:
		sys.exit(f'experiment table write error {rows} {vals}')
	
	# run table
	for run in data['runs']:
		rid = run['run_id']
		nts = run['nts']
		rlen = int(run['nts'] / run['seqs'])
		rows = '(srr_id, nts, rlen, srx_id, aligned, locked)'
		vals = f'("{rid}", {nts}, {rlen}, "{srx}", 0, 0)'
		try:
			cur.execute(f'INSERT OR IGNORE INTO run {rows} VALUES {vals}')
			con.commit()
		except:
			sys.exit(f'run table write error {rows} {vals}')

con.close()
