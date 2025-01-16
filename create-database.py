import argparse
import glob
import json
import requests
import sqlite3
import sys
import time

import sraxml

def get_gsm_text(gfile, arg):
	with open(gfile) as fp: text = fp.read()
	lines = text.split('!')

	d = {}
	for line in lines[1:]:
		f = line.split()
		if len(f) == 0: continue
		f = line.split()
		k = f[0]
		v = ' '.join(f[2:])
		if k not in d: d[k] = v
		else: d[k] += f'; {v}'

	# difficulties
	if d['Sample_channel_count'] != '1': return None
	if ',' in d['Sample_taxid_ch1']: return None

	# sometimes present
	some = ('Sample_description', 'Sample_growth_protocol_ch1',
		'Sample_treatment_protocol_ch1')
	for s in some:
		if s not in d: d[s] = ''

	# construct blurb
	output = []
	output.append(f"Source: {d['Sample_source_name_ch1']}")
	output.append(f"Title: {d['Sample_title']}")
	output.append(f"Description: {d['Sample_description']}")
	output.append(f"Molecule: {d['Sample_molecule_ch1']}")
	output.append(f"Selection: {d['Sample_library_selection']}")
	output.append(f"Strategy: {d['Sample_library_strategy']}")
	output.append(f"Characteristics: {d['Sample_characteristics_ch1']}")
	output.append(f"Growth: {d['Sample_growth_protocol_ch1']}")
	output.append(f"Treatment: {d['Sample_treatment_protocol_ch1']}")
	return '. '.join(output).replace('"', '')


parser = argparse.ArgumentParser(description='create intronomicon database')
parser.add_argument('name', help='database name (including .db suffix)')
parser.add_argument('dir', help='directory of sra/geo files, (e.g. build)')
parser.add_argument('taxid', help='taxid (e.g. 6239 is C. elegans)')
parser.add_argument('--debug', action='store_true')
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
for filename in glob.glob(f'{arg.dir}/sra/*'):
	if arg.debug and n >= 10: break

	# read sra and geo
	with open(filename) as fp: data, status = sraxml.read(fp)
	if data is None: continue
	if data['taxid'] != arg.taxid:
		print(data['taxid'], arg.taxid)
		sys.exit(f'taxid mismatch ')
	if not data['geo_id']: continue
	gtx = get_gsm_text(f'{arg.dir}/geo/{data["geo_id"]}.txt', arg)
	if gtx is None: continue
	n += 1

	# experiment table
	texts = [f'{k}: {v}' for k, v in data['info'].items()]
	srx = data['srx_id']
	gsm = data['geo_id']
	stx = '; '.join(texts).replace('"', '')
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
