import argparse
import glob
import sqlite3
import sys

import sraxml

parser = argparse.ArgumentParser(description='create intronomicon database')
parser.add_argument('name', help='database name (including .db suffix)')
parser.add_argument('xml', help='directory of xml files, (e.g. build)')
parser.add_argument('taxid', help='taxid (e.g. 6239 is C. elegans)')
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
		exp_id TEXT,
		FOREIGN KEY (exp_id) REFERENCES experiment (exp_id))""",
	"""CREATE TABLE experiment(
		exp_id TEXT PRIMARY KEY,
		geo_id TEXT,
		platform TEXT,
		model TEXT,
		paired INTEGER CHECK (paired in (0, 1)),
		info TEXT)""",
	"""CREATE TABLE runs(
		run_id TEXT PRIMARY KEY,
		nts INTEGER,
		seqs INTEGER,
		exp_id TEXT,
		FOREIGN KEY (exp_id) REFERENCES experiment (exp_id))"""
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
	xid = data['srx_id']
	gid = data['geo_id']
	plt = data['platform']
	mod = data['model']
	par = data['paired']
	freetext = []
	for tag, txt in data['info'].items(): freetext.append(f'{tag}: {txt}')
	info = '\n'.join(freetext)
	info = info.replace('"', "'")
	rows = '(exp_id, geo_id, platform, model, paired, info)'
	vals = f'("{xid}", "{gid}", "{plt}", "{mod}", {par}, "{info}")'
	try:
		cur.execute(f'INSERT OR IGNORE INTO experiment {rows} VALUES {vals}')
	except:
		sys.exit(f'experiment table write error {rows} {vals}')

	# runs table
	for run in data['runs']:
		rid = run['run_id']
		nts = run['nts']
		seqs = run['seqs']
		rows = '(run_id, nts, seqs, exp_id)'
		vals = f'("{rid}", {nts}, {seqs}, "{xid}")'
		try:
			cur.execute(f'INSERT OR IGNORE INTO runs {rows} VALUES {vals}')
		except:
			sys.exit(f'runs table write error {rows} {vals}')

con.commit()
con.close()
