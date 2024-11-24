import argparse
import glob
import sqlite3
import sys

import intronomicon

parser = argparse.ArgumentParser(description='create intronomicon database')
parser.add_argument('--name', default='intronomicon.db',
	help='database name [%(default)s]')
parser.add_argument('--xml', default='build',
	help='directory of xml files [%(default)s]')
parser.add_argument('--taxid', default='6239',
	help='verify the correct taxon [%(default)s]')
parser.add_argument('--test', action='store_true')
arg = parser.parse_args()

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

	with open(filename) as fp: data, status = intronomicon.read_sra_xml(fp)
	if data is None: continue
	if data['taxid'] != arg.taxid:
		print(data['taxid'], arg.taxid)
		sys.exit(f'taxid mismatch ')
	n += 1

	# experiment table
	xid = data['srx_id']
	plt = data['platform']
	mod = data['model']
	par = data['paired']
	freetext = []
	for tag, txt in data['info'].items(): freetext.append(f'{tag}: {txt}')
	info = '\n'.join(freetext)
	info = info.replace('"', "'")
	rows = '(exp_id, platform, model, paired, info)'
	vals = f'("{xid}", "{plt}", "{mod}", {par}, "{info}")'
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
