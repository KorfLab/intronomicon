import argparse
import glob
import json
import re
import sys

import sraxml
import korflab

parser = argparse.ArgumentParser(description='xml reader test')
parser.add_argument('dir', help='directory of xml files')
parser.add_argument('--test', action='store_true', help='limit to 999 files')
arg = parser.parse_args()

log = {'count': 0, 'RNA-Seq': 0, 'GEO': 0}
printed = set()
for filename in glob.glob(f'{arg.dir}/*'):
	if arg.test and log['count'] >= 999: break
	log['count'] += 1

	with open(filename) as fp: data, status = sraxml.read(fp)
	with open(filename) as fp: raw = korflab.read_xml(fp)
	if data is None:
		if status not in log: log[status] = 0
		log[status] += 1
	else:
		log['RNA-Seq'] += 1
		stuff = []
		for tag, val in data['info'].items():
			stuff.append(f'{tag}: {val}')
		if data['srx_id'] in printed: continue
		if data['geo_id'] is None: continue
		log['GEO'] += 1
		printed.add(data['srx_id'])

		year = int(data['runs'][0]['date'][:4])
		txt = ('; '.join(stuff))
		rsize = int(data['runs'][0]['nts'] / data['runs'][0]['seqs'])

		if len(txt) < 100: continue
		if year < 2014: continue
		if rsize < 100: continue

		print('\t'.join((data['srx_id'], data['geo_id'], txt)))

print(log)
