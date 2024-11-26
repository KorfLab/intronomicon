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

log = {'count': 0, 'RNA-Seq': 0}
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
		printed.add(data['srx_id'])
		print(data['srx_id'], end='\t')
		print('; '.join(stuff))
		#if 'strain' in data['info']:
		#	print(data['info']['strain'])
		#else:
		#	pass
		#	#for tag in data['info']:
			#	if re.search('tissue', tag, re.IGNORECASE):
			#		print('found related', tag)
			#print(data['info'].keys())

		#print(data['runs'][0]['date'])
		#print(json.dumps(data, indent=4))
		#print(json.dumps(raw, indent=2))
#print(json.dumps(log, indent=4), file=sys.stderr)
