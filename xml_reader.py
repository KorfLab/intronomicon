import argparse
import glob
import json
import sys

import intronomicon

parser = argparse.ArgumentParser(description='xml reader test')
parser.add_argument('dir', help='directory of xml files')
parser.add_argument('--test', action='store_true', help='limit to 999 files')
arg = parser.parse_args()

log = {'count': 0, 'RNA-Seq': 0}
for filename in glob.glob(f'{arg.dir}/*'):
	if arg.test and log['count'] >= 999: break
	log['count'] += 1

	with open(filename) as fp: data, status = intronomicon.read_sra_xml(fp)
	if data is None:
		if status not in log: log[status] = 0
		log[status] += 1
	else:
		log['RNA-Seq'] += 1
		print(json.dumps(data, indent=4))
print(json.dumps(log, indent=4), file=sys.stderr)
