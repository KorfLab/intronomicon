import glob
import json
import sys

import intronomicon

log = {'count': 0}
for filename in glob.glob('build/*'):
	#if log['count'] >= 999: break # testing
	log['count'] += 1

	with open(filename) as fp: data, status = intronomicon.read_sra_xml(fp)
	if data is None:
		if status not in log: log[status] = 0
		log[status] += 1
	else:
		#print(json.dumps(data, indent=4))
		for run in data['runs']: print(run['srr_id'])
#print(json.dumps(log, indent=4))
