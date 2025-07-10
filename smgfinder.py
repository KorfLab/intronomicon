import glob
import json
from ncbi_reader import soft_read, sraxml_read
import sys


for filename in glob.glob('build/gse_soft/*'):
	soft = soft_read(filename)
	if not soft: continue
	smg = False
	for k1 in soft:
		for k2 in soft[k1]:
			for text in soft[k1][k2]:
				if 'smg' in text:
					smg = True
					break
			if smg: break
		if smg: break
	
	if smg:
		skey = [x for x in soft.keys() if x.startswith('SERIES')][0]
		print(skey)
