import argparse
import glob
import gzip
import os
import sys
import time
import requests
import json

from ncbi_reader import soft_read

def findkeys(d, prefix):
	tags = []
	for key in d:
		if key.startswith(prefix): tags.append(key)
	return tags

def should_be_one(things):
	if len(things) == 1: return things[0]
	if len(things) == 0: sys.exit('should have existed')
	sys.exit('more than one thing found')

def retrieve(url, arg):
	for _ in range(arg.retry + 1):
		if arg.verbose: print(f'requesting: {url}', file=sys.stderr)
		response = requests.get(url)
		response.encoding = 'utf-8'
		time.sleep(arg.delay)
		if response.status_code != 200:
			if arg.verbose: print(f'ERROR {response.status_code}, will retry in {arg.delay}')
			continue
		return response.text
	sys.exit(f'download ({url}) failed...')


parser = argparse.ArgumentParser(description='Get SRA files from NCBI')
parser.add_argument('dir', help='download directory (e.g. build)')
parser.add_argument('--taxid', default='6239',
	help='NCBI taxonomy identifier [e.g. C. elegans is %(default)s]')
parser.add_argument('--delay', type=float, default=0.5,
	help='delay between download requests [%(default).2f]')
parser.add_argument('--retry', type=int, default=3,
	help='number of times to retry after download failure [%(default)i]')
parser.add_argument('--verbose', action='store_true',
	help='print status messages')
parser.add_argument('--debug', action='store_true')
arg = parser.parse_args()

subs = ('gsm', 'srx')
for sub in subs: os.system(f'mkdir -p {arg.dir}/{sub}')

for path in glob.glob(f'{arg.dir}/soft/*'):
	if os.path.getsize(path) == 0: continue # damn micorarry blobs
	d = soft_read(path)
		
	# remove anything that isn't sequencing
	plats = findkeys(d, 'PLATFORM')
	techs = set()
	taxids = set()
	for plat in plats:
		techs.add(d[plat]['Platform_technology'][0])
		for taxid in d[plat]['Platform_taxid']: taxids.add(taxid)
	if len(techs) > 1: continue
	if 'high-throughput sequencing' not in techs: continue
	if len(taxids) > 1: continue
	if arg.taxid not in taxids: sys.exit(f'taxid mismatch')

	# skip anything that isn't normal RNA-Seq
	srcs = set()
	strats = set()
	bios = []
	sras = []
	weird = False
	for k, sdata in d.items():
		if not k.startswith('SAMPLE'): continue
		gsm_id = should_be_one(sdata['Sample_geo_accession'])
		src = should_be_one(sdata['Sample_library_source'])
		strat = should_be_one(sdata['Sample_library_strategy'])
		srcs.add(src)
		strats.add(strat)

		b = []
		s = []
		for text in sdata['Sample_relation']:
			if text.startswith('BioSample'): b.append(text[10:-1])
			if text.startswith('SRA'): s.append(text[4:-1])
		
		if len(b) == 1: bios.append(b[0])
		else: weird = True
		if len(s) == 1: sras.append(s[0])
		else: weird = True
	
	if weird: continue
	if len(srcs) != 1 or len(strats) != 1: continue # multi-omic something
	if 'transcriptomic' not in srcs: continue
	if 'RNA-Seq' not in strats: continue
	if len(bios) != len(sras): sys.exit('wtf')
	
	# find files that may be gene expression counts
	not_wanted = ('RAW.tar', '.bw', '.bigwig', '.pdf', '.txt')
	keep = []
	skey = findkeys(d, 'SERIES')[0]
	if 'Series_supplementary_file' not in d[skey]: continue
	for filename in d[skey]['Series_supplementary_file']:
		skip = False
		for ext in not_wanted:
			if filename.endswith(ext):
				skip = True
				break
		if skip: continue
		keep.append(filename)
	if len(keep) == 0: continue
	
	# retrieve GSMs
	for k, sdata in d.items():
		if not k.startswith('SAMPLE'): continue
		gsm_id = should_be_one(sdata['Sample_geo_accession'])
		gsm_file = f'{arg.dir}/gsm/{gsm_id}.txt'
		if os.path.exists(gsm_file):
			with open(gsm_file) as fp: gsm_txt = fp.read()
		else:
			base = f'https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc='
			query = f'{gsm_id}&targ=self&view=brief&form=text'
			gsm_txt = retrieve(f'{base}{query}', arg)
			with open(gsm_file, 'w') as fp: fp.write(gsm_txt)
		
