import argparse
import io
import json
import os
import re
import requests
from subprocess import run
import sys
import time

from ncbi_reader import sraxml_read

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


parser = argparse.ArgumentParser(description='Get GEO files from NCBI')
parser.add_argument('dir', help='output directory (e.g. build)')
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

subs = ('gds', 'soft')
for sub in subs: os.system(f'mkdir -p {arg.dir}/{sub}')

base = 'https://eutils.ncbi.nlm.nih.gov/entrez/eutils'
taxid = f'term=txid{arg.taxid}[Organism]'

# initial request to get number of geo dataset records
url = f'{base}/esearch.fcgi?db=gds&{taxid}&retmax=1'
txt = retrieve(url, arg)
n = int(re.search(r'<Count>(\d+)</Count>', txt).group(1))

# full request of all records
url = f'{base}/esearch.fcgi?db=gds&{taxid}&retmax={n}'
txt = retrieve(url, arg)

# download each GSE family soft file
done = 0
failures = []
for m in re.finditer(r'<Id>(\d+)</Id>', txt):

	# get the uid for each gds record that match (2nnnnnnnn)
	uid = m.group(1)
	if len(uid) != 9: continue
	if not uid.startswith('2'): continue
	gds_file = f'{arg.dir}/gds/{uid}.txt'
	if os.path.exists(gds_file):
		if arg.verbose: print('skipping', gds_file)
		continue
	
	# get the ftp location for the GSE from the gds record
	url = f'{base}/efetch.fcgi?db=gds&id={uid}'
	gds_text = retrieve(url, arg)
	m = re.search(r'(ftp:\S+)', gds_text)
	if not m: # no ftp file listed
		failures.append(('no-ftp-file', url))
		continue
	ftp_url = m.group(1)
	m = re.search(r'n\/(GSE\d+)', gds_text)
	if not m: # no GSE file listed
		failures.append(('ftp-fail', ftp_url))
		continue
	
	gse_id = m.group(1)
	soft_url = f'{ftp_url}soft/{gse_id}_family.soft.gz'
	file = f'{arg.dir}/soft/{gse_id}_family.soft.gz'
	try:
		run(f'wget -P {arg.dir}/soft "{soft_url}"', check=True, shell=True)
		if os.path.getsize(file) > 100000: # has embedded microarray shit
			os.unlink(file)
			os.system(f'touch {file}')
	except: # wget fails (URL has no soft file)
		failures.append(('soft-fail', soft_url))
		continue
	
	# save the gds
	with open(gds_file, 'w') as fp: print(gds_text, file=fp)	
	if arg.debug: break

print('gds failures', failures)
