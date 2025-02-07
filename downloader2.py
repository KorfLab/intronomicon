import argparse
import io
import json
import os
import re
import requests
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

def have(file, arg):
	if os.path.exists(file) and os.path.getsize(file) > 0:
		if arg.verbose: print(f'skipping previous download {file}')
		return True
	return False

parser = argparse.ArgumentParser(description='Get SRA and GEO files from NCBI')
parser.add_argument('dir', help='output directory (e.g. build)')
parser.add_argument('--taxid', default='6239',
	help='NCBI taxonomy identifier [e.g. C. elegans is %(default)s]')
parser.add_argument('--delay', type=float, default=0.5,
	help='delay between download requests [%(default).2f]')
parser.add_argument('--retry', type=int, default=3,
	help='number of times to retry after download failure [%(default)i]')
parser.add_argument('--verbose', action='store_true',
	help='print status messages')
arg = parser.parse_args()

subs = ('soft', 'srx')
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
for m in re.finditer(r'<Id>(\d+)</Id>', txt):
	uid = m.group(1)
	url = f'{base}/efetch.fcgi?db=gds&id={uid}'
	text = retrieve(url, arg)
	m = re.search(r'(ftp:\S+)', text)
	ftp_url = m.group(1)
	m = re.search(r'n\/(GSE\d+)', text)
	gse_id = m.group(1)
	soft_url = f'{ftp_url}soft/{gse_id}_family.soft.gz'
	file = f'{arg.dir}/soft/{gse_id}_family.soft.gz'
	if not have(file, arg): os.system(f'curl "{url}" --output {file}')

	sys.exit()
# next up: sample and run files
# can I get these from an alternate URL or do I need to use efetch?