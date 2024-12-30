import argparse
import os
import re
import requests
import sys
import time

parser = argparse.ArgumentParser(description='get SRA XML from NCBI')
parser.add_argument('dir', help='output directory')
parser.add_argument('--taxid', default='6239',
	help='NCBI taxonomy identifier [%(default)s]')
parser.add_argument('--delay', type=float, default=0.5,
	help='delay between requests [%(default).2f]')
parser.add_argument('--retry', type=int, default=3,
	help='number of times to retry after download failure [%(default)i]')
parser.add_argument('--limit', type=int,
	help='limit number of new downloads')
parser.add_argument('--verbose', action='store_true',
	help='print status messages')
arg = parser.parse_args()

base = 'https://eutils.ncbi.nlm.nih.gov/entrez/eutils'
taxid = f'term=txid{arg.taxid}[Organism]'

# initial request to get number of records
url = f'{base}/esearch.fcgi?db=sra&{taxid}&retmax=1'
if arg.verbose: print(f'requesting: {url}', file=sys.stderr)
response = requests.get(url)
response.encoding = 'utf-8'
n = int(re.search(r'<Count>(\d+)</Count>', response.text).group(1))
time.sleep(arg.delay)

# full request of all records
url = f'{base}/esearch.fcgi?db=sra&{taxid}&retmax={n}'
if arg.verbose: print(f'requesting: {url}', file=sys.stderr)
response = requests.get(url)
response.encoding = 'utf-8'
time.sleep(arg.delay)

# download each XML file
done = 0
fail = []
for m in re.finditer(r'<Id>(\d+)</Id>', response.text):
	uid = m.group(1)
	path = f'{arg.dir}/{uid}.xml'
	if os.path.exists(path) and os.path.getsize(path) > 0:
		if arg.verbose: print(f'skipping previous download {path}')
		continue

	# download file
	url = f'{base}/efetch.fcgi?db=sra&id={uid}&rettype=xml&retmode=text'
	success = False
	for _ in range(arg.retry + 1):
		if arg.verbose: print(f'requesting: {url}', file=sys.stderr)
		response = requests.get(url)
		if response.status_code != 200:
			print(f'ERROR {response.status_code}, will retry in {arg.delay}')
		else:
			with open(path, 'w') as fp: fp.write(response.text)
			success = True
			break
		time.sleep(arg.delay)
	if success: done += 1
	else:       fail.append(uid)
	if arg.limit and done >= arg.limit: break

if arg.verbose:
	print(f'downloaded {done}, failed {len(fail)}', file=sys.stderr)
	if fail: print(fail, file=sys.stderr)
