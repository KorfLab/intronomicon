import argparse
import io
import json
import os
import re
import requests
import sys
import time

from ncbi_reader import sraxml

def download(url, arg):
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
parser.add_argument('--debug', action='store_true',
	help='limit number of new downloads for testing purposes')
parser.add_argument('--verbose', action='store_true',
	help='print status messages')
arg = parser.parse_args()

os.system(f'mkdir -p {arg.dir}/sra {arg.dir}/gsm {arg.dir}/gse {arg.dir}/raw')

base = 'https://eutils.ncbi.nlm.nih.gov/entrez/eutils'
taxid = f'term=txid{arg.taxid}[Organism]'

# initial request to get number of records
url = f'{base}/esearch.fcgi?db=sra&{taxid}&retmax=1'
txt = download(url, arg)
n = int(re.search(r'<Count>(\d+)</Count>', txt).group(1))

# full request of all records
url = f'{base}/esearch.fcgi?db=sra&{taxid}&retmax={n}'
txt = download(url, arg)

# download each XML file and corresponding GEO files (if available)
done = 0
for m in re.finditer(r'<Id>(\d+)</Id>', txt):

	# download sra xml file
	uid = m.group(1)
	sfile = f'{arg.dir}/sra/{uid}.xml'
	sra = None # to be filled
	if have(sfile, arg):
		with open(sfile) as fp: sra = fp.read()
	else:
		url = f'{base}/efetch.fcgi?db=sra&id={uid}&rettype=xml&retmode=text'
		sra = download(url, arg)
		if sra is None: continue
		with open(sfile, 'w') as fp: fp.write(sra)

	# download gsm txt file, also get GSE id
	data, status = sraxml.read(io.StringIO(sra))
	if data is None: continue
	if not data['gsm_id']: continue
	gsm_file = f'{arg.dir}/gsm/{data["gsm_id"]}.txt'
	if have(gsm_file, arg):
		with open(gsm_file) as fp: gsm = fp.read()
	else:
		url = f'https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc={data["gsm_id"]}&targ=self&view=brief&form=text'
		gsm = download(url, arg)
		with open(gsm_file, 'w') as fp: fp.write(gsm)
	m = re.search(r'Sample_series_id = (\S+)', gsm)
	gse_id = m.group(1)

	# download gse txt file
	gse_file = f'{arg.dir}/gse/{gse_id}.txt'
	if not have(gse_file, arg):
		url = f'https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc={data["gsm_id"]}&targ=gse&view=brief&form=text'
		gse = download(url, arg)
		with open(gse_file, 'w') as fp: fp.write(gse)

	# download RAW - sometimes it doesn't exist...
	#raw_file = f'{arg.dir}/raw/{gse_id}'
	#tar_file = raw_file + '.tar'
	#if have(raw_file, arg): continue
	#url = f'https://www.ncbi.nlm.nih.gov/geo/download/?acc={gse_id}&format=file'
	#verbose = '' if arg.verbose else '--silent'
	#os.system(f'mkdir -p {raw_file}')
	#os.system(f'curl {verbose} "{url}" --output {tar_file}')
	#os.system(f'tar -xf {tar_file} -C {raw_file}')
	#os.system(f'rm -f {tar_file}')
	#time.sleep(arg.delay)

	# debug
	done += 1
	if arg.debug and done >= 10: sys.exit('debugging')
