import argparse
import glob
import io
import json
import os
import re
import requests
from subprocess import run
import sys
import time

from ncbi_reader import sraxml_read, soft_read

def remote_filesize(url, arg):
	try:
		if arg.verbose: print(f'checking filesize at {url}')
		text = run(f'wget --spider "{url}"', check=True, shell=True,
			capture_output=True).stderr.decode()
	except:
		print(f'remote_filesize wget failed on {url}', file=sys.stderr)
		return None
		
	m = re.search(r'Length: (\d+) ', text)
	if not m:
		print(f'remote_filesize grep failed on {text}', file=sys.stderr)
		return None
		
	size = int(m.group(1))
	if arg.verbose: print(f'size: {size} bytes')	
	
	time.sleep(arg.delay)
	return size

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

def findkeys(d, prefix):
	tags = []
	for key in d:
		if key.startswith(prefix): tags.append(key)
	return tags

def should_be_one(things):
	if len(things) == 1: return things[0]
	if len(things) == 0: sys.exit('should have existed')
	sys.exit('more than one thing found')

##############################################################################
# CLI
##############################################################################

parser = argparse.ArgumentParser(description='Get GEO & SRA files from NCBI')
parser.add_argument('dir', help='output directory (e.g. build)')
parser.add_argument('--taxid', default='6239',
	help='NCBI taxonomy identifier [e.g. C. elegans is %(default)s]')
parser.add_argument('--delay', type=float, default=0.5,
	help='delay between download requests [%(default).2f]')
parser.add_argument('--retry', type=int, default=5,
	help='number of times to retry after download failure [%(default)i]')
parser.add_argument('--max-gse-size', type=float, default=1e7,
	help='maximum size of gse.soft file [%(default).0f]')
parser.add_argument('--max-supp-size', type=float, default=1e8,
	help='maximum size of supplementary file [%(default).0f]')
parser.add_argument('--verbose', action='store_true',
	help='print status messages')
parser.add_argument('--debug', action='store_true')
arg = parser.parse_args()

subs = ('gds_txt', 'gse_soft', 'gse_supp', 'gsm_soft', 'sra_xml')
for sub in subs: os.system(f'mkdir -p {arg.dir}/{sub}')

##############################################################################
# PART 0: Retrieve all GEO data series txt files for specific taxid
##############################################################################

base = 'https://eutils.ncbi.nlm.nih.gov/entrez/eutils'
taxid = f'term=txid{arg.taxid}[Organism]'

# initial request to get number of geo dataset records
url = f'{base}/esearch.fcgi?db=gds&{taxid}&retmax=1'
txt = retrieve(url, arg)
n = int(re.search(r'<Count>(\d+)</Count>', txt).group(1))

# full request of all records
url = f'{base}/esearch.fcgi?db=gds&{taxid}&retmax={n}'
txt = retrieve(url, arg)

# download all gds text files
for m in re.finditer(r'<Id>(\d+)</Id>', txt):

	# get the uid for each gds record that match (2nnnnnnnn)
	uid = m.group(1)
	if len(uid) != 9: continue
	if not uid.startswith('2'): continue
	gds_file = f'{arg.dir}/gds_txt/{uid}.txt'
	if os.path.exists(gds_file):
		if arg.verbose: print('skipping', gds_file)
		continue
	gds_text = retrieve(url, arg)
	with open(gds_file, 'w') as fp: print(gds_text, file=fp)
	
	if arg.debug:
		print('stopping to debug')
		break
	
##############################################################################
# PART 1: Retrieve all GEO series soft files
##############################################################################

for filename in glob.glob(f'{arg.dir}/gds_txt/*'):
	with open(filename) as fp: gds_text = fp.read()
	m = re.search(r'ftp:..(\S+)/', gds_text)
	if not m:
		if arg.verbose: print(f'no ftp in {filename}')
		continue
	
	ftp_url = m.group(1)
	m = re.search(r'n\/(GSE\d+)', gds_text)
	gse_id = m.group(1)
	soft_url = f'{ftp_url}/soft/{gse_id}_family.soft.gz'
	gse_file = f'{arg.dir}/gse_soft/{gse_id}_family.soft.gz'
	
	if os.path.exists(gse_file):
		if arg.verbose: print('skipping', gse_file)
		continue
	
	size = remote_filesize(soft_url, arg)
	if size > arg.max_gse_size:
		if arg.verbose: print(f'gse too large {size}')
		os.system(f'touch {gse_file}')
		continue
	
	try:
		cli = f'wget -P {arg.dir}/gse_soft "{soft_url}"'
		if arg.verbose: print(cli, flush=True)
		run(cli, check=True, shell=True)
		time.sleep(arg.delay)
	except:
		if arg.verbose: print(f'wget failed on {soft_url}')
		continue
		
	if arg.debug:
		print('stopping to debug')
		break

##############################################################################
# PART 2: Download GSE supps and GSMs but filter
#  - RNA-Seq
#  - Not mixtures of RNA-Seq with other stuff
#  - Have supplementary gene expression files
##############################################################################

for path in glob.glob(f'{arg.dir}/gse_soft/*'):
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
	weird = False
	for k, sdata in d.items():
		if not k.startswith('SAMPLE'): continue
		gsm_id = should_be_one(sdata['Sample_geo_accession'])
		src = should_be_one(sdata['Sample_library_source'])
		strat = should_be_one(sdata['Sample_library_strategy'])
		srcs.add(src)
		strats.add(strat)
	if 'transcriptomic' not in srcs: continue
	if 'RNA-Seq' not in strats: continue

	# find files that may be gene expression counts
	not_wanted = ('RAW.tar', '.bw', '.bigwig', '.pdf', 'pdf.gz', '.txt',
		'.bedgraph.gz', '.fasta.gz', '.wig.gz')
	keep = []
	skey = findkeys(d, 'SERIES')[0]
	if 'Series_supplementary_file' not in d[skey]:
		if arg.verbose: print(f'no supps in {path}')
		continue
		
	for filename in d[skey]['Series_supplementary_file']:
		skip = False
		for ext in not_wanted:
			if filename.endswith(ext):
				skip = True
				break
		if skip: continue		
		keep.append(filename)
	if len(keep) == 0: continue
	
	# download supps
	gse_id = skey.split()[2]
	gse_dir = f'{arg.dir}/gse_supp/{gse_id}'
	os.system(f'mkdir -p {gse_dir}')
	for url in keep:
		filename = url.split('/')[-1]
		filepath = f'{gse_dir}/{filename}'
		
		if os.path.exists(filepath):
			if arg.verbose: print('skipping', filepath)
			continue

		size = remote_filesize(url, arg)
		if size is None or size > arg.max_supp_size:
			os.system(f'touch {filepath}') # don't attempt download again
			continue
			
		run(f'wget -P {gse_dir} "{url}"', check=True, shell=True)
		time.sleep(arg.delay)
		
	# retrieve GSMs
	for k, sdata in d.items():
		if not k.startswith('SAMPLE'): continue
		gsm_id = should_be_one(sdata['Sample_geo_accession'])
		gsm_file = f'{arg.dir}/gsm_soft/{gsm_id}.txt'
		if os.path.exists(gsm_file):
			if arg.verbose: print('skipping', gsm_file)
			continue
		else:
			base = f'https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc='
			query = f'{gsm_id}&targ=self&view=brief&form=text'
			gsm_txt = retrieve(f'{base}{query}', arg)
			with open(gsm_file, 'w') as fp: fp.write(gsm_txt)

##############################################################################
# PART 3: Download XML from SRA
##############################################################################

base = 'https://eutils.ncbi.nlm.nih.gov/entrez/eutils'
taxid = f'term=txid{arg.taxid}[Organism]'

# initial request to get number of records
url = f'{base}/esearch.fcgi?db=sra&{taxid}&retmax=1'
txt = retrieve(url, arg)
n = int(re.search(r'<Count>(\d+)</Count>', txt).group(1))

# full request of all records
url = f'{base}/esearch.fcgi?db=sra&{taxid}&retmax={n}'
txt = retrieve(url, arg)

# download each XML file
done = 0
for m in re.finditer(r'<Id>(\d+)</Id>', txt):

	# download sra xml file
	uid = m.group(1)
	sfile = f'{arg.dir}/sra_xml/{uid}.xml'
	sra = None # to be filled
	if os.path.exists(sfile):
		if arg.verbose: print('skipping', sfile)
		continue
	else:
		url = f'{base}/efetch.fcgi?db=sra&id={uid}&rettype=xml&retmode=text'
		sra = retrieve(url, arg)
		if sra is None: continue
		with open(sfile, 'w') as fp: fp.write(sra)

