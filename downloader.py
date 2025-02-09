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


###########################################################################

parser = argparse.ArgumentParser(description='Get GEO & SRA files from NCBI')
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

subs = ('gds_txt', 'gse_soft', 'gsm_soft', 'sra_xml')
for sub in subs: os.system(f'mkdir -p {arg.dir}/{sub}')

##############################################################################
# PART 1: Retrieve GEO data for taxid
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

# download each GSE family soft file
done = 0
failures = []
for m in re.finditer(r'<Id>(\d+)</Id>', txt):

	# get the uid for each gds record that match (2nnnnnnnn)
	uid = m.group(1)
	if len(uid) != 9: continue
	if not uid.startswith('2'): continue
	gds_file = f'{arg.dir}/gds_txt/{uid}.txt'
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
	
	# fetch the gse soft file
	gse_id = m.group(1)
	soft_url = f'{ftp_url}soft/{gse_id}_family.soft.gz'
	gse_file = f'{arg.dir}/gse_soft/{gse_id}_family.soft.gz'
	try:
		run(f'wget -P {arg.dir}/gse_soft "{soft_url}"', check=True, shell=True)
		if os.path.getsize(gse_file) > 100000: # has embedded microarray shit
			os.unlink(gse_file)
			os.system(f'touch {gse_file}')
	except: # wget fails (URL has no soft file)
		failures.append(('soft-fail', soft_url))
		continue
	
	# save the gds
	with open(gds_file, 'w') as fp: print(gds_text, file=fp)	
	if arg.debug: break

print('gds failures', failures)

#sys.exit('end of part 1')


##############################################################################
# PART 2: Select the GSEs worth exploring
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
	
	# retrieve GSMs (just the good ones, is this even necessary?)
	for k, sdata in d.items():
		if not k.startswith('SAMPLE'): continue
		gsm_id = should_be_one(sdata['Sample_geo_accession'])
		gsm_file = f'{arg.dir}/gsm_soft/{gsm_id}.txt'
		if os.path.exists(gsm_file):
			with open(gsm_file) as fp: gsm_txt = fp.read()
		else:
			base = f'https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc='
			query = f'{gsm_id}&targ=self&view=brief&form=text'
			gsm_txt = retrieve(f'{base}{query}', arg)
			with open(gsm_file, 'w') as fp: fp.write(gsm_txt)
		

sys.exit('end of part 2')



##############################################################################
# PART 3: Download XML from SRA
##############################################################################


base = 'https://eutils.ncbi.nlm.nih.gov/entrez/eutils'
taxid = f'term=txid{arg.taxid}[Organism]'

# initial request to get number of records
url = f'{base}/esearch.fcgi?db=sra&{taxid}&retmax=1'
txt = download(url, arg)
n = int(re.search(r'<Count>(\d+)</Count>', txt).group(1))

# full request of all records
url = f'{base}/esearch.fcgi?db=sra&{taxid}&retmax={n}'
txt = download(url, arg)

# download each XML file
done = 0
for m in re.finditer(r'<Id>(\d+)</Id>', txt):

	# download sra xml file
	uid = m.group(1)
	sfile = f'{arg.dir}/sra_xml/{uid}.xml'
	sra = None # to be filled
	if os.path.exists(sfile):
		with open(sfile) as fp: sra = fp.read()
	else:
		url = f'{base}/efetch.fcgi?db=sra&id={uid}&rettype=xml&retmode=text'
		sra = download(url, arg)
		if sra is None: continue
		with open(sfile, 'w') as fp: fp.write(sra)



##############################################################################
# PART 4: Connect best GSEs to SRAs
##############################################################################