
import argparse
import os
import sqlite3
import requests
import time



parser = argparse.ArgumentParser(description='metadata stuff')
parser.add_argument('database', help='database name (e.g. something.db)')
parser.add_argument('dir', help='output directory')
parser.add_argument('--delay', type=float, default=0.5,
	help='delay between requests [%(default).2f]')
parser.add_argument('--retry', type=int, default=3,
	help='number of times to retry after download failure [%(default)i]')
parser.add_argument('--verbose', action='store_true',
	help='print status messages')
arg = parser.parse_args()

os.system(f'mkdir -p {arg.dir}')

# connect to db
con = sqlite3.connect(arg.database)
cur = con.cursor()

# get a record that is not labeled and not locked
cur.execute('SELECT srx_id, gsm_id from experiment limit 25')
for srx, gsm in cur.fetchall():
	if not gsm: continue
	if os.path.exists(f'{arg.dir}/{gsm}.txt'): continue

	url = f'https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc={gsm}&targ=self&view=brief&form=text'

	success = False
	for _ in range(arg.retry + 1):
		if arg.verbose: print(f'requesting: {url}', file=sys.stderr)
		response = requests.get(url)
		if response.status_code != 200:
			print(f'ERROR {response.status_code}, will retry in {arg.delay}')
		else:
			success = True
			break
		time.sleep(arg.delay)
	if not success: sys.exit('problems with downloading GEO records')


	with open(f'{arg.dir}/{gsm}.txt', 'w') as fp:
		print(response.text, file=fp)
