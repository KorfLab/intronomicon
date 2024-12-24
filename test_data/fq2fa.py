import gzip
import sys

if len(sys.argv) == 1: sys.exit('error: provide fq or fq.gz files')

for file in sys.argv[1:]:
	if   file.endswith('fq'): fp = open(file)
	elif file.endswith('fq.gz'): fp = gzip.open(file, 'rt')
	while True:
		h = fp.readline()
		if h == '': break
		s = fp.readline()
		p = fp.readline()
		q = fp.readline()
		print('>', h[1:], sep='', end='')
		print(s, end='')
	fp.close()
