import glob
import gzip
import sys

if len(sys.argv) == 1: sys.exit(f'usage: {sys.argv[0]} <dir>')

kill = ('.beg.gz', '.h5', '.bw')

for gsedir in glob.glob(f'{sys.argv[1]}/*'):
	for filename in glob.glob(f'{gsedir}/*'):
		for target in kill:
			if filename.endswith(target):
				print('removing', filename)
				break
		# so many names for legal files
		# probably have to look inside to examine gene names
