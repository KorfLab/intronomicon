import glob
import gzip
import os
import sys

import json

from ncbi_reader import soft_read


for filename in glob.glob(f'build/gse/*'):
	d = soft_read(filename)
	print(d['Series_type'])
		
	