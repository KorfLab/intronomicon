import glob
import gzip
import os
import sys

import h5py

if len(sys.argv) == 1: sys.exit(f'usage: {sys.argv[0]} <dir>')

kill = (
	'.bw',      # bigwig is not gene-based
	'.bed.gz',  # bed is not gene-based
	'.hic',     # Hi-C isn't RNA-Seq
)

"""
Various file formats

- h5 GSE150135
- jam.txt.gz sam with gene names, would need to do the counting GSE187348
- mtx.gz GSE208154
- xlsx GSE184415
- txt.gz GSE151828
- sf.gz GSE199717
- gtf.trk.txt.gz GSE178335
- tsv.gz GSE130730
- _RSEM.txt.gz GSE200912
- bam.tdf GSE137267
- gtf.gz GSE143599
- out.tab.gz GSE252593
- bigwig GSE147401
- gtf.track.gz
- narrowPeak.gz
- bedgraph.gz GSE98919
"""

for gsedir in glob.glob(f'{sys.argv[1]}/*'):
	for filename in glob.glob(f'{gsedir}/*'):
		for target in kill:
			if filename.endswith(target):
				os.unlink(filename)
				break
		# so many names for legal files
		# probably have to look inside to examine gene names
