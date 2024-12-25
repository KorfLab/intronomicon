import argparse
import gzip
import sys

def readfasta(filename):
	name = None
	seqs = []

	fp = None
	if   filename.endswith('.gz'): fp = gzip.open(filename, 'rt')
	elif filename == '-':          fp = sys.stdin
	else:                          fp = open(filename)

	while True:
		line = fp.readline()
		if line == '': break
		line = line.rstrip()
		if line.startswith('>'):
			if len(seqs) > 0:
				seq = ''.join(seqs)
				yield(name, seq)
				name = line[1:]
				seqs = []
			else:
				name = line[1:]
		else:
			seqs.append(line)
	yield(name, ''.join(seqs))
	fp.close()

parser = argparse.ArgumentParser('convert sam to pairwise alignment')
parser.add_argument('ref', help='reference genome')
parser.add_argument('sam', help='sam file')
arg = parser.parse_args()

def readsam(filename):
	fp = None
	if   filename.endswith('.gz'): fp = gzip.open(filename, 'rt')
	elif filename == '-':          fp = sys.stdin
	else:                          fp = open(filename)

	while True:
		line = fp.readline()
		if line == '': break
		if line.startswith('@'): continue
		f = line.split('\t')
		qname = f[0]
		bflag = int(f[1])
		pos   = int(f[3])
		cigar = f[5]
		seq   = f[9]
		yield qname, seq, pos, cigar, bflag


chrom = {}
for name, seq in readfasta(arg.ref):
	chrom[name] = seq

for query, seq, pos, cigar, bflag in readsam(arg.sam):
	print(cigar)
