# bakeoff.py
# various functions for performing the spliced alignment bakeoff

import gzip

def getfp(filename):
	"""gets file pointer from file, compressed file, or stdin"""
	fp = None
	if   filename.endswith('.gz'): fp = gzip.open(filename, 'rt')
	elif filename == '-':          fp = sys.stdin
	else:                          fp = open(filename)
	return fp

def fxcompose(chrom, uid, strand, exons, info=None):
	"""creates fx string from components, including exons"""
	estr = ','.join([f'{beg+1}-{end+1}' for beg, end in exons]) # 1-based
	return '|'.join((chrom, uid, strand, estr, info))

def fxparse(text):
	"""returns structred object from fx string"""
	chrom, uid, strand, estr, info = text.split('|')
	exons = []
	for s in estr.split(','):
		beg, end = s.split('-')
		exons.append((int(beg)-1, int(end)-1))
	return chrom, uid, strand, exons, info

def readfasta(filename):
	"""yields name, seq from fasta file"""
	name = None
	seqs = []
	fp = getfp(filename)
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

