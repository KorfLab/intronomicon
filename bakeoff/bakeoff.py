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
	stuff = [chrom, uid, strand, estr]
	if info: stuff.append(info)
	return '|'.join(stuff)

def fxparse(text):
	"""returns structred object from fx string"""
	f = text.split('|')
	chrom, uid, strand, estr = f[:4]
	if len(f) > 4: xtra = '|'.join(f[4:])
	else:          xtra = None
	exons = []
	for s in estr.split(','):
		beg, end = s.split('-')
		exons.append((int(beg)-1, int(end)-1))
	return chrom, uid, strand, exons, xtra

def fxread(filename):
	"""returns dictionary of chromosomes with lists of transcripts"""
	fp = getfp(filename)
	data = {}
	for line in fp:
		chrom, uid, st, exons, info = fxparse(line)
		if chrom not in data: data[chrom] = []
		data[chrom].append({'id':uid, 'st':st, 'info':info, 'exons':exons})
	return data

def fxoverlap(fx1, fx2):
	"""compares 2 fx objects"""

	# compare chromosomes - should be unneccesary
	if fx1.chrom != fx2.chrom: return False

	# compare transcript extent
	b1 = fx1.exons[0][0]
	e1 = fx1.exons[-1][1]
	b2 = fx2.exons[0][0]
	e2 = fx2.exons[-1][1]

	# compare interiors

	#if b2 <= e1 and e2 >= b1: return True

	# is this just a true/false?
	# is this even necessary?

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

