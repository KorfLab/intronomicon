import gzip
import re
import sys
from ftx import FTX

class BitFlag:
	def __init__(self, val):
		i = int(val)
		b = f'{i:012b}'
		self.read_unmapped = True if b[-3] == '1' else False
		self.read_reverse_strand = True if b[-5] == '1' else False
		self.not_primary_alignment = True if b[-9] == '1' else False
		self.supplementary_alignment = True if b[-12] == '1' else False
		self.otherflags = []
		for i in (1, 2, 4, 6, 7, 8, 10, 11):
			if b[-i] == '1': self.otherflags.append(i)

def cigar_to_exons(cigar, pos):
	exons = []
	beg = 0
	end = 0
	for match in re.finditer(r'(\d+)([A-Z])', cigar):
		n = int(match.group(1))
		op = match.group(2)
		if   op == 'M': end += n
		elif op == 'D': pass
		elif op == 'I': end += n
		#elif op == 'S': beg += n # is this right?
		#elif op == 'H': pass
		elif op == 'N':
			exons.append((pos+beg-1, pos+end-2))
			beg = end + n
			end = beg
	exons.append((pos+beg-1, pos+end-2))
	return exons

def sam_to_ftx(filename, out, name):
	n = 0
	with open(filename) as fp:
		for line in fp:
			if line == '': break
			if line.startswith('@'): continue
			f = line.split('\t')
			qname = f[0]
			bf = BitFlag(f[1])
			chrom = f[2]
			pos   = int(f[3])
			cigar = f[5]

			st = '-' if bf.read_reverse_strand else '+'
			if bf.read_unmapped: continue
			if bf.not_primary_alignment: continue
			if bf.supplementary_alignment: continue
			if bf.otherflags:
				print(bf.otherflags)
				sys.exit('unexpected flags found, debug time')
			n += 1
			exons = cigar_to_exons(cigar, pos)
			ftx = FTX(chrom, str(n), st, exons, f'{name} ref:{qname}')
			print(ftx, file=out)

if __name__ == '__main__':

	if len(sys.argv) != 3: sys.exit(f'usage: {sys.argv[0]} <sam file> <name>')
	sam_to_ftx(sys.argv[1], sys.stdout, sys.argv[2])
