import argparse
import random
import sys

def random_seq(n):
	return ''.join(random.choices('ACGT', k=n))

def random_intron(n, canonical=True):
	if canonical: return 'GT' + random_seq(n-4) + 'AG'
	else:         return 'AA' + random_seq(n-4) + 'TT'


parser = argparse.ArgumentParser()
parser.add_argument('basename', help='directory for files')
parser.add_argument('--flank', type=int, default=100, metavar='<int>',
	help='flanking sequence on either side of gene [%(default)i]')
parser.add_argument('--min', type=int, default=5, metavar='<int>',
	help='minimum exon or intron length when variable [%(default)i]')
parser.add_argument('--max', type=int, default=50, metavar='<int>',
	help='maximum exon or intron length when variable [%(default)i]')
parser.add_argument('--fixed', type=int, default=100, metavar='<int>',
	help='exon or intron length when fixed [%(default)i]')
parser.add_argument('--chroms', type=int, default=10, metavar='<int>',
	help='number of chromosomes in experiment [%(default)i]')
parser.add_argument('--seed', type=int, metavar='<int>', help='random seed')
arg = parser.parse_args()
if arg.seed: random.seed(arg.seed)

fa = open(f'{arg.basename}.fa', 'w')
fx = open(f'{arg.basename}.fx', 'w')

for chrom in range(arg.chroms):
	print(f'>{chrom}', file=fa)
	off = 0
	# f1-[e1]-i1-[e2]-v.i2-[e3]-i3-[v.e4]-i4-[e5]-f2
	for vi2 in range(arg.min, arg.max+1):
		for ve4 in range(arg.min, arg.max+1):
			f1 = random_seq(arg.flank)
			e1 = random_seq(arg.fixed)
			i1 = random_intron(arg.fixed)
			e2 = random_seq(arg.fixed)
			i2 = random_intron(vi2)
			e3 = random_seq(arg.fixed)
			i3 = random_intron(arg.fixed)
			e4 = random_seq(ve4)
			i4 = random_intron(arg.fixed)
			e5 = random_seq(arg.fixed)
			f2 = random_seq(arg.flank)
			seq = ''.join((f1, e1, i1, e2, i2, e3, i3, e4, i4, e5, f2))
			e1a = off + len(f1)
			e1b = e1a + len(e1)
			e2a = e1b + len(i1)
			e2b = e2a + len(e2)
			e3a = e2b + len(i2)
			e3b = e3a + len(e3)
			e4a = e3b + len(i3)
			e4b = e4a + len(e4)
			e5a = e4b + len(i4)
			e5b = e5a + len(e5)
			exons = [f'{a}-{b}' for a, b in ((e1a,e1b), (e2a,e2b), (e3a,e3b), (e4a,e4b), (e5a,e5b))]
			print(f'{chrom}', ','.join(exons), sep='|', file=fx)
			print(seq, file=fa)
			off += len(seq)
	
fa.close()
fx.close()

