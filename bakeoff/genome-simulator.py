import argparse
import random
import sys

COMPLEMENT = str.maketrans('ACGTRYMKWSBDHV', 'TGCAYRKMWSVHDB')

def test_flank(n): return 'f' * n
def test_exon(n): return 'e' * n

def random_seq(n):
	return ''.join(random.choices('ACGT', k=n))

def random_intron(n, site, strand):
	if strand == '-': site = site.translate(COMPLEMENT)[::-1]
	return site[0:2] + random_seq(n-4) + site[2:]

def write_gff(fp, chrom, uid, exons):
	print('\t'.join((chrom, 'bakeoff', 'gene', str(exons[0][0]+1),
		str(exons[-1][1]+1), strand, '.', f'ID={uid}')), file=fp)
	print('\t'.join((chrom, 'bakeoff', 'mRNA', str(exons[0][0]+1),
		str(exons[-1][1]+1), strand, '.', f'ID=tx:{uid};Parent={uid}')),
		file=fp)
	for beg, end in exons:
		print('\t'.join((chrom, 'bakeoff', 'exon', str(beg+1),
			str(end+1), strand, '.', f'Parent=tx:{uid}')), file=fp)


parser = argparse.ArgumentParser()
parser.add_argument('basename', help='directory for files')
parser.add_argument('--flank', type=int, default=100, metavar='<int>',
	help='flanking sequence on either side of gene [%(default)i]')
parser.add_argument('--emin', type=int, default=5, metavar='<int>',
	help='minimum variable exon length [%(default)i]')
parser.add_argument('--emax', type=int, default=25, metavar='<int>',
	help='maximum variable exon length [%(default)i]')
parser.add_argument('--estep', type=int, default=1, metavar='<int>',
	help='variable exon step size [%(default)i]')
parser.add_argument('--imin', type=int, default=5, metavar='<int>',
	help='minimum variable intron length [%(default)i]')
parser.add_argument('--imax', type=int, default=50, metavar='<int>',
	help='maximum variable intron length [%(default)i]')
parser.add_argument('--istep', type=int, default=5, metavar='<int>',
	help='variable intron step size [%(default)i]')
parser.add_argument('--fixed', type=int, default=100, metavar='<int>',
	help='exon or intron length when fixed length [%(default)i]')
parser.add_argument('--chroms', type=int, default=10, metavar='<int>',
	help='number of chromosomes in experiment [%(default)i]')
parser.add_argument('--seed', type=int, metavar='<int>',
	help='set random seed')
parser.add_argument('--double', action='store_true',
	help='create genes on both strands')
parser.add_argument('--noncanonical', action='store_true',
	help='include non-canonical splice sites')
arg = parser.parse_args()
if arg.seed: random.seed(arg.seed)

sites = ['GTAG']
if arg.noncanonical: sites.extend(['GCAG', 'ATAC', 'AATT'])
strands = ['+']
if arg.double: strands.append('-')


fa = open(f'{arg.basename}.fa', 'w')
gff = open(f'{arg.basename}.gff', 'w')

for c in range(arg.chroms):
	chrom = f'chr{c+1}'
	print(f'>{chrom}', file=fa)
	off = 0

	# variable intron: ~~~[exon]--v.intron--[exon]~~~
	# intron may have variable splice sites
	for i in range(arg.imin, arg.imax+1, arg.istep):
		for strand in strands:
			for site in sites:
				f1 = random_seq(arg.flank)
				e1 = random_seq(arg.fixed)
				vi = random_intron(i, site, strand)
				e2 = random_seq(arg.fixed)
				f2 = random_seq(arg.flank)
				e1b = len(f1) + off
				e1e = e1b + len(e1)
				e2b = e1e + len(vi)
				e2e = e2b + len(e2)
				exons = ( (e1b,e1e-1), (e2b,e2e-1) )
				seq = f1 + e1 + vi + e2 + f2
				off += len(seq)
				print(seq, file=fa)
				uid = f'{chrom}:vi:{i}:{site}{strand}'
				write_gff(gff, chrom, uid, exons)

	# variable exon: ~~~[exon]--intron--[v.exon]--intron--[exon]~~~
	for i in range(arg.emin, arg.emax+1, arg.estep):
		for strand in strands:
			f1 = random_seq(arg.flank)
			e1 = random_seq(arg.fixed)
			i1 = random_intron(arg.fixed, 'GTAG', strand)
			ve = random_seq(i)
			i2 = random_intron(arg.fixed, 'GTAG', strand)
			e2 = random_seq(arg.fixed)
			f2 = random_seq(arg.flank)
			e1b = len(f1) + off
			e1e = e1b + len(e1)
			e2b = e1e + len(i1)
			e2e = e2b + len(e2)
			e3b = e2e + len(i2)
			e3e = e3b + len(e2)
			exons = ( (e1b,e1e-1), (e2b,e2e-1), (e3b,e3e-1) )
			seq = f1 + e1 + i1 + ve + i2 + e2 + f2
			off += len(seq)
			print(seq, file=fa)
			uid = f'{chrom}:ve:{i}:{strand}'
			write_gff(gff, chrom, uid, exons)

fa.close()
gff.close()
