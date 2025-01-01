import argparse
import random

COMPLEMENT = str.maketrans('ACGTRYMKWSBDHV', 'TGCAYRKMWSVHDB')

def random_exon(n, debug):
	if debug: return '[' + '#' * (n-2) + ']'
	return random_seq(n)

def random_flank(n, debug):
	if debug: return '~' * n
	return random_seq(n)

def random_intron(n, site, strand, debug):
	if strand == '-': site = site.translate(COMPLEMENT)[::-1]
	if debug: return site[0:2] + '-' * (n-4) + site[2:]
	return site[0:2] + random_seq(n-4) + site[2:]

def random_seq(n):
	return ''.join(random.choices('ACGT', k=n))

def write_gff(fp, chrom, uid, exons):
	print('\t'.join((chrom, 'bakeoff', 'gene', str(exons[0][0]+1),
		str(exons[-1][1]+1), '.', strand, '.', f'ID={uid}')), file=fp)
	print('\t'.join((chrom, 'bakeoff', 'mRNA', str(exons[0][0]+1),
		str(exons[-1][1]+1), '.', strand, '.', f'ID=tx.{uid};Parent={uid}')),
		file=fp)
	for beg, end in exons:
		print('\t'.join((chrom, 'bakeoff', 'exon', str(beg+1),
			str(end+1), '.', strand, '.', f'Parent=tx.{uid}')), file=fp)

# CLI

parser = argparse.ArgumentParser(description='Create synthetic fasta & gff',
	epilog='structure: ~~~[exon]--intron--[var.exon]--intron--[exon]~~~')
parser.add_argument('name', metavar='<name>',
	help='results in <name>.fa, <name>.gff')
parser.add_argument('--flank', type=int, default=50, metavar='<int>',
	help='flanking sequence on either side of gene [%(default)i]')
parser.add_argument('--emin', type=int, default=5, metavar='<int>',
	help='minimum middle exon length [%(default)i]')
parser.add_argument('--emax', type=int, default=50, metavar='<int>',
	help='maximum middle exon length [%(default)i]')
parser.add_argument('--estep', type=int, default=1, metavar='<int>',
	help='middle exon step size [%(default)i]')
parser.add_argument('--exon', type=int, default=100, metavar='<int>',
	help='outer exon lengths [%(default)i]')
parser.add_argument('--intron', type=int, default=50, metavar='<int>',
	help='intron lengths [%(default)i]')
parser.add_argument('--chroms', type=int, default=10, metavar='<int>',
	help='number of chromosomes in experiment [%(default)i]')
parser.add_argument('--seed', type=int, metavar='<int>',
	help='set random seed')
parser.add_argument('--double', action='store_true',
	help='create genes on both strands')
parser.add_argument('--noncanonical', action='store_true',
	help='include non-canonical splice sites')
parser.add_argument('--debug', action='store_true')
arg = parser.parse_args()

# SETUP

if arg.debug:
	arg.flank = 5
	arg.emin = 5
	arg.emax = 6
	arg.exon = 10
	arg.intron = 10
	arg.chroms = 2
	arg.noncanonical = True
	arg.double = True

if arg.seed: random.seed(arg.seed)

sites = ['GTAG']
if arg.noncanonical: sites.extend(['GCAG', 'ATAC', 'AATT'])
strands = ['+']
if arg.double: strands.append('-')

fa = open(f'{arg.name}.fa', 'w')
gff = open(f'{arg.name}.gff', 'w')

# MAIN LOOP

for c in range(arg.chroms):
	chrom = f'c{c+1}'
	print(f'>{chrom}', file=fa)
	off = 0
	gene = 0

	for i in range(arg.emin, arg.emax+1, arg.estep):
		for strand in strands:
			for site in sites:
				gene += 1
				f1 = random_flank(arg.flank, arg.debug)
				e1 = random_exon(arg.exon, arg.debug)
				i1 = random_intron(arg.intron, site, strand, arg.debug)
				ve = random_exon(i, arg.debug)
				i2 = random_intron(arg.intron, site, strand, arg.debug)
				e2 = random_exon(arg.exon, arg.debug)
				f2 = random_flank(arg.flank, arg.debug)
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
				uid = f'{chrom}.g{gene}.e{i}{strand}{site}'
				write_gff(gff, chrom, uid, exons)

fa.close()
gff.close()
