import argparse
import random
import sys

from ftx import FTX
from grimoire.genome import Reader
from grimoire.toolbox import revcomp_str

def generate_reads(tx, size):
	chrom = tx.dna

	# create linear indexes
	dna = [] # dna positional index
	rna = [] # rna sequence
	for exon in tx.exons:
		for i in range(exon.length):
			coor = i + exon.beg -1 # zero-based coordinates
			dna.append(coor)
			rna.append(chrom.seq[coor])

	# read generator
	for i in range(len(rna) - size + 1):
		coor = [dna[i+j] for j in range(size)]
		exons = []
		beg = coor[0]
		seen = 0
		for j in range(size -1):
			d = coor[j+1] - coor[j]
			if d > 1:
				end = beg + j -seen
				exons.append( (beg, end) )
				seen += end - beg + 1
				beg = coor[j+1]
		exons.append( (beg, beg+j -seen +1) )
		ftx = FTX(chrom.name, tx.id, tx.strand, exons, 'r')
		read = ''.join([chrom.seq[beg:end+1] for beg, end in exons])

		yield ftx, read, len(exons)

##############################################################################

parser = argparse.ArgumentParser()
parser.add_argument('fasta', help='fasta file, compressed ok')
parser.add_argument('gff', help='gff file, compressed ok')
parser.add_argument('--readlength', type=int, default=100, metavar='<int>',
	help='[%(default)i]')
parser.add_argument('--samplegenes', type=float, default=1.0, metavar='<p>',
	help='downsample genes [%(default).3f]')
parser.add_argument('--samplereads', type=float, default=1.0, metavar='<p>',
	help='downsample reads [%(default).3f]')
parser.add_argument('--seed', type=int, default=0, metavar='<int>',
	help='set random seed')
parser.add_argument('--double', action='store_true',
	help='produce reads from both strands')
parser.add_argument('--spliced', action='store_true',
	help='only produce reads that splice')
arg = parser.parse_args()

if arg.seed != 0: random.seed(arg.seed)

genes = 0
reads = 0
bases = 0
genome = Reader(fasta=arg.fasta, gff=arg.gff)
for chrom in genome:
	for gene in chrom.ftable.build_genes():
		if not gene.is_coding(): continue
		if gene.issues: continue
		if random.random() > arg.samplegenes: continue
		genes += 1
		tx = gene.transcripts()[0] # 1 transcript per gene

		for ftx, seq, exons in generate_reads(tx, arg.readlength):
			if arg.spliced and exons == 1: continue
			if random.random() < arg.samplereads:
				print('>', ftx, '+', sep='')
				print(seq)
				reads += 1
				bases += arg.readlength
			if arg.double and random.random() < arg.samplereads:
				seq = revcomp_str(seq)
				print('>', ftx, '-', sep='')
				print(seq)
				reads += 1
				bases += arg.readlength

print(f'genes: {genes}', f'reads: {reads}', f'bases: {bases}',
	sep='\n', file=sys.stderr)

