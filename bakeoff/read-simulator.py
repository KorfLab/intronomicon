import argparse
import random
import sys

from grimoire.genome import Reader
from grimoire.toolbox import revcomp_str

def show_introns(genome, exons):
	for i in range(len(exons) -1):
		ib = exons[i][1]
		ie = exons[i+1][0]
		print('>', genome[ib:ie], '<')

def read_seq(genome, exons):
	seqs = []
	for beg, end in exons: seqs.append(genome[beg-1:end])
	return ''.join(seqs)

def read_name(exons):
	name = []
	for beg, end in exons: name.append(f'{beg}-{end}')
	return ','.join(name)


def read_coords(tx, pos, rlen):
	# collect exons in rna and dna coordinates
	rna_exons = []
	dna_exons = []
	used = 0
	sizes = []
	for exon in tx.exons:
		size = exon.end - exon.beg + 1
		sizes.append(len(exon.seq_str()))
		rna_exons.append((used, used+size))
		dna_exons.append((exon.beg-1, exon.end))
		used += size

	# find the exon where the read starts and its offset within the exon
	idx = None
	tot = 0
	off = None
	for i, (rbeg, rend) in enumerate(rna_exons):
		if pos >= rbeg and pos < rend:
			idx = i
			off = pos - tot
			break
		tot += rend - rbeg
	assert(idx is not None)

	# create human-readable (1-offset) alignment coordinates
	align = []
	avail = rlen
	for (beg, end) in dna_exons[idx:]:
		elen = end - beg
		if len(align) == 0:
			a = beg + off + 1
			b = min(end, beg + off + avail)
		else:
			if avail > elen:  a, b = beg+1, end
			else:             a, b = beg+1, beg + avail
		align.append((a, b))
		avail -= b - a + 1
		if avail < 1: break

	return align

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
parser.add_argument('--paired', action='store_true',
	help='produce paired reads (not yet supported)')
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
		txs = tx.tx_str()
		for i in range(0, len(txs) - arg.readlength -1):
			exons = read_coords(tx, i, arg.readlength)
			if len(exons) == 1 and arg.spliced: continue
			seq = read_seq(tx.dna.seq, exons)
			name = '|'.join((str(i+1), tx.id, tx.strand, read_name(exons)))

			if random.random() < arg.samplereads:
				print(f'>{name}|+', seq, sep='\n')
				reads += 1
				bases += arg.readlength
			if arg.double and random.random() < arg.samplereads:
				seq = revcomp_str(seq)
				print(f'>{name}|-', seq, sep='\n')
				reads += 1
				bases += arg.readlength

print(f'genes: {genes}', f'reads: {reads}', f'bases: {bases}',
	sep='\n', file=sys.stderr)

