import argparse
import random
import sys

from grimoire.genome import Reader
from grimoire.toolbox import revcomp_str

def read_coords(tx, pos, rlen):
	rna_exons = []
	dna_exons = []
	#rseq = [] # not needed eventually
	used = 0
	for exon in tx.exons:
		size = (exon.end - exon.beg) + 1
		rna_exons.append((used, used+size))
		dna_exons.append((exon.beg-1, exon.end-1))
		#rseq.append(tx.dna.seq[exon.beg-1:exon.end]) # remove later
		used += size
	#print(rna_exons)
	#print(dna_exons)
	#print(rseq)

	# find the exon where the read starts and its offset within the exon
	idx = None
	tot = 0
	off = None
	for i, (rbeg, rend) in enumerate(rna_exons):
		if pos >= rbeg and pos <= rend:
			idx = i
			off = pos - tot
			tot += rbeg
			break
	assert(idx is not None)
	#print(f'exon-index:{idx} offset:{off}')

	# create alignment coordinates
	align = []
	avail = rlen
	for (beg, end) in dna_exons[idx:]:
		elen = end - beg + 1
		if len(align) == 0:
			a = beg + off
			b = min(end, beg + off + avail)
		else:
			if avail > elen:  a, b = beg, end
			else:             a, b = beg, beg + avail
		align.append((a, b))
		avail -= b - a
		if avail < 1: break

	return align


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
		if gene.strand == '-': continue # DEBUG
		if not gene.is_coding(): continue
		if gene.issues: continue
		if random.random() > arg.samplegenes: continue
		genes += 1
		tx = gene.transcripts()[0]
		txs = tx.seq_str()
		for i in range(0, len(txs) - arg.readlength -1):

			exons = read_coords(tx, i, arg.readlength)
			if len(exons) == 1 and arg.spliced:
			#	print(f'skipping {i}')
				continue
			print('using coordinate', i)
			for f in tx.exons: print('exon', f.beg-1, f.end, f. strand, f.end-f.beg)
			for f in tx.introns: print('intron', f.beg-1, f.end, f.strand, f.end-f.beg)

			print(exons)
			sys.exit()

			if random.random() < arg.samplereads:
				print(f'@{tx.id}|{i}+', seq, '+', 'J'*arg.readlength, sep='\n')
				reads += 1
				bases += arg.readlength
			if arg.double and random.random() < arg.samplereads:
				seq = revcomp_str(seq)
				print(f'@{tx.id}|{i}-', seq, '+', 'J'*arg.readlength, sep='\n')
				reads += 1
				bases += arg.readlength

print(f'genes: {genes}', f'reads: {reads}', f'bases: {bases}',
	sep='\n', file=sys.stderr)

