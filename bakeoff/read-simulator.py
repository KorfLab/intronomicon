import argparse
import random
import sys

from grimoire.genome import Reader
from grimoire.toolbox import revcomp_str

def genomic_coords(tx, seq, pos):
	rna_off = []
	dna_off = []
	used = 0
	for exon in tx.exons:
		size = (exon.end - exon.beg) + 1
		rna_off.append( (used, used+size))
		dna_off.append( (exon.beg, exon.beg + size))
		used += size
	print(rna_off)
	print(dna_off)

	remain = len(seq)
	rexons = []
	gexons = []
	beg = pos
	end = pos + remain
	for (rbeg, rend), (dbeg, dend) in zip(rna_off, dna_off):
		#print('considering', beg, rbeg, rend)
		if beg >= rbeg and beg <= rend:
			#print('found overlapping exon starting at', rbeg, rend)
			avail = rend - beg
			if remain > avail:
				end = beg + avail
				rexons.append( (beg, end) )
				gexons.append( (beg+dbeg, end+dbeg) )
				beg = end + 1
				remain -= avail
				#print('created exon, sequence remains')
			else:
				end = beg + remain
				rexons.append( (beg, end) )
				gexons.append( (beg+dbeg, end+dbeg) )
				#print('created exon, done')
				break

	rna = tx.seq_str()
	dna = tx.dna.seq.lower()

	print(len(rna), len(dna))
	print(rexons)
	print(gexons)
	for (rbeg, rend), (gbeg, gend) in zip(rexons, gexons):
		print('lengths', rend-rbeg, gend-gbeg)
		rseq = rna[rbeg:rend]
		gseq = dna[gbeg:gend]
		print(rbeg, rend, gbeg, gend)
		print(rseq)
		print(gseq)

	sys.exit('here')
	return exons


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
		tx = gene.transcripts()[0]
		txs = tx.seq_str()
		start = 0
		size = 250
		seq = txs[start:start+size]
		exons = genomic_coords(tx, seq, start)
		print(exons)
		sys.exit()





		txs = tx.seq_str()
		for i in range(0, len(txs) - arg.readlength -1):
			seq = txs[i:i+arg.readlength]
			exons = genomic_coords(tx, seq, i)
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

