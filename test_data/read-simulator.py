import argparse
import random
import sys

from grimoire.genome import Reader
from grimoire.toolbox import revcomp_str

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
		for i in range(0, len(txs) - arg.readlength -1):
			seq = txs[i:i+arg.readlength]
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

