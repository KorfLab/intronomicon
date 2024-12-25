import argparse
import gzip
import os
import sys

def fqgz2fa(file, out):
	fp = gzip.open(file, 'rt')
	f2 = open(out, 'w')
	while True:
		h = fp.readline()
		if h == '': break
		s = fp.readline()
		p = fp.readline()
		q = fp.readline()
		print('>', h[1:], sep='', end='', file=f2)
		print(s, end='', file=f2)
	fp.close()
	f2.close()

def run(cli, verbose):
	cli = cli.replace('\t', ' ')
	if verbose: print(cli, file=sys.stderr)
	os.system(cli)

def readfasta(filename):
	name = None
	seqs = []

	fp = None
	if   filename.endswith('.gz'): fp = gzip.open(filename, 'rt')
	elif filename == '-':          fp = sys.stdin
	else:                          fp = open(filename)

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

def readsam(filename):
	fp = None
	if   filename.endswith('.gz'): fp = gzip.open(filename, 'rt')
	elif filename == '-':          fp = sys.stdin
	else:                          fp = open(filename)

	while True:
		line = fp.readline()
		if line == '': break
		if line.startswith('@'): continue
		f = line.split('\t')
		qname = f[0]
		bflag = int(f[1])
		strand = '+' if f'{bflag:012b}'[7] == '0' else '-'
		pos   = int(f[3])
		cigar = f[5]
		seq   = f[9]
		yield qname, seq, pos, strand, cigar

##############################################################################

parser = argparse.ArgumentParser(description='read alignment wrapper',
	epilog='programs: blat bt2 bwa hs2 mm2 star')
parser.add_argument('genome', help='genome file in FASTA format')
parser.add_argument('reads', help='reads file in fq.gz format')
parser.add_argument('prog', help='program name')
parser.add_argument('--threads', type=int, default=1,
	help='number of cores/threads [%(default)i]')
parser.add_argument('--debug', action='store_true', help='keep SAM file')
parser.add_argument('--verbose', action='store_true')
arg = parser.parse_args()

# Run the aligner

out = f'temp-{os.getpid()}'
sam = f'{out}.sam'

if arg.prog == 'blat':
	sys.exit('blat is acutally not supported')
	# create fa, no index
	fa = arg.reads[:-5] + 'fa'
	if not os.path.exists(fa):
		fqgz2fa(arg.reads, fa)
	# align
	cli = f'blat {arg.genome} {fa} {out} -out -q=rna'
	if not arg.verbose: cli += ' > /dev/null'
	run(cli, arg.verbose)
	# convert
	cli = f'psl2sam.pl...'
	# consider adding -fine -q=rna (they are used for full length mRNAs)
	# -maxIntron=N  defaults to 750000
	# -out=wublast blast blast8 blast9 sim4 maf axt pslx
elif arg.prog == 'bt2':
	# index
	if not os.path.exists(f'{arg.genome}.1.bt2'):
		cli = f'bowtie2-build {arg.genome} {arg.genome}'
		if not arg.verbose: cli += ' 2> /dev/null'
		run(cli, arg.verbose)
	# align, no clean
	cli = f'bowtie2 -x {arg.genome} -U {arg.reads} > {sam}'
	if not arg.verbose: cli += ' 2> /dev/null'
	run(cli, arg.verbose)
elif arg.prog == 'bwa':
	# index
	if not os.path.exists(f'{arg.genome}.bwt'):
		cli = f'bwa index {arg.genome}'
		if not arg.verbose: cli += ' 2> /dev/null'
		run(cli, arg.verbose)
	# align, no clean
	cli = f'bwa mem {arg.genome} {arg.reads} > {sam}'
	if not arg.verbose: cli += ' 2> /dev/null'
	run(cli, arg.verbose)
elif arg.prog == 'hs2':
	# index
	if not os.path.exists(f'{arg.genome}.1.ht2'):
		cli = f'hisat2-build -f {arg.genome} {arg.genome}'
		if not arg.verbose: cli += ' > /dev/null 2> /dev/null'
		run(cli, arg.verbose)
	# align
	cli = f'hisat2 -x {arg.genome} -U {arg.reads} > {sam}'
	if not arg.verbose: cli += ' 2> /dev/null'
	run(cli, arg.verbose)
elif arg.prog == 'mm2':
	# no index, no clean
	cli = f'minimap2 -ax splice {arg.genome} {arg.reads} > {sam}'
	run(cli, arg.verbose)
elif arg.prog == 'star':
	# index
	idx = f'{arg.genome}-star'
	if not os.path.exists(idx):
		cli = f'STAR --runMode genomeGenerate --genomeDir {idx}\
			--genomeFastaFiles {arg.genome} --genomeSAindexNbases 8\
			--runThreadN {arg.threads}'
		if not arg.verbose: cli += ' > /dev/null'
		run(cli, arg.verbose)
	# align
	cli = f'STAR --genomeDir {idx}\
		--readFilesIn {arg.reads} --readFilesCommand "gunzip -c"\
		--runThreadN {arg.threads} --outFileNamePrefix {out}'
	if not arg.verbose: cli += ' > /dev/null'
	run(cli, arg.verbose)
	# 2nd pass...
	# clean
	os.rename(f'{out}Aligned.out.sam', f'{out}.sam')
	for x in ('Log.final.out', 'Log.out', 'Log.progress.out', 'SJ.out.tab'):
		os.unlink(f'{out}{x}')
else:
	sys.exit(f'ERROR: unknown program: {arg.prog}')

# Explore alignments...
chrom = {}
for name, seq in readfasta(arg.genome):
	chrom[name] = seq

for query, seq, pos, strand, cigar in readsam(f'{sam}'):
	print(cigar, strand)


# Clean up

if not arg.debug: os.unlink(f'{sam}')
