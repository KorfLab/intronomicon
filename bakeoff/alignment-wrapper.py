import argparse
import gzip
import os
import sys

import bakeoff

def run(cli, verbose):
	cli = cli.replace('\t', ' ')
	if verbose: print(cli, file=sys.stderr)
	os.system(cli)

def needfastq(arg):
	fastq = f'{arg.reads}.fq.gz'
	if not os.path.exists(fastq):
		if arg.verbose: print('creating fastq file', file=sys.stderr)
		with gzip.open(fastq, 'wt') as fp:
			for name, seq in bakeoff.readfasta(arg.reads):
				print('@', name, file=fp, sep='')
				print(seq, file=fp)
				print('+', file=fp)
				print('J' * len(seq), file=fp)
	arg.reads = fastq

def sam_to_fx(filename):
	! need to sort out sam bitwise flags
	! need to finish function
	with open(filename) as fp:
		for line in fp:
			line = fp.readline()
			if line == '': break
			if line.startswith('@'): continue
			f = line.split('\t')
			if    f[1] == '0': st = '+'
			elif f[1] == '16': st = '-'
			else:              sys.exit(f'unexpected bflag {f[1]}')
			pos   = int(f[3])
			cigar = f[5]
			print(pos, cigar, st)

def sim4_to_fx(filename):
	fp = open(filename)
	chrom = None
	strand = None
	exons = []
	n = 0

	for line in fp:
		if line.startswith('seq1 ='):
			if chrom is not None:
				print(bakeoff.fxcompose(chrom, str(n), strand, exons, info='blat'))
			chrom = None
			strand = None
			exons = []
			n += 1
			continue
		elif line.startswith('seq2 ='):
			f = line.split()
			chrom = f[2][:-1]
			continue
		f = line.split()
		if len(f) != 4: continue
		beg, end = f[1][1:-1].split('-')
		exons.append((int(beg), int(end)))
		st = '+' if f[3] == '->' else '-'
		if strand is None: strand = st
		else: assert(strand == st)
	print(bakeoff.fxcompose(chrom, str(n), strand, exons, info='blat'))
	fp.close()

parser = argparse.ArgumentParser(description='bakeoff alignment wrapper',
	epilog='programs: blat bt2 bwa hs2 mm2 star')
parser.add_argument('genome', help='genome file in FASTA format, gz ok')
parser.add_argument('reads', help='reads file in FASTA format, gz ok')
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
	# indexing: no
	# reads: fa.gz
	# output: sim4 converted to fx
	# notes:
	#   consider adding -fine -q=rna (they are used for full length mRNAs)
	#   -maxIntron=N  defaults to 750000
	#   -out=wublast blast blast8 blast9 sim4 maf axt pslx (using sim4)
	cli = f'blat {arg.genome} {arg.reads} {out} -out=sim4'
	if not arg.verbose: cli += ' > /dev/null'
	run(cli, arg.verbose)
	sim4_to_fx(out)
elif arg.prog == 'bt2':
	# indexing: yes
	# reads: fq.gz
	# output: sam
	if not os.path.exists(f'{arg.genome}.1.bt2'):
		cli = f'bowtie2-build {arg.genome} {arg.genome}'
		if not arg.verbose: cli += ' 2> /dev/null'
		run(cli, arg.verbose)
	needfastq(arg)
	cli = f'bowtie2 -x {arg.genome} -U {arg.reads} > {out}'
	if not arg.verbose: cli += ' 2> /dev/null'
	run(cli, arg.verbose)
	sam_to_fx(out)
elif arg.prog == 'bwa':
	# indexing: yes
	# reads: fa.gz or fq.gz
	# output: sam
	if not os.path.exists(f'{arg.genome}.bwt'):
		cli = f'bwa index {arg.genome}'
		if not arg.verbose: cli += ' 2> /dev/null'
		run(cli, arg.verbose)
	cli = f'bwa mem {arg.genome} {arg.reads} > {out}'
	if not arg.verbose: cli += ' 2> /dev/null'
	run(cli, arg.verbose)
	sam_to_fx(out)
elif arg.prog == 'hs2':
	# indexing: yes
	# reads: fq.gz
	# output: sam
	if not os.path.exists(f'{arg.genome}.1.ht2'):
		cli = f'hisat2-build -f {arg.genome} {arg.genome}'
		if not arg.verbose: cli += ' > /dev/null 2> /dev/null'
		run(cli, arg.verbose)
	needfastq(arg)
	cli = f'hisat2 -x {arg.genome} -U {arg.reads} > {out}'
	if not arg.verbose: cli += ' 2> /dev/null'
	run(cli, arg.verbose)
	sam_to_fx(out)
elif arg.prog == 'mm2':
	# indexing: no
	# reads: fa.gz or fq.gz
	cli = f'minimap2 -ax splice {arg.genome} {arg.reads} > {out}'
	run(cli, arg.verbose)
	sam_to_fx(out)
elif arg.prog == 'star':
	# indexing: yes
	# reads: fa.gz or fq.gz
	# output: sam
	# notes
	#   twopassMode Basic ?
	#   requires temp file cleanup
	idx = f'{arg.genome}-star'
	if not os.path.exists(idx):
		cli = f'STAR --runMode genomeGenerate --genomeDir {idx}\
			--genomeFastaFiles {arg.genome} --genomeSAindexNbases 8\
			--runThreadN {arg.threads}'
		if not arg.verbose: cli += ' > /dev/null'
		run(cli, arg.verbose)
	cli = f'STAR --genomeDir {idx}\
		--readFilesIn {arg.reads} --readFilesCommand "gunzip -c"\
		--runThreadN {arg.threads} --outFileNamePrefix {out}'
	if not arg.verbose: cli += ' > /dev/null'
	run(cli, arg.verbose)
	os.rename(f'{out}Aligned.out.sam', f'{out}')
	for x in ('Log.final.out', 'Log.out', 'Log.progress.out', 'SJ.out.tab'):
		os.unlink(f'{out}{x}')
	sam_to_fx(out)
else:
	sys.exit(f'ERROR: unknown program: {arg.prog}')

# Clean up

if not arg.debug: os.unlink(f'{out}')
