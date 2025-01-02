import argparse
import gzip
import os
import re
import sys

from ftx import FTX
from sam2ftx import sam_to_ftx

def run(cli, arg):
	cli = cli.replace('\t', ' ')
	if arg.verbose: print(cli, file=sys.stderr)
	os.system(cli)

def getfp(filename):
	fp = None
	if   filename.endswith('.gz'): fp = gzip.open(filename, 'rt')
	elif filename == '-':          fp = sys.stdin
	else:                          fp = open(filename)
	return fp

def readfasta(filename):
	name = None
	seqs = []
	fp = getfp(filename)
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

def needfastq(arg):
	fastq = f'{arg.reads[0:arg.reads.find(".")]}.fq.gz'
	if not os.path.exists(fastq):
		if arg.verbose: print('creating fastq file', file=sys.stderr)
		with gzip.open(fastq, 'wt') as fp:
			for name, seq in readfasta(arg.reads):
				print('@', name, file=fp, sep='')
				print(seq, file=fp)
				print('+', file=fp)
				print('J' * len(seq), file=fp)
	return fastq

def needfasta(arg):
	fasta = arg.reads[:-3]
	if not os.path.exists(fasta): os.system(f'gunzip -k {arg.reads}')
	return fasta

def sim4_to_ftx(filename, out, name):
	chrom = None
	strand = None
	exons = []
	ref = None
	n = 0
	with open(filename) as fp:
		for line in fp:
			if line.startswith('seq1 ='):
				if chrom is not None:
					ftx = FTX(chrom, str(n), strand, exons,
						f'{arg.program} ref:{ref}')
					print(ftx, file=out)
				chrom = None
				strand = None
				exons = []
				ref = line[7:].split(' ')[0][:-1]
				n += 1
				continue
			elif line.startswith('seq2 ='):
				f = line.split()
				chrom = f[2][:-1]
				continue
			f = line.split()
			if len(f) != 4: continue
			beg, end = f[1][1:-1].split('-')
			exons.append((int(beg) -1, int(end) -1))
			st = '+' if f[3] == '->' else '-'
			if strand is None: strand = st
			else: assert(strand == st)
		ftx = FTX(chrom, str(n), strand, exons, f'{name} ref:{ref}')
		print(ftx, file=out)

#######
# CLI #
#######

parser = argparse.ArgumentParser(description='spliced alignment evaluator: '
	+ ', '.join(('blat', 'bowtie2', 'bwa-mem', 'gmap', 'hisat2',
	'magicblast', 'minimap2', 'star', 'tophat2')))
parser.add_argument('genome', help='genome file in FASTA format')
parser.add_argument('reads', help='reads file in FASTA format')
parser.add_argument('program', help='program name')
parser.add_argument('--threads', type=int, default=1,
	help='number of threads if changeable [%(default)i]')
parser.add_argument('--optimize', action='store_true',
	help='run additional parameters rather than defaults')
parser.add_argument('--debug', action='store_true',
	help='keep temporary files (e.g. SAM)')
parser.add_argument('--verbose', action='store_true',
	help='show various informative messages')
arg = parser.parse_args()

###############
# Run Aligner #
###############

out = f'temp-{os.getpid()}' # temporary output file
ftx = f'ftx-{os.getpid()}'  # temporary ftx file
fp = open(ftx, 'w')

if arg.program == 'blat':
	# indexing: no
	# reads: fa.gz
	# output: sam not available, so using/converting from sim4
	# notes:
	#   no options for threads (there is a different pblat for that)
	#   adding -fine or -q=rna is suggested for long mRNA
	#   -maxIntron=N  defaults to 750000, changing had no effect
	#   -out=wublast blast blast8 blast9 sim4 maf axt pslx (using sim4)
	cli = f'blat {arg.genome} {arg.reads} {out} -out=sim4'
	if arg.optimize: cli += ' -fine -q=rna'
	if not arg.verbose: cli += ' > /dev/null'
	run(cli, arg)
	sim4_to_ftx(out, fp, arg.program)
elif arg.program == 'bowtie2':
	# indexing: yes
	# reads: fq.gz
	# output: sam
	# notes: does not do spliced alignment, so no point in optimizing
	if not os.path.exists(f'{arg.genome}.1.bt2'):
		cli = f'bowtie2-build {arg.genome} {arg.genome}'
		if not arg.verbose: cli += '>/dev/null 2> /dev/null'
		run(cli, arg)
	fastq = needfastq(arg)
	cli = f'bowtie2 -x {arg.genome} -U {fastq} > {out}'
	if not arg.verbose: cli += ' 2> /dev/null'
	run(cli, arg)
	sam_to_ftx(out, fp, arg.program)
elif arg.program == 'bwa-mem':
	# indexing: yes
	# reads: fa.gz or fq.gz
	# output: sam
	# notes: does not do spliced alignment, so no point in optimizing
	if not os.path.exists(f'{arg.genome}.bwt'):
		cli = f'bwa index {arg.genome}'
		if not arg.verbose: cli += ' 2> /dev/null'
		run(cli, arg)
	cli = f'bwa mem {arg.genome} {arg.reads} > {out}'
	if not arg.verbose: cli += ' 2> /dev/null'
	run(cli, arg)
	sam_to_ftx(out, fp, arg.program)
elif arg.program == 'gmap':
	# indexing: yes
	# reads: fa (not compressed)
	# output: sam
	# --microexon-spliceprob
	if not os.path.exists(f'{arg.genome}-gmap'):
		cli = f'gmap_build -d {arg.genome}-gmap -D . {arg.genome}'
		if not arg.verbose: cli += '>/dev/null 2>/dev/null'
		run(cli, arg)
	cli = f'gunzip -c {arg.reads} | gmap -d {arg.genome}-gmap -D . -f samse -t {arg.threads}'
	if arg.optimize: cli += f' -k 15 --min-intronlength 25 --microexon-spliceprob 0.05'
	cli += f' > {out}'
	if not arg.verbose: cli += ' 2>/dev/null'
	run(cli, arg)
	sam_to_ftx(out, fp, arg.program)
elif arg.program == 'hisat2':
	# indexing: yes
	# reads: fq.gz or fa.gz with -f
	# output: sam
	if not os.path.exists(f'{arg.genome}.1.ht2'):
		cli = f'hisat2-build -f {arg.genome} {arg.genome}'
		if not arg.verbose: cli += ' >/dev/null 2> /dev/null'
		run(cli, arg)
	cli = f'hisat2 -x {arg.genome} -U {arg.reads} -f -p {arg.threads}'
	if arg.optimize: cli += f' --min-intronlen 25'
	cli += f' > {out}'
	if not arg.verbose: cli += ' 2> /dev/null'
	run(cli, arg)
	sam_to_ftx(out, fp, arg.program)
elif arg.program == 'magicblast':
	# indexing: yes
	# reads: fa.gz
	# output: sam
	# notes: renames the chromosomes as numbers, maybe fix?
	if not os.path.exists(f'{arg.genome}.nsq'):
		cli = f'makeblastdb -dbtype nucl -in {arg.genome}'
		if not arg.verbose: cli += ' > /dev/null'
		run(cli, arg)
	cli = f'magicblast -db {arg.genome} -query {arg.reads} -num_threads {arg.threads}'
	if arg.optimize: cli += f' -word_size 15 -limit_lookup F'
	cli += f' > {out}'
	run(cli, arg)
	sam_to_ftx(out, fp, arg.program)
elif arg.program == 'minimap2':
	# indexing: no
	# reads: fa.gz or fq.gz
	# notes:
	#   -u can be f, b, n for matching transcript, both, or not match GT-AG
	cli = f'minimap2 -ax splice {arg.genome} {arg.reads} -t {arg.threads}'
	if arg.optimize: cli += f' -k 15 -u b'
	cli += f' > {out}'
	if not arg.verbose: cli += ' 2> /dev/null'
	run(cli, arg)
	sam_to_ftx(out, fp, arg.program)
elif arg.program == 'pblat':
	# see blat except reads need to be uncompressed (seekable)
	fasta = needfasta(arg)
	cli = f'pblat {arg.genome} {fasta} {out} -threads={arg.threads} -out=sim4'
	if arg.optimize: cli += ' -fine -q=rna'
	if not arg.verbose: cli += ' > /dev/null'
	run(cli, arg)
	sim4_to_ftx(out, fp, arg.program)
elif arg.program == 'star':
	# indexing: yes
	# reads: fa.gz or fq.gz
	# output: sam
	# notes
	#   twopassMode Basic ?
	#   requires temp file cleanup
	idx = f'{arg.genome}-star'
	if not os.path.exists(idx):
		cli = f'STAR --runMode genomeGenerate --genomeDir {idx}\
			--genomeFastaFiles {arg.genome} --genomeSAindexNbases 8'
		if not arg.verbose: cli += ' > /dev/null'
		run(cli, arg)
	cli = f'STAR --genomeDir {idx} --readFilesIn {arg.reads} --readFilesCommand "gunzip -c" --outFileNamePrefix {out} --runThreadN 1'
	if arg.optimize: cli += f' --alignIntronMin 25'
	if not arg.verbose: cli += ' > /dev/null'
	run(cli, arg)
	os.rename(f'{out}Aligned.out.sam', f'{out}')
	for x in ('Log.final.out', 'Log.out', 'Log.progress.out', 'SJ.out.tab'):
		os.unlink(f'{out}{x}')
	sam_to_ftx(out, fp, arg.program)
elif arg.program == 'tophat2':
	# indexing: yes, uses bowtie2
	# reads: fa.gz
	# output: sam
	# notes: must clean up tophat_out directory
	if not os.path.exists(f'{arg.genome}.1.bt2'):
		cli = f'bowtie2-build {arg.genome} {arg.genome}'
		if not arg.verbose: cli += '>/dev/null 2> /dev/null'
		run(cli, arg)
	opt = f'-i 25 --coverage-search --microexon-search --min-coverage-intron 25 --b2-very-sensitive' if arg.optimize else ''
	cli = f'tophat2 -p {arg.threads} {opt} {arg.genome} {arg.reads} '
	if not arg.verbose: cli += ' 2> /dev/null'
	run(cli, arg)
	cli = f'samtools view -h tophat_out/accepted_hits.bam > {out}'
	run(cli, arg)
	sam_to_ftx(out, fp, arg.program)
	if not arg.debug: os.system('rm -rf tophat_out')
else:
	sys.exit(f'ERROR: unknown program: {arg.program}')

fp.close()

#####################
# Report Alignments #
#####################

refs = [name for name, seq in readfasta(arg.reads)]
aligned = {}
with open(ftx) as fp:
	for line in fp:
		ali, ref = line.rstrip().split(' ref:')
		if ref not in aligned: aligned[ref] = ali
		# keeping only first match (blat sometimes has more than 1)

with gzip.open(f'{arg.program}.ftx.gz', 'wt') as fp:
	for ref in refs:
		if ref in aligned: print(ref, aligned[ref], sep='\t', file=fp)
		else:              print(ref, 'None', file=fp)

############
# Clean up #
############

if not arg.debug and os.path.exists(out): os.unlink(f'{out}')
if not arg.debug and os.path.exists(ftx): os.unlink(f'{ftx}')
