import argparse
import gzip
import os
import re
import sys

import bakeoff

def run(cli, arg):
	cli = cli.replace('\t', ' ')
	if arg.verbose: print(cli, file=sys.stderr)
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

class BitFlag:
	def __init__(self, val):
		i = int(val)
		b = f'{i:012b}'
		self.read_reverse_strand = True if b[-5] == '1' else False
		self.not_primary_alignment = True if b[-9] == '1' else False
		self.supplementary_alignment = True if b[-12] == '1' else False
		self.unexpected_flag = False
		for i in (1, 2, 3, 4, 6, 7, 8, 10, 11):
			if b[-i] == '1': self.unexpected_flag = True
		self.binary = b

def cigar_to_exons(s, off):
	print(s, off)
	exons = []
	beg = off
	end = off
	for match in re.finditer(r'(\d+)([A-Z])', s):
		n = int(match.group(1))
		op = match.group(2)
		if   op == 'M': end += n
		elif op == 'D': pass
		elif op == 'I': end += n
		elif op == 'S': beg += n
		elif op == 'H': pass
		elif op == 'N':
			exons.append((beg, end))
			beg += n
			end = beg
	exons.append((beg, end))


	return exons


def sam_to_fx(filename):
	with open(filename) as fp:
		for line in fp:
			line = fp.readline()
			if line == '': break
			if line.startswith('@'): continue
			f = line.split('\t')
			bf = BitFlag(f[1])
			st = '-' if bf.read_reverse_strand else '+'
			# do something with supplementary, not primary, unexpected
			pos   = int(f[3])
			cigar = f[5]
			print(pos, cigar, st)
			#cigar_to_exons(f[5], int(f[3]))
			#print(pos, cigar, st)
			#stuff = cigar_to_exons(f[5], int(f[3]))
			#print(stuff)

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
	epilog='prog abbr: blat bt2 bwa gmap hisat magic mini star tophat')
parser.add_argument('genome', help='genome file in FASTA format')
parser.add_argument('reads', help='reads file in FASTA format')
parser.add_argument('prog', help='program name')
parser.add_argument('--threads', type=int, default=1,
	help='number of cores/threads [%(default)i]')
parser.add_argument('--debug', action='store_true', help='keep SAM file')
parser.add_argument('--verbose', action='store_true')
arg = parser.parse_args()

out = f'temp-{os.getpid()}'

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
	run(cli, arg)
	sim4_to_fx(out)
elif arg.prog == 'bt2':
	# indexing: yes
	# reads: fq.gz
	# output: sam
	# notes: does not do spliced alignment
	if not os.path.exists(f'{arg.genome}.1.bt2'):
		cli = f'bowtie2-build {arg.genome} {arg.genome}'
		if not arg.verbose: cli += '>/dev/null 2> /dev/null'
		run(cli, arg)
	needfastq(arg)
	cli = f'bowtie2 -x {arg.genome} -U {arg.reads} > {out}'
	if not arg.verbose: cli += ' 2> /dev/null'
	run(cli, arg)
	sam_to_fx(out)
elif arg.prog == 'bwa':
	# indexing: yes
	# reads: fa.gz or fq.gz
	# output: sam
	# notes: does not do spliced alignment
	if not os.path.exists(f'{arg.genome}.bwt'):
		cli = f'bwa index {arg.genome}'
		if not arg.verbose: cli += ' 2> /dev/null'
		run(cli, arg)
	cli = f'bwa mem {arg.genome} {arg.reads} > {out}'
	if not arg.verbose: cli += ' 2> /dev/null'
	run(cli, arg)
	sam_to_fx(out)
elif arg.prog == 'gmap':
	# indexing: yes
	# reads: fa (not compressed)
	# output:
	# --splice --microexon-spliceprob
	if not os.path.exists(f'{arg.genome}-gmap'):
		cli = f'gmap_build -d {arg.genome}-gmap -D . {arg.genome}'
		if not arg.verbose: cli += '>/dev/null 2>/dev/null'
		run(cli, arg)
	cli = f'gunzip -c {arg.reads} | gmap -d {arg.genome}-gmap -D . -f samse > {out}'
	if not arg.verbose: cli += ' 2>/dev/null'
	run(cli, arg)
	sam_to_fx(out)
elif arg.prog == 'hisat':
	# indexing: yes
	# reads: fq.gz
	# output: sam
	if not os.path.exists(f'{arg.genome}.1.ht2'):
		cli = f'hisat2-build -f {arg.genome} {arg.genome}'
		if not arg.verbose: cli += ' >/dev/null 2> /dev/null'
		run(cli, arg)
	needfastq(arg)
	cli = f'hisat2 -x {arg.genome} -U {arg.reads} > {out}'
	if not arg.verbose: cli += ' 2> /dev/null'
	run(cli, arg)
	sam_to_fx(out)
elif arg.prog == 'magic':
	# indexing: yes
	# reads: fa.gz
	# output: sam
	if not os.path.exists(f'{arg.genome}.nsq'):
		cli = f'makeblastdb -dbtype nucl -in {arg.genome}'
		if not arg.verbose: cli += ' > /dev/null'
		run(cli, arg)
	cli = f'magicblast -db {arg.genome} -query {arg.reads} > {out}'
	run(cli, arg)
	sam_to_fx(out)
elif arg.prog == 'mini':
	# indexing: no
	# reads: fa.gz or fq.gz
	cli = f'minimap2 -ax splice {arg.genome} {arg.reads} > {out}'
	if not arg.verbose: cli += ' 2> /dev/null'
	run(cli, arg)
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
		run(cli, arg)
	cli = f'STAR --genomeDir {idx}\
		--readFilesIn {arg.reads} --readFilesCommand "gunzip -c"\
		--runThreadN {arg.threads} --outFileNamePrefix {out}'
	if not arg.verbose: cli += ' > /dev/null'
	run(cli, arg)
	os.rename(f'{out}Aligned.out.sam', f'{out}')
	for x in ('Log.final.out', 'Log.out', 'Log.progress.out', 'SJ.out.tab'):
		os.unlink(f'{out}{x}')
	sam_to_fx(out)
elif arg.prog == 'tophat':
	# indexing: yes, bowttie2
	# reads: fq.gz probably
	# output: sam
	# notes: does not do spliced alignment
	if not os.path.exists(f'{arg.genome}.1.bt2'):
		cli = f'bowtie2-build {arg.genome} {arg.genome}'
		if not arg.verbose: cli += '>/dev/null 2> /dev/null'
		run(cli, arg)
	cli = f'tophat2 {arg.genome} {arg.reads}'
	if not arg.verbose: cli += ' 2> /dev/null'
	run(cli, arg)
	cli = f'samtools view -h tophat_out/accepted_hits.bam > {out}'
	run(cli, arg)
	sam_to_fx(out)
	if not arg.debug: os.system('rm -rf tophat_out')
else:
	sys.exit(f'ERROR: unknown program: {arg.prog}')

# Clean up

if not arg.debug and os.path.exists(out): os.unlink(f'{out}')
