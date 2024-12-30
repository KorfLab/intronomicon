import argparse
import gzip
import os
import re
import sys

from ftx import FTX

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
	fastq = f'{arg.reads}.fq.gz'
	if not os.path.exists(fastq):
		if arg.verbose: print('creating fastq file', file=sys.stderr)
		with gzip.open(fastq, 'wt') as fp:
			for name, seq in readfasta(arg.reads):
				print('@', name, file=fp, sep='')
				print(seq, file=fp)
				print('+', file=fp)
				print('J' * len(seq), file=fp)
	return fastq

class BitFlag:
	def __init__(self, val):
		i = int(val)
		b = f'{i:012b}'
		self.read_unmapped = True if b[-3] == '1' else False
		self.read_reverse_strand = True if b[-5] == '1' else False
		self.not_primary_alignment = True if b[-9] == '1' else False
		self.supplementary_alignment = True if b[-12] == '1' else False
		self.otherflags = []
		for i in (1, 2, 4, 6, 7, 8, 10, 11):
			if b[-i] == '1': self.otherflags.append(i)

def cigar_to_exons(cigar, pos):
	exons = []
	beg = 0
	end = 0
	for match in re.finditer(r'(\d+)([A-Z])', cigar):
		n = int(match.group(1))
		op = match.group(2)
		if   op == 'M': end += n
		elif op == 'D': pass
		elif op == 'I': end += n
		#elif op == 'S': beg += n # is this right?
		#elif op == 'H': pass
		elif op == 'N':
			exons.append((pos+beg-1, pos+end-2))
			beg = end + n
			end = beg
	exons.append((pos+beg-1, pos+end-2))
	return exons

def sam_to_ftx(filename, out, arg):
	n = 0
	with open(filename) as fp:
		for line in fp:
			if line == '': break
			if line.startswith('@'): continue
			f = line.split('\t')
			qname = f[0]
			bf = BitFlag(f[1])
			chrom = f[2]
			pos   = int(f[3])
			cigar = f[5]

			st = '-' if bf.read_reverse_strand else '+'
			if bf.read_unmapped: continue
			if bf.not_primary_alignment: continue
			if bf.supplementary_alignment: continue
			if bf.otherflags:
				print(bf.otherflags)
				sys.exit('unexpected flags found, debug time')
			n += 1
			exons = cigar_to_exons(cigar, pos)
			ftx = FTX(chrom, str(n), st, exons, f'{arg.program} ref:{qname}')
			print(ftx, file=out)

def sim4_to_ftx(filename, out, arg):
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
		ftx = FTX(chrom, str(n), strand, exons, f'{arg.program} ref:{ref}')
		print(ftx, file=out)

#######
# CLI #
#######

parser = argparse.ArgumentParser(description='spliced alignment evaluator: '
	+ ', '.join(('blat', 'bowtie2', 'bwa-mem', 'gmap', 'hisat2',
	'magic-blast', 'minimap2', 'star', 'tophat2')))
parser.add_argument('genome', help='genome file in FASTA format')
parser.add_argument('reads', help='reads file in FASTA format')
parser.add_argument('program', help='program name')
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
	# output: sim4 converted to fx
	# notes:
	#   consider adding -fine -q=rna (they are used for full length mRNAs)
	#   -maxIntron=N  defaults to 750000
	#   -out=wublast blast blast8 blast9 sim4 maf axt pslx (using sim4)
	cli = f'blat {arg.genome} {arg.reads} {out} -out=sim4'
	if not arg.verbose: cli += ' > /dev/null'
	run(cli, arg)
	sim4_to_ftx(out, fp, arg)
elif arg.program == 'bowtie2':
	# indexing: yes
	# reads: fq.gz
	# output: sam
	# notes: does not do spliced alignment
	if not os.path.exists(f'{arg.genome}.1.bt2'):
		cli = f'bowtie2-build {arg.genome} {arg.genome}'
		if not arg.verbose: cli += '>/dev/null 2> /dev/null'
		run(cli, arg)
	fastq = needfastq(arg)
	cli = f'bowtie2 -x {arg.genome} -U {fastq} > {out}'
	if not arg.verbose: cli += ' 2> /dev/null'
	run(cli, arg)
	sam_to_ftx(out, fp, arg)
elif arg.program == 'bwa-mem':
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
	sam_to_ftx(out, fp, arg)
elif arg.program == 'gmap':
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
	sam_to_ftx(out, fp, arg)
elif arg.program == 'hisat2':
	# indexing: yes
	# reads: fq.gz
	# output: sam
	if not os.path.exists(f'{arg.genome}.1.ht2'):
		cli = f'hisat2-build -f {arg.genome} {arg.genome}'
		if not arg.verbose: cli += ' >/dev/null 2> /dev/null'
		run(cli, arg)
	fastq = needfastq(arg)
	cli = f'hisat2 -x {arg.genome} -U {fastq} > {out}'
	if not arg.verbose: cli += ' 2> /dev/null'
	run(cli, arg)
	sam_to_ftx(out, fp, arg)
elif arg.program == 'magic-blast':
	# indexing: yes
	# reads: fa.gz
	# output: sam
	if not os.path.exists(f'{arg.genome}.nsq'):
		cli = f'makeblastdb -dbtype nucl -in {arg.genome}'
		if not arg.verbose: cli += ' > /dev/null'
		run(cli, arg)
	cli = f'magicblast -db {arg.genome} -query {arg.reads} > {out}'
	run(cli, arg)
	sam_to_ftx(out, fp, arg)
elif arg.program == 'minimap2':
	# indexing: no
	# reads: fa.gz or fq.gz
	cli = f'minimap2 -ax splice {arg.genome} {arg.reads} > {out}'
	if not arg.verbose: cli += ' 2> /dev/null'
	run(cli, arg)
	sam_to_ftx(out, fp, arg)
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
	cli = f'STAR --genomeDir {idx}\
		--readFilesIn {arg.reads} --readFilesCommand "gunzip -c"\
		--outFileNamePrefix {out}'
	if not arg.verbose: cli += ' > /dev/null'
	run(cli, arg)
	os.rename(f'{out}Aligned.out.sam', f'{out}')
	for x in ('Log.final.out', 'Log.out', 'Log.progress.out', 'SJ.out.tab'):
		os.unlink(f'{out}{x}')
	sam_to_ftx(out, fp, arg)
elif arg.program == 'tophat2':
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
	sam_to_ftx(out, fp, arg)
	if not arg.debug: os.system('rm -rf tophat_out')
else:
	sys.exit(f'ERROR: unknown program: {arg.program}')

fp.close()

######################
# Evalute Alignments #
######################

hitcount = {name:0 for name, seq in readfasta(arg.reads)}

with open(ftx) as fp:
	for line in fp:
		f1, f2 = line.rstrip().split(' ref:')
		alignment = FTX.parse(f1)
		reference = FTX.parse(f2)
		shared, total = reference.compare_coordinates(alignment)
	#	print(reference)
	#	print(alignment)
		print(shared, total)

		hitcount[f'{reference}'] += 1

missed = 0
extra = 0
unique = 0
for name, n in hitcount.items():
	if   n == 0: missed += 1
	elif n > 1:  extra += 1
	else:        unique += 1

print(f'missed: {missed}, extra: {extra}, unique: {unique}')



############
# Clean up #
############

if not arg.debug and os.path.exists(out): os.unlink(f'{out}')
if not arg.debug and os.path.exists(ftx): os.unlink(f'{ftx}')
