import argparse
import gzip
import os
import sys


def run(cli, verbose):
	cli = cli.replace('\t', ' ')
	if verbose: print(cli, file=sys.stderr)
	os.system(cli)

# CLI

parser = argparse.ArgumentParser(description='read alignment wrapper',
	epilog='programs: bt2 bwa hs2 mm2 star')
parser.add_argument('genome', help='genome file in FASTA format')
parser.add_argument('reads', help='reads file in fq.gz format')
parser.add_argument('prog', help='program name')
parser.add_argument('--threads', type=int, default=1,
	help='number of cores/threads [%(default)i]')
parser.add_argument('--debug', action='store_true', help='keep temp file')
parser.add_argument('--verbose', action='store_true')
arg = parser.parse_args()

# Run the aligner

out = f'temp-{os.getpid()}'

if arg.prog == 'bt2':
	# index
	if not os.path.exists(f'{arg.genome}.1.bt2'):
		cli = f'bowtie2-build {arg.genome} {arg.genome}'
		if not arg.verbose: cli += '2> /dev/null'
		run(cli, arg.verbose)
	# align, no clean
	cli = f'bowtie2 -x {arg.genome} -U {arg.reads} > {out}.sam'
	if not arg.verbose: cli += '2> /dev/null'
	run(cli, arg.verbose)
elif arg.prog == 'bwa':
	# index
	if not os.path.exists(f'{arg.genome}.bwt'):
		cli = f'bwa index {arg.genome}'
		if not arg.verbose: cli += '2> /dev/null'
		run(cli, arg.verbose)
	# align, no clean
	cli = f'bwa mem {arg.genome} {arg.reads} > {out}.sam'
	if not arg.verbose: cli += '2> /dev/null'
	run(cli, arg.verbose)
elif arg.prog = 'hs2':
	sys.exit('not ready yet')
elif arg.prog == 'mm2':
	# no index, no clean
	cli = f'minimap2 -ax splice {arg.genome} {arg.reads} > {out}.sam'
	run(cli, arg.verbose)
elif arg.prog == 'star':
	# index
	idx = f'{arg.genome}-star'
	if not os.path.exists(idx):
		cli = f'STAR --runMode genomeGenerate --genomeDir {idx}\
			--genomeFastaFiles {arg.genome} --genomeSAindexNbases 8\
			--runThreadN {arg.threads}'
		if not arg.verbose: cli += '> /dev/null'
		run(cli, arg.verbose)
	# align
	cli = f'STAR --genomeDir {idx}\
		--readFilesIn {arg.reads} --readFilesCommand "gunzip -c"\
		--runThreadN {arg.threads} --outFileNamePrefix {out}'
	if not arg.verbose: cli += '> /dev/null'
	run(cli, arg.verbose)
	# 2nd pass?
	# clean
	os.rename(f'{out}Aligned.out.sam', f'{out}.sam')
	for x in ('Log.final.out', 'Log.out', 'Log.progress.out', 'SJ.out.tab'):
		os.unlink(f'{out}{x}')


# Convert the SAM to something else...

# Clean up

#if not arg.debug: os.unlink(out)
