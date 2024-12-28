import argparse
import bakeoff
import json

parser = argparse.ArgumentParser()
parser.add_argument('genome', help='genome in fasta format')
parser.add_argument('exons', help='reference annotation in fx format')
parser.add_argument('reads', help='read alignments in fx format')
parser.add_argument('--introns', metavar='<file>',
	help='save introns to <file>')
arg = parser.parse_args()

real = bakeoff.fxread(arg.exons)
test = bakeoff.fxread(arg.reads)

print(json.dumps(test, indent=4))
