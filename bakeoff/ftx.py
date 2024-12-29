# ftx.py
import sys

class FTX:
	"""Class to represent transcripts with one-line formatting"""

	def __init__(self, chrom, name, strand, exons, info):
		# sanity checks
		assert(' ' not in chrom)
		assert(' ' not in name)
		assert(strand == '+' or strand == '-')
		for beg, end in exons: assert(beg <= end)
		for i in range(len(exons) -1): assert(exons[i][0] < exons[i+1][0])

		self.chrom = chrom
		self.beg = exons[0][0]
		self.end = exons[-1][1]
		self.name = name
		self.strand = strand
		self.exons = exons
		self.info = info

	def overlaps(f1, f2, strand_sensitive=True):
		if f1.chrom != f2.chrom: return False
		if f1.strand != f2.strand and strand_sensitive: return False
		if f1.beg <= f2.beg and f1.end >= f2.beg: return True
		return False

	def text(self):
		"""text-based version of ftx, 1-based"""
		estr = ','.join([f'{beg+1}-{end+1}' for beg, end in self.exons])
		return '|'.join((self.chrom, self.name, self.strand, estr, self.info))

	def __str__(self):
		return self.text()

	@classmethod
	def parse(self, text):
		"""returns ftx object from ftx string"""
		chrom, name, strand, estr, info = text.split('|')
		exons = []
		for s in estr.split(','):
			beg, end = s.split('-')
			exons.append((int(beg)-1, int(end)-1))
		return FTX(chrom, name, strand, exons, info)
