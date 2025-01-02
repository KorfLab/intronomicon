Bakeoff
=======

This is the home of the alignment bakeoff sub-project that seeks to answer the
questions:

- Which is the most accurate aligner for RNA-seq data?
- What kinds of exon-intron structures create the most problems?
- What kinds of improvements need to made in this space?"

## Quickstart ##

1. Install conda (e.g. miniforge3)
2. Clone the intronomicon repo
3. Create the bakeoff environment `conda env create -f bakeoff.yml`
4. Run the `bakeoff` program (see usage)

Notes:

- you will need `grimore` to simulate reads
- you may want `datacore2024` for minimal datasets
- x86 (only hisat2, minimap2, and pblat on osx-arm)
- dragen needs evaluation outside this scope

## Flattened Transcript Format ##

Oh yay, a new non-standard standard. This file format is used internally for
the bakeoff project only. The reason this is used is to embed the meta-data for
alignment coordinates directly in the FASTA/FASTQ deflines.

- file extension: `.ftx` (not an official file extension)
- field delimiter: `|`
- 5 fields

1. chromosome identifier, should match some fasta file
2. name of transcript, should be unique for each line
3. strand indicator `+` or `-` for the transcript
4. exon structure:
	- hyphen separated coordinates
	- comma separated exons
	- must be sorted left to right, low to high
	- numbers are 1-based
5. information: extra free text

Example: Plus-strand transcript with introns inferred at 201-299 and 401-499.

```
chr1|name|+|100-200,300-400,500-600|whatever you like
```

## Programs & Libraries ##

- `genome-simulator.py` creates an experimental genome and annotation
- `read-simulator.py` creates synthetic RNA-seq reads from FASTA + GFF
- `run-aligner.py` provides a consistent interface to multiple programs
- `compare-alignments.py` evaluates performance of aligners
- `ftx.py` defines the ftx class and methods
- `sam2ftx.py` does what it says and is also a library for `run-aligner.py`

## Progress ##

- Stage 1: hello world
	- install spliced aligners (done)
	- create minimal synthetic data (done)
	- determine cli for each program (done)
	- develop input wrapper for each program (done)
	- develop output wrapper for each program (done)
	- develop comparison program to check accuracy (done for now)
-Stage 2: synthetic data with adversarial properties
	- create synthetic genomes with variable exon and splice sites (done)
	- figure out best way to run each program (maybe)
	- full reads on 1% data
- Stage 3: full study
	- create full C. elegans and synthetic data sets
	- perform alignments
	- compare performance and resources
	- determine best practices
	- determine areas to improve

## 1pct_full Experiment ##

- genes: 131
- reads: 357,568
- bases: 35,756,800
- runtime (approximate, including indexing and 2 cpus):
	- hisat2: 5.5
	- star: 6
	- pblat: 7
	- bwa-mem 8.5
	- blat: 9
	- minimap2: 10
	- bowtie2 12
	- magicblast: 22
	- tophat2: 26
	- gmap: 490 (no, that's not a typo)

sort of like 2 min (all except gmap) and 8 min (gmap) for the 1% genome, so
expecting 200 min and 800 min for full genome, or ~17 hours.

## Synthetic Experiment ##

- genes: 3680
- reads: 1,479,360
- bases: 147,936,000

The synthetic experiment with 10 chromosomes is about 4 times larger than the
1pct_full (8 + 32 min)

magicblast uses more than 4G memory and kills itself on the vm
