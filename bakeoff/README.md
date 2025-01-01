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

## Programs ##

- `genome-simulator.py` creates an experimental genome and annotation
- `read-simulator.py` creates synthetic RNA-seq reads from FASTA + GFF
- `run-aligner.py` provides a consistent interface to multiple programs
- `compare-alignments.py` evaluates performance of aligners

## Progress ##

- Stage 1: hello world
	- install spliced aligners (done)
	- create minimal synthetic data (done)
	- determine cli for each program (done)
	- develop input wrapper for each program (done)
	- develop output wrapper for each program (done)
	- develop comparison program to check accuracy (done for now)
-Stage 2: synthetic data with adversarial properties
	- create synthetic genomes with
		- variable sized exons
		- non-canonical splice sites
	- figure out best way to run each program
- Stage 3: full study
	- create full C. elegans and synthetic data sets
	- perform alignments
	- compare performance and resources
	- determine best practices
	- determine areas to improve
