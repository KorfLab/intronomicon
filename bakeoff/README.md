Bakeoff
=======

This is the home of the alignment bakeoff sub-project that seeks to answer the
questions "which is the best aligner for RNA-seq data?" and "what kinds of
improvements need to made in this space?"

## Quickstart ##

1. Install conda (e.g. miniforge3)
2. Clone the intronomicon repo
3. Create the bakeoff environment `conda env create -f bakeoff.yml`
4. Run the `bakeoff` program (see usage)

Notes:

- not all software built for ARM, may need an alternate environment
- need to do testing on non-x86 hardware
- dragen needs evaluation outside this scope

## Flattened Transcript Format ##

Oh yay, a new non-standard standard. This file format is used internally for
the bakeoff project only.

- file extension: `.ftx` (not an official file extension)
- field delimiter: `|`
- 5 fields

1. chromosome identifier, should match some fasta file
2. name of transcript, should be unique for each line
3. strand indicator `+` or `-`
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

## Minimal synthetic data ##

`read-simulator.py` is used to create a synthetic reads across the entirety of
a gene. The FASTA header describes the coordinates of the reads in ftx format.

For testing purposes, something like this is useful.

```
python3 read-simulator.py  ~/Code/datacore2024/genome_celegans/1* --double  --samplegenes 0.1 --samplereads 0.1 --seed 13 > mini.fa
```

## Alignment Runner ##

`run-aligner.py` provides a consistent interface for every program. Many
things are not optimized yet.

- Resources
- Alignment and splicing parameters

Wrappers for the following exist so far.

- blat
- bowtie2
- bwa
- gmap
- hisat2
- magicblast
- minimap2
- star
- tophat

## Genome Simulator ##

The `genome-simulator.py` program creates synthetic genomes and corresponding
annotation.


## Progress ##

- Stage 1: hello world
	- install spliced aligners (done)
	- create minimal synthetic data (done)
	- determine cli for each program (done)
	- develop input wrapper for each program (done)
	- develop output wrapper for each program (done)
	- develop comparison program to check accuracy (in progress)
		- **may require debugging previous steps**
	- figure out best way to run each program
		- accuracy
		- time, memory, cpus
-Stage 2: synthetic data with adversarial properties
	- create synthetic genomes
	- test alignments
- Stage 3: full study
	- create full C. elegans and synthetic data sets
	- perform alignments
	- compare performance and resources
	- determine best practices
	- determine areas to improve
