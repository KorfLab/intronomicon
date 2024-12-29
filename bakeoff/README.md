Bakeoff
=======

This is the home of the alignment bakeoff sub-project that seeks to answer the
questions "which is the best aligner of RNA-seq data?" and "what kinds of
improvements need to made in this space?"

- Stage 1: hello world
	- install spliced aligners (done)
	- create minimal synthetic data (done)
	- determine cli for each program (done)
	- develop input wrapper for each program (done)
	- develop output wrapper for each program (in progress for sam)
	- develop comparison program to check accuracy
- Stage 2: parameter optimization
	- figure out best way to run each program
		- accuracy
		- time, memory, cpus
	- create synthetic data with adversarial properties
		- C. elegans sampled
		- synthetic genome
- Stage 3: full study
	- create full data sets
	- perform alignments
	- compare performance
	- determine best practices
	- determine areas to improve

## Environment ##

```
cd bakeoff
conda env create -f bakeoff.yml
conda activate bakeoff
```

## Flattened Transcript Format ##

Oh yay, a new non-standard standard. This file format is used for the bakeoff
project only.

- file extension: `.fx` (not an official file extension)
- field delimiter: `|`
- 5 fields

1. chromosome identifier, should match some fasta file
2. name of transcript, should be unique for each line
3. strand indicator `+` or `-`
4. exon structure:
	- hyphen separated coordinates
	- comma separated exons
	- must be sorted left to right, low to high
5. information: extra free text

Example: Plus-strand transcript with introns inferred at 201-299 and 401-499.

```
chr1|name|+|100-200,300-400,500-600|whatever you like
```

## Alignment Programs ##

Created a conda environment `bakeoff.yml` for several aligners:

- blat
- gmap
- hisat2
- magicblast
- minimap2
- star
- tophat

Originally, bowtie2 and bwa were included, but they do not perform spliced
alignment, so they are not going to be tested. They are still included in the
alignment wrapper (see below).

Aligners not yet considered here.

- dragen

## Minimal synthetic data ##

`read-simulator.py` is used to create a synthetic set of reads with complete
coverage for each canonical transcript of every protein-coding gene in the 1pct
sets in datacore2024.

For minimal testing purposes, an even smaller dataset is useful. Use
`--samplegenes 0.1` and `--samplereads 0.1` to reduce the set size

```
python3 read-simulator.py  ~/Code/datacore2024/genome_celegans/1* --double  --samplegenes 0.1 --samplereads 0.1 --seed 13 > mini.fa
```

The output of `read-simulator.py` has headers that look like the following:

```
>I|Transcript:F53G12.5b.1|+|127305-127336,127385-127436,128697-128712|-
```

- "I" means it was generated from chromosome I
- "Transcript:F53G12.5b.1" is the name of the transcript
- "+" means the transcript is on the plus strand
- "27305-127336,127385-127436,128697-128712" are genome coordinates
- "-" means the sequence is generated from the reverse-complement

------------------------------------------------------------------------------

## Alignment Wrapper ##

`alignment-wrapper.py` provides a consistent interface for every program. Input
files are fasta, output is .fx.

Output is not yet complete for SAM files.

Many things are not optimized yet.

- Resources
- Alignment and splicing parameters

Wrappers for the following

- blat
- bowtie2
- bwa
- gmap
- hisat2
- magicblast
- minimap2
- star
- tophat

DRAGEN needs evaluation: @dragen-1.hpc.genomecenter.ucdavis.edu

------------------------------------------------------------------------------

## Genome Simulator ##

The `genome-simulator.py` program creates synthetic genomes and corresponding
annotation.

- Every gene is flanked on either side by a fixed amount of sequence (100)
- Every gene has 5 exons
- The standard exon length is 100
- The standard intron length is 100
- Intron 2 is variable length
- Exon 4 is variable length
- `f1-[e1]-i1-[e2]-v.i2-[e3]-i3-[v.e4]-i4-[e5]-f2`

Todo...

- Genes are currently only on the positive strand
	- meaning introns are all GT..AG
	- there is no coding sequence implied
- No non-canonical splice sites yet
- There should be a companion read mutator
	- substitution
	- indels
	- trans-spliced leaders
	- poly-A tails




## Comparison ##

Comparing truth to aligner

- we know the location of the read, it's in the fx
- aligners report locations of the read
- what if the alignment reports multiple locations?
	- reporting too many should lower aggregate performance

correct - both sites found
partial - one site found
missed - no sites found
fabricated - extra sites found


when a program is wrong, why did it get the wrong answer?


pie-charts for each algorithm?
