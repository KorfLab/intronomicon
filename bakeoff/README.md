Bakeoff
=======

This is the home of the alignment bakeoff sub-project that seeks to answer the
questions:

- Which is the most accurate aligner for RNA-seq data?
- What kinds of exon-intron structures create the most problems?
- What kinds of improvements need to made in this space?"

This project is moving to its own repo shortly.

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

might want to use a masked genome here
and omit all transcripts with Ns in the exons

## Synthetic Experiment ##

Ideally, there should be a wrapper for the entire experiment that collects the
resource usage as well.


### Setup

```
python3 genome-simulator.py exp1 --seed 1 --double --noncanonical
./bakeoff -x exp1.fa exp1.gff exp1
```

- genes: 3680
- reads: 1,479,360
- bases: 147,936,000

### blat

blat doesn't do separate indexing and only uses 1 cpu

```
/usr/bin/time -f "%P %M %U %S %E" ./bakeoff -x exp1.fa exp1.gff exp1 blat
101% 512952 46.17 2.84 0:48.53
```

### bowtie2

1 cpu with indexing

```
/usr/bin/time -f "%P %M %U %S %E" ./bakeoff -x exp1.fa exp1.gff exp1 bowtie2
98% 285928 62.47 1.84 1:05.52
```

2 cpus with previous indexing

```
/usr/bin/time -f "%P %M %U %S %E" ./bakeoff -xp2 exp1.fa exp1.gff exp1 bowtie2
98% 285820 53.88 1.39 0:56.34
```

### bwa-mem

```
/usr/bin/time -f "%P %M %U %S %E" ./bakeoff -x exp1.fa exp1.gff exp1 bwa-mem
95% 495716 38.37 2.19 0:42.41
```

2 cpus with previous indexing

```
/usr/bin/time -f "%P %M %U %S %E" ./bakeoff -xp2 exp1.fa exp1.gff exp1 bwa-mem
99% 495832 39.76 2.04 0:41.82
```

### gmap

```
/usr/bin/time -f "%P %M %U %S %E" ./bakeoff -x exp1.fa exp1.gff exp1 gmap
100% 507772 6720.62 721.45 2:03:18
```

2 cpus with previous indexing

```
/usr/bin/time -f "%P %M %U %S %E" ./bakeoff -xp2 exp1.fa exp1.gff exp1 gmap
200% 507800 7238.28 1053.59 1:08:48
```

### hisat2

```
usr/bin/time -f "%P %M %U %S %E" ./bakeoff -x exp1.fa exp1.gff exp1 hisat2
97% 453324 25.57 1.69 0:27.86
```

2 cpus with previous indexing

```
/usr/bin/time -f "%P %M %U %S %E" ./bakeoff -xp2 exp1.fa exp1.gff exp1 hisat2
131% 444808 28.49 4.70 0:25.26
```

### magicblast

```
/usr/bin/time -f "%P %M %U %S %E" ./bakeoff -x exp1.fa exp1.gff exp1 magicblast
98% 3838932 48.34 8.08 0:57.04
```

2 cpus with previous indexing

```
/usr/bin/time -f "%P %M %U %S %E" ./bakeoff -xp2 exp1.fa exp1.gff exp1 magicblast
139% 7063188 52.20 20.25 0:51.99
```

### minimap2

```
/usr/bin/time -f "%P %M %U %S %E" ./bakeoff -x exp1.fa exp1.gff exp1 minimap2
99% 632196 39.52 1.75 0:41.29
```

2 cpus with previous indexing

```
/usr/bin/time -f "%P %M %U %S %E" ./bakeoff -xp2 exp1.fa exp1.gff exp1 minimap2
138% 632196 32.42 8.94 0:29.86
```

### pblat

```
/usr/bin/time -f "%P %M %U %S %E" ./bakeoff -x exp1.fa exp1.gff exp1 pblat
99% 514032 45.01 1.03 0:46.04
```

2 cpus with previous indexing

```
/usr/bin/time -f "%P %M %U %S %E" ./bakeoff -xp2 exp1.fa exp1.gff exp1 pblat
148% 514184 43.18 0.90 0:29.76
```


### star

```
/usr/bin/time -f "%P %M %U %S %E" ./bakeoff -x exp1.fa exp1.gff exp1 star
99% 511572 31.46 1.20 0:32.92
```

2 cpus with previous indexing

```
/usr/bin/time -f "%P %M %U %S %E" ./bakeoff -xp2 exp1.fa exp1.gff exp1 star
99% 511528 31.02 0.96 0:32.23
```

### tophat2

```
/usr/bin/time -f "%P %M %U %S %E" ./bakeoff -x exp1.fa exp1.gff exp1 tophat2
146% 442372 222.51 233.12 5:10.03
```

2 cpus with previous indexing

```
/usr/bin/time -f "%P %M %U %S %E" ./bakeoff -xp2 exp1.fa exp1.gff exp1 tophat2
177% 880140 217.00 246.21 4:20.46
```


