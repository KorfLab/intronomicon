Bakeoff
=======

## WARNING ##

EVERYTHING IS CURRENTLY BROKEN. THERE IS A NEW BAKEOFF LIBRARY. MUST CHANGE ALL
PREVIOUS CODE TO USE LIBRARY.


This is the home of the alignment bakeoff project.

- Stage 1: hello world
	- install a few aligners
	- create minimal synthetic data
	- determine cli for each program
	- develop input wrapper for each program
	- develop output wrapper for each program
	- develop comparison program to check accuracy
- Stage 2: parameter optimization
	- figure out best way to run each program
	- create synthetic data with adversarial properties
		- C. elegans sampled
		- synthetic genome
- Stage 3: more programs
	- install other aligners
	- add input/output wrappers
	- optimize parameters
- Stage 4: best practices
	- multiple stages?

## Flattened Transcript Format ##

Oh yay, a new non-standard standard. This file format is used for the bakeoff
project only.

- file extension: `.fx` (not an official file extension)
- field delimiter: `|`
- 4 or 5 fields

1. chromosome identifier, should match some fasta file
2. name of transcript, should be unique for each line
3. strand indicator `+` or `-`
4. exon structure:
	- hyphen separated coordinates
	- comma separated exons
	- must be sorted left to right, low to high
5. information: optional free text

Example: Plus-strand transcript with introns inferred at 201-299 and 401-499.

```
chr1|name|+|100-200,300-400,500-600
chr1|name|+|100-200,300-400,500-600|with extra stuff
```

## Install aligners ##

Created a simplified conda environment `intronomicon-align.yml` for 6 aligners:

- blat
- bowtie2
- bwa
- hisat2
- minimap2
- star

There are several other aligners to test, but 6 is plenty to get the comparison
framework underway.

## Minimal synthetic data ##

`read-simulator.py` is used to create a synthetic set of reads with complete
coverage for each canonical transcript of every protein-coding gene in the 1pct
sets in datacore2024.

For minimal testing purposes, an even smaller dataset is useful. Use
`--samplegenes 0.1` and `--samplereads 0.1` to reduce the set size

```
python3 read-simulator.py  ~/Code/datacore2024/genome_celegans/1* --double  --samplegenes 0.1 --samplereads 0.1 --seed 29 > mini.fa
```

The output of `read-simulator.py` has headers that look like the following:

```
>81|Transcript:C23H3.3b.1|+|57750-57846,58307-58309|-
```

- "81" means the read was generated at position 81 of the transcript
- "Transcript:C23H3.3b.1" is the name of the transcript
- "+" means the transcript is on the plus strand
- "57750-57846,58307-58309" describes the position of the read in the genome
- "-" means the sequence is generated from the reverse-complement

------------------------------------------------------------------------------

## Alignment Wrapper ##

Need to modify this for fasta input as that is now what I'm using

All programs should result in SAM and are then processed to _something_ else.

- The `--threads` parameter isn't implemented for most.
- Alignment and splicing parameters not optimized

- blat - alignments look suspect and no SAM
- bowtie2 - done
- bwa - done
- hisat2 - done
- minimap2 - done
- star - done

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
- No non-canonical splice sites yet
- Reads should also be varied:
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