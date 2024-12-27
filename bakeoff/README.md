Bakeoff
=======

This is the home of the alignment bakeoff project.

- Stage 1: hello world
	- install a few aligners (done)
	- create minimal synthetic data (done)
	- determine cli for each program (done)
	- develop input wrapper for each program (done)
	- all programs result in SAM output (done except blat)
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
