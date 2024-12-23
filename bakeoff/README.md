Bakeoff
=======

This is the home of the alignment bakeoff project.

- Stage 1: hello world
	- install a few aligners (done)
	- create minimal synthetic data (done)
	- determine cli for each program (in progress)
	- develop input/output wrapper for each program
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

```
python3 read-simulator.py ~/Code/datacore2024/genome_celegans/1* --double > ce.fq
```

For minimal testing purposes, an even smaller dataset is useful. The command
below generates 1000 reads covering 7 of the genes. There is a mixture of plus
and minus strand genes with variable intron sizes.

```
python3 read-simulator.py ~/Code/datacore2024/genome_celegans/1* --samplegenes 0.1 --samplereads 0.0321 --seed 1 --double > ce-min.fq
```

------------------------------------------------------------------------------

## Alignment Wrapper ##

All programs should result in SAM.

- The `--threads` parameter isn't implemented for most.
- Alignment and splicing parameters not optimized

- blat - alignments look suspect and no SAM
- bowtie2 - done
- bwa - done
- hisat2 -
- minimap2 - done
- star - done

DRAGEN needs evaluation: @dragen-1.hpc.genomecenter.ucdavis.edu
