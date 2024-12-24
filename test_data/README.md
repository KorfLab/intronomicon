Test Data
=========



## C. elegans ##

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

## Run Examples ##

Many tutorials at https://docs.tinybio.cloud/docs/tinybio-intern-user-guide


blat

```
blat 1pct.fa ce-min.fa foo
```

bowtie2

```
bowtie2-build 1pct.fa 1pct
bowtie2 -x 1pct -U ce-min.fq.gz > foo
```

bwa

```
bwa index 1pct.fa
bwa mem 1pct.fa ce-min.fq.gz > foo
```

minimap2

```
minimap2 -ax splice 1pct.fa ce-min.fq.gz > foo
```

Others

- star
- hisat2
- rsem
- salmon
- kallisto
- tpmcalculator

@dragen-1.hpc.genomecenter.ucdavis.edu
