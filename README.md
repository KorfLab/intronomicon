INTRONOMICON
============

The Book of Introns

## Quick Start ##

Download the files. If you run this again later, it skips over files already
present.

```
python3 downloader build --verbose
```

Build the database for worm introns.

```
python3 database-builder.py worm.db build
```




## Manifest ##

In the order of use (-ish)

- `downloader.py` retrieve GEO and SRA files for a taxid
- `database-builder.py` read GEO and SRA files into database
- `intron-injector.py` add introns to database (via sequence alignment)
- `ncbi_reader.py` shared library for reading various files from NCBI
- `cultist.py` something for debugging and exploring database



## Notes ##

- completely new downloader.py
- database-builder.py needs to be refactored to fit
- database will change substantially
	- might not store meta-data
	- more intron focused
