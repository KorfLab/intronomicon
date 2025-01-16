intronomicon
============

Database of introns

## Manifest ##

In the order of how you use them

- `downloader.py` downloads SRA and GEO files from NCBI for taxid
- `create-database.py` sets up sqlite db, reads SRA/GEO, imports data
- `intron-injector.py` requests run, aligns, post-processes, injects introns
- `meta-injector.py` requests run, parses SRA/GEO, injects meta
- `meta-tabler.py` creates some training data
- `sraxml.py` library for interacting with sra xml files
- `xmltoy.py` Ian's toy for figuring shit out
- `notes` random thoughts of the moment
