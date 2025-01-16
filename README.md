intronomicon
============

Database of introns

## Manifest ##

In the order of how you use them

- `downloader.py` downloads SRA and GEO files from NCBI for taxid
- `create_database.py` sets up sqlite db, reads SRA/GEO, imports data
- `intron-injector.py` requests run, aligns, post-processes, injects introns
- `meta-injector.py` requests run, parses SRA/GEO, injects meta
- `sraxml.py` library for interacting with sra xml files
- `xml_reader.py` toy for figuring shit out
- `notes` random thoughts of the moment
