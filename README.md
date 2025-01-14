intronomicon
============

Database of introns

## Manifest ##

In the order of how you use them

- `sra_xml_download.py` downloads xml files from NCBI for taxid
- `create_database.py` sets up sqlite db, reads sra xml, imports data
- `intron-injector.py` requests run, aligns, post-processes, injects introns
- `meta-injector.py` requests run, parses SRA/GEO, injects meta
- `sraxml.py` library for interacting with sra xml files
- `xml_reader.py` toy for figuring shit out
- `notes` random thoughts of the moment
