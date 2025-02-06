INTRONOMICON
============

The Book of Introns

## Manifest ##

In the order of use (-ish)

- `downloader.py` retrieve GEO and SRA files for a taxid
- `database-builder.py` read GEO and SRA files into database
- `intron-injector.py` add introns to database (via sequence alignment)
- `sraxml.py` shared library for reading XML files from SRA
- `cultist.py` something for debugging and exploring database
- `xmltoy.py` something for debugging and exploring xml files
- `notes.md` random stuff the devs are currently thinking about


## Notes ##

- Don't go after the raw files
- Need to use the FTP site
	- Check regularity of the files there
- Refactored the soft2text and this breaks lots of things
