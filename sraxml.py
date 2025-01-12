import re
import sys
import xml.etree.ElementTree as ET

def descend_tree(node, prev):
	if len(node) == 0: return prev
	objects = []
	for item in node:
		obj = {'tag': item.tag}
		if item.text and re.match(r'\S',  item.text): obj['txt'] = item.text
		if item.attrib: obj['att'] = item.attrib
		contents = descend_tree(item, [])
		if len(contents) > 0: obj['has'] = contents
		objects.append(obj)
	return objects

def read_xml(fp):
	tree = ET.parse(fp)
	root = tree.getroot()
	data = {'tag': root.tag}
	if re.search(r'\S', root.text): data['txt'] = root.text
	if root.attrib: data['att'] = root.attrib
	contents = descend_tree(root, [])
	if contents: data['has'] = contents
	return data

def get(src, tag, att=None):
	if src['tag'] != tag: return None
	if att:
		return src['att'][att]
	else:
		if 'txt' not in src: return None
		else:                return src['txt']
	sys.exit('impossible')

def read(fp):
	obj = {
		'srx_id': None,
		'srp_id': None,
		'sam_id': None,
		'geo_id': None,
		'paired': None,
		'platform': None,
		'model': None,
		'taxid': None,
		'info': None,
		'runs': [],
	}

	eps = read_xml(fp)
	pkg = eps['has'][0]
	if len(pkg['has']) != 7: return None, 'incomplete experiment package'

	######################
	# EXPERIMENT SECTION #
	######################
	exp = pkg['has'][0]
	obj['srx_id'] = get(exp['has'][0]['has'][0], 'PRIMARY_ID')
	obj['srp_id'] = get(exp['has'][2], 'STUDY_REF', att='accession')

	# get library data
	des = exp['has'][3]
	obj['sam_id'] = get(des['has'][1], 'SAMPLE_DESCRIPTOR', att='accession')
	lib = {}
	for thing in des['has'][2]['has']:
		if 'tag' in thing and 'txt' in thing:
			tag = thing['tag']
			val = thing['txt']
		elif 'tag' in thing and 'has' in thing:
			if len(thing['has']) == 1 and 'tag' in thing['has'][0]:
				val = thing['has'][0]['tag']
			else:
				sys.exit('wtf 1')
		elif 'tag' in thing:
			tag = thing['tag']
			val = None # sometimes the tag has no value...
		else:
			print(json.dumps(thing, indent=4))
			sys.exit('wtf 2')
		lib[tag] = val

	# library-level checks
	if lib['LIBRARY_STRATEGY'] != 'RNA-Seq':
		return None, lib['LIBRARY_STRATEGY']
	if lib['LIBRARY_SELECTION'] is None:
		return None, 'unknown strategy'
	if lib['LIBRARY_SOURCE'] != 'TRANSCRIPTOMIC':
		return None, lib['LIBRARY_SOURCE']
	obj['paired'] = 0 if lib['LIBRARY_SELECTION'] == 'SINGLE' else 1

	# platform
	obj['platform'] = exp['has'][4]['has'][0]['tag']
	obj['model'] = exp['has'][4]['has'][0]['has'][0]['txt']

	# geo (may be buried)
	for sid in des['has'][1]['has'][0]['has']:
		if 'att' in sid and 'namespace' in sid['att'] \
			and sid['att']['namespace'] == 'GEO': obj['geo_id'] = sid['txt']

	##################
	# SAMPLE SECTION #
	##################
	samp = pkg['has'][4]

	# saying sample is organized like garbage is a disservice to garbage
	taxid = None
	satt = {} # sample attributes (free text tags and values...)
	for thing in samp['has']:
		if thing['tag'] == 'SAMPLE_NAME':
			taxid = get(thing['has'][0], 'TAXON_ID')
		if thing['tag'] == 'SAMPLE_ATTRIBUTES':
			for x in thing['has']:
				if 'txt' not in x['has'][1]:
					print(f'warning, sample section', file=sys.stderr)
					continue
				tag = x['has'][0]['txt']
				val = x['has'][1]['txt']
				satt[tag] = val
	if taxid is None: sys.exit('wtf taxid')
	obj['taxid'] = taxid
	obj['info'] = satt

	##################
	# RUNSET SECTION #
	##################
	runs = pkg['has'][6]
	for r in runs['has']:
		srr_id = get(r, 'RUN', att='accession')
		spots = get(r, 'RUN', att='total_spots')
		bases = get(r, 'RUN', att='total_bases')
		size = get(r, 'RUN', att='size')
		date = get(r, 'RUN', att='published').split(' ')[0]
		for val in (srr_id, spots, bases, size):
			if val is None: sys.exit('wtf run')
			obj['runs'].append({
				'run_id': srr_id,
				'nts': int(bases),
				'seqs': int(spots),
				'date': date})

	###################
	# UNUSED SECTIONS #
	###################
	sub = pkg['has'][1]
	org = pkg['has'][2]
	stdy = pkg['has'][3]
	pool = pkg['has'][5]

	return obj, 'ok'
