import korflab

def get(src, tag, att=None):
	if src['tag'] != tag: return None
	if att:
		return src['att'][att]
	else:
		if 'txt' not in src: return None
		else:                return src['txt']
	sys.exit('impossible')

def read_sra_xml(fp):
	obj = {
		'srx_id': None,
		'srp_id': None,
		'sam_id': None,
		'paired': None,
		'platform': None,
		'model': None,
		'taxid': None,
		'info': None,
		'runs': [],
	}

	eps = korflab.read_xml(fp)
	pkg = eps['has'][0]
	if len(pkg['has']) != 7: return None, 'incomplete experiment package'

	######################
	# EXPERIMENT SECTION #
	######################
	exp = pkg['has'][0]
	obj['srx_id'] = get(exp['has'][0]['has'][0], 'PRIMARY_ID')
	obj['srp_id'] = get(exp['has'][2], 'STUDY_REF', att='accession')
	des = exp['has'][3]
	obj['sam_id'] = get(des['has'][1], 'SAMPLE_DESCRIPTOR', att='accession')

	# library info is not always in the same order
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
	obj['paired'] = False if lib['LIBRARY_SELECTION'] == 'SINGLE' else True

	# platform
	obj['platform'] = exp['has'][4]['has'][0]['tag']
	obj['model'] = exp['has'][4]['has'][0]['has'][0]['txt']

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
				tag = x['has'][0]['txt']
				val = x['has'][1]['txt']
				satt[tag] = val
	if taxid is None: sys.exit('wtf taxid')
	obj['taxid'] = taxid

	##################
	# RUNSET SECTION #
	##################
	runs = pkg['has'][6]
	for r in runs['has']:
		srr_id = get(r, 'RUN', att='accession')
		spots = get(r, 'RUN', att='total_spots')
		bases = get(r, 'RUN', att='total_bases')
		size = get(r, 'RUN', att='size')
		for val in (srr_id, spots, bases, size):
			if val is None: sys.exit('wtf run')
			obj['runs'].append({
				'srr_id': srr_id,
				'seqs': int(spots),
				'length': int(bases)/int(spots)})

	###################
	# UNUSED SECTIONS #
	###################
	sub = pkg['has'][1]
	org = pkg['has'][2]
	stdy = pkg['has'][3]
	pool = pkg['has'][5]

	return obj, 'ok'
