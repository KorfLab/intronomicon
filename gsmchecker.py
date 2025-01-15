import os
import sys

for path in os.listdir(sys.argv[1]):
	with open(f'{sys.argv[1]}/{path}') as fp: text = fp.read()
	lines = text.split('!')

	d = {}
	for line in lines[1:]:
		f = line.split()
		if len(f) == 0: continue
		f = line.split()
		k = f[0]
		v = ' '.join(f[2:])
		if k not in d: d[k] = v
		else: d[k] += f'; {v}'

	# sanity checks
	if d['Sample_channel_count'] != '1': sys.exit('unexpected channel count')
	if d['Sample_geo_accession'] != path[:-4]:
		sys.exit('filename does not match geo accession')

	# difficulties
	if ',' in d['Sample_taxid_ch1']: continue

	# sometimes present
	some = ('Sample_description', 'Sample_growth_protocol_ch1',
		'Sample_treatment_protocol_ch1')
	for s in some:
		if s not in d: d[s] = ''

	# always present
	output = []
	output.append(f"Source: {d['Sample_source_name_ch1']}")
	output.append(f"Title: {d['Sample_title']}")
	output.append(f"Description: {d['Sample_description']}")
	output.append(f"Molecule: {d['Sample_molecule_ch1']}")
	output.append(f"Selection: {d['Sample_library_selection']}")
	output.append(f"Strategy: {d['Sample_library_strategy']}")
	output.append(f"Characteristics: {d['Sample_characteristics_ch1']}")
	output.append(f"Growth: {d['Sample_growth_protocol_ch1']}")
	output.append(f"Treatment: {d['Sample_treatment_protocol_ch1']}")

	print('\t'.join(output))
