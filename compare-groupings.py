import glob
import json
import re
import sys


if len(sys.argv) < 3: sys.exit(f'usage: {sys.argv[0]} <dir> <dir> ...')

def read_html(filename):
	with open(filename) as fp:
		groups = []
		group = []
		tags = None
		for line in fp:
			line = line.lstrip()
			if line.startswith('<dl'):
				if group: groups.append({'group': group, 'tags': tags})
				group = []
			if line.startswith('<dd'):
				m = re.search(r'(GSM\d+)', line)
				group.append(m.group(1))
			if line.startswith('<p>Tags:'):
				tags = line[9:-5].split(', ')
		if group: groups.append({'group': group, 'tags': tags})
	return groups


data = {}
for curator in sys.argv[1:]:
	if curator not in data: data[curator] = {}
	for file in glob.glob(f'{curator}/*.html'):
		m = re.search(r'(GSE\d+)', file)
		gse = m.group(1)
		data[curator][gse] = read_html(file)

print(json.dumps(data, indent=4))
