import os
import glob
import sys
import re
import pprint
import argparse

parser = argparse.ArgumentParser(description='GSE compare grouping')
parser.add_argument('meta_training_dir', help='Path to the meta-training directory')
parser.add_argument('-d','--data', help='Print out the info', action= 'store_true')
parser.add_argument('-i','--info', help='Print out the data', action= 'store_true')
parser.add_argument('-t','--table', help='Print out the table', action= 'store_true')
arg = parser.parse_args()

meta_training_dir = arg.meta_training_dir

if len(sys.argv) < 2:
    sys.exit(f"usage: {sys.argv[0]} <meta-training dir>")

data = {}
info = {}

def key_addition(key, name):

    key = tuple(key)
    # if there is no such key, create such key
    if key and key not in data[gse]:
        data[gse][key] = set()
        data[gse][key].add(name)
    # if there is such key, assign the key owner
    elif key:
        data[gse][key].add(name)

def read_html(filename, gse, name):
    with open(filename) as fp:
        key = []
        for line in fp:
            line = line.lstrip()
            
            #############################################
            ''' area for extracting those GSM from html'''
            #############################################
            
            # google chrome in pulling GSM; situation1
            # found in ian
            if line.startswith('<dt'):
                m_list = re.findall(r'<b>(GSM\d+)</b>', line)
                if m_list: key.extend(m_list)
                else: continue
                info[gse]['Total_gsms']        += len(m_list)
                info[gse][f'{name}_GSM_found'] += len(m_list)
                
            # google chrome in pulling GSM; situation2
            # found in ian
            if line.startswith('<dd'):
                m_list = re.findall(r'<b>(GSM\d+)</b>', line)
                if m_list: key.extend(m_list)
                else: continue
                info[gse]['Total_gsms']        += len(m_list)
                info[gse][f'{name}_GSM_found'] += len(m_list)
            
            # google chrome in pulling GSM; situation3
            # found in ian
            if line.startswith('<dl'):
                m_list = re.findall(r'<b>(GSM\d+)</b>', line)
                if m_list: key.extend(m_list)
                else: continue                
                info[gse]['Total_gsms']        += len(m_list)
                info[gse][f'{name}_GSM_found'] += len(m_list)
            
            # if there is no GSM yet, continue search GSM
            if not key: 
                line = line.strip()
            
            #############################################
            ''' area for extracting those tags from html'''
            #############################################
            
            if line.endswith('</p>\n'):
                tags = line[9:-5].split(', ')
            elif line.endswith('</p></dl>\n'):
                tags = line[9:-10].split(', ')
            elif line.endswith('</p></dl></body></html>'):
                tags = line[9:-23].split(', ')
            elif line.startswith('<p>Tags:'):
                tags = line[9:-23].split(', ')
            else:
                tags = None

            if tags is not None:
                info[gse][f'{name}_GSM_tagged'] += len(key)
                key.sort(key=lambda x: int(x[3:]))
                key.extend(tags)
                key_addition(key, name)
                key = []
            
        if key: 
            key_addition(key, name)
            key = []

meta_training_dir = sys.argv[1]
# exclude the src out from the dir
subdirs = [
    os.path.join(meta_training_dir, d)
    for d in os.listdir(meta_training_dir)
    # remove the source directory
    if os.path.isdir(os.path.join(meta_training_dir, d)) and d != "src"

]

if not subdirs:
    sys.exit(f"No valid subdirectories found in {meta_training_dir}")

for curator in subdirs:

    name = os.path.basename(curator.rstrip(os.sep))

    for file in glob.glob(os.path.join(curator, "*.html")):

        m = re.search(r'(GSE\d+)', file)

        gse = m.group(1)

        if gse not in info:
            info[gse] = {}
            info[gse]['Crews_count'] = 0
            info[gse]['Total_gsms']  = 0
            
        if gse not in data:
            data[gse] = {}

        if f'{name}_GSM_found' not in info[gse]:
            info[gse][f'{name}_GSM_found'] = 0

        if f'{name}_GSM_tagged' not in info[gse]:
            info[gse][f'{name}_GSM_tagged'] = 0
            
        info[gse]['Crews_count'] += 1

        if 'Crews_name' not in info[gse]:
            info[gse]['Crews_name'] = [name]
        else:
            info[gse]['Crews_name'].append(name)

        read_html(file, gse, name)

# if we wanna print out the raw data or the raw info

if arg.data:
    pprint.pprint(data)
if arg.info:
    pprint.pprint(info)

###################
''' error check '''
###################

# consensus on GSM number check
for gse in info:
    crews = info[gse]['Crews_name']
    set1  = set()
    set2  = set()
    # iterate to get different crew data
    for crew in crews:
        gsm        = info[gse][f'{crew}_GSM_found']
        gsm_tagged = info[gse][f'{crew}_GSM_tagged']
        try:
            if not set1:
                set1.add(gsm)
            else: 
                set1.remove(gsm)
                set1.add(gsm)
        # if anything wrong with the gsm part
        except:
            print(f"{gse}-----crews don't have form consensus on it")
        
        # check whether number of gsm found is match with number of gsm tagged
        # not working for source data
        if gsm + gsm_tagged != gsm * 2:
            print(f'{gse}-----data extraction in {crew} have a bug')
            
###################
'''formal output'''
###################

def scoring( number_of_gsms, number_of_crews, number_of_keys):
    assert( type(number_of_gsms)  == int )
    assert( type(number_of_crews) == int )
    
    gsms = number_of_gsms
    
    # if there's only one person did that so far
    if number_of_crews == 1:
        score = 1
        return score
    
    # ideal situation; for low gsms, that ain't sum of 2 or 3, it shall be between 1 and 2
    if gsms < 6:    ideal_key_number = gsms/1.5
    # ideal situation; for high gsms, that prob are sum of 2 or 3, it shall be between 2 and 3
    elif 15 > gsms >= 6: ideal_key_number = gsms/2.5
    # there is one with 20gsms that 5 in group
    elif gsms >= 15: ideal_key_number = gsms/4.5
    
    # calculate score
    adjustment = 0 # parameter for future thoughts 
    score = abs ( (ideal_key_number - number_of_keys)/(ideal_key_number) ) + adjustment + 1 #plus-one-for-better-output-format
    return score


if arg.table:
    rows = []

    for gse in data:
        keys  =    len( data[gse] )
        crews =    int ( info[gse]['Crews_count'] )
        names =    info[gse]['Crews_name']
        names =    ", ".join(names)
        all_gsms = int ( info[gse]['Total_gsms']  )
        gsms  = int(all_gsms/crews)
        score = scoring( all_gsms, crews, keys)
        rows.append((gse, gsms, score, crews, keys, names))
    
    rows = sorted(rows, key=lambda row: (row[2], row[3], row[4]), reverse=True)

    # header
    header = f"{'GSE':<10} {'GSMs':<10} {'Score':<10} {'Crews':<7} {'Keys':<7} {'Crew Names':<40}"
    print(header)
    print("-" * len(header))

    for gse, gsms, score, crews, keys, names in rows:
        print(f"{gse:<10} {gsms:<10} {score:<10.4g} {crews:<7} {keys:<7} {names:<40}")