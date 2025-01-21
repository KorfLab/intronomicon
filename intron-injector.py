import argparse
import sqlite3
import sys

parser = argparse.ArgumentParser(description='aligns sequences, talks to db')
parser.add_argument('database', help='database name (e.g. something.db)')
parser.add_argument('--minlen', type=int, default=100,
	help='minimum read length [%(default)i]')
parser.add_argument('--platform', default='ILLUMINA',
	help='sequencing platform [%(default)s]')
arg = parser.parse_args()

# connect to db
con = sqlite3.connect(arg.database)
cur = con.cursor()

# get an SRR to align
cur.execute(f'SELECT srr_id, rlen, nts, platform FROM run INNER JOIN experiment ON experiment.srx_id = run.srx_id WHERE rlen >= {arg.minlen} and platform = "{arg.platform}" and aligned = 0 and locked = 0 limit 1')


x = cur.fetchall()
if len(x) == 0: sys.exit('no records to process')

print(x)

# create lock
#print('working on', srr, nts, rlen, exp_id)
#cur.execute(f'UPDATE runs SET locked = 1 WHERE run_id = "{srr}"')
#con.commit()


# retreive sequences
# do alignments
# QC the alignments
# inject introns

con.close()
