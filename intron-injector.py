import argparse
import sqlite3
import sys

parser = argparse.ArgumentParser(description='aligns sequences, talks to db')
parser.add_argument('database', help='database name (e.g. something.db)')
arg = parser.parse_args()

# connect to db
con = sqlite3.connect(arg.database)
cur = con.cursor()

# get a record that is not aligned and not locked
cur.execute('SELECT * from run where aligned = 0 and locked = 0 limit 1')
x = cur.fetchall()
if len(x) == 0: sys.exit('no records to process')
srr, nts, rlen, exp_id, aligned, labeled, locked = x[0]

# create lock
print('working on', srr, nts, rlen, exp_id)
#cur.execute(f'UPDATE runs SET locked = 1 WHERE run_id = "{srr}"')
#con.commit()

# do alignments

con.close()
