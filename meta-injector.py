
import argparse
import sqlite3

parser = argparse.ArgumentParser(description='metadata stuff')
parser.add_argument('database', help='database name (e.g. something.db)')
arg = parser.parse_args()

# connect to db
con = sqlite3.connect(arg.database)
cur = con.cursor()

# get a record that is not labeled and not locked
cur.execute('SELECT * from runs where labeled = 0 and locked = 0 limit 1')
x = cur.fetchall()
if len(x) == 0: sys.exit('no records to process')
srr, nts, seqs, exp_id, aligned, labeled, locked = x[0]

# create lock 
print('working on', srr)
cur.execute(f'UPDATE runs SET locked = 1 WHERE run_id = "{srr}"')
con.commit()

# get SRA metadata
# get GEO metadata
# do _something_ with them

con.close()