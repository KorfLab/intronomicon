
import argparse
import sqlite3

parser = argparse.ArgumentParser(description='metadata table maker')
parser.add_argument('database', help='database name (e.g. something.db)')
arg = parser.parse_args()

# connect to db
con = sqlite3.connect(arg.database)
cur = con.cursor()

# get all experiment
exp = {}
cur.execute('SELECT * from experiment')
for srx, gsm, stx, gtx, plt, mod, par in cur.fetchall():
	exp[srx] = (gsm, stx, gtx, plt, mod)

# get some runs from different experiments
used = set()
for srr, bases, nts, srx in cur.execute('SELECT * from runs limit 5'):
	if srx in used: continue


	# some minimum seq length?
	# some limit to platforms?
	used.add(srx)


