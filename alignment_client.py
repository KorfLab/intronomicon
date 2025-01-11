import argparse
import sqlite3

parser = argparse.ArgumentParser(description='aligns sequences, talks to db')
parser.add_argument('database', help='database name (e.g. something.db)')
arg = parser.parse_args()

# connect to db
# loop
	# ask for sra - query status table SRR -> [done, free, working]
		# it hasn't been done yet
		# it isn't being worked on
	# check out sra - prevent others from working on this
	# do the alignments
	# do whatever post-processing is required (e.g. exon overlap length)
	# push intron data to db


con = sqlite3.connect(arg.database)
cur = con.cursor()

cur.execute('SELECT * from runs ORDER BY bases')
for x in cursor.fetchall():
	pass
