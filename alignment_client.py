import argparse
import sqlite3

parser = argparse.ArgumentParser(description='aligns sequences, talks to db')
parser.add_argument('--database', default='intronomicon.db',
	help='database name [%(default)s]')
arg = parser.parse_args()

# connect to db
# loop
	# ask for sra - query status table SRR -> [done, free, working]
		# it hasn't been done yet
		# it isn't being worked on
	# check out sra - prevent others from working on this
	# do the alignments
	# push intron data to db


con = sqlite3.connect(arg.database)
cur = con.cursor()

cur.execute('SELECT * from runs ORDER BY bases')
for x in cursor.fetchall():
	pass
