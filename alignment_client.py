import argparse
import sqlite3

parser = argparse.ArgumentParser(description='aligns sequences, talks to db')
parser.add_argument('--database', default='intronomicon.db',
	help='database name [%(default)s]')
#parser.add_argument('--min
arg = parser.parse_args()

# connect to db
# loop
	# ask for sra - query status table SRR -> [done, free, working]
		# it hasn't been done yet
		# it isn't being worked on
	# check out sra - prevent others from working on this
	# do the alignments
	# push intron data to db


conn = sqlite3.connect(arg.database)
cursor = conn.cursor()

cursor.execute('SELECT * from runs ORDER BY bases')
for rid, eid, nts, seqs, filesize in cursor.fetchall():
	if seqs != 0:
		print(nts//seqs)

"""
CREATE TABLE runs(
        run_ID TEXT PRIMARY KEY,
        exp_ID TEXT,
        bases INTEGER,
        spots INTEGER,
        file_size INTEGER,
        FOREIGN KEY (exp_ID) REFERENCES experiment (exp_ID)
"""
