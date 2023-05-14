#!/usr/bin/env python

## get_protein_sql.py -- 
## get a protein sequence from a local Uniprot or NCBI/Refseq mySQL database using the accession
##

## modified to work with mysql.connector 7-Nov-2022

import sys
import re
import textwrap
import mysql.connector

db_r = mysql.connector.connect(db='seqdb_demox', host='wrpxdb.bioch.virginia.edu', user='web_user', passwd='fasta_www')
db_u = mysql.connector.connect(db='uniprot', host='wrpxdb.bioch.virginia.edu', user='web_user', passwd='fasta_www')

cur1_r = db_r.cursor(dictionary=True, buffered=True)
cur1_u = db_u.cursor(dictionary=True, buffered=True)

sql_get_uniprot='select db, acc, id, descr, seq from annot2 join protein using(acc) where acc="%s"'
sql_get_refseq ='select db, acc, "" as id, descr, seq from seqdb_demox.annot join seqdb_demox.protein using(prot_id) where acc="%s"'

sub_range = ''
for acc in sys.argv[1:]:

  if (re.search(r':',acc)):
    (acc, sub_range) = acc.split(':')

  if (re.match(r'^(sp|tr|iso|ref)\|',acc)):
      acc=acc.split('|')[1]

  if (re.match(r'[A-Z]P_\d+',acc)):
    sql_get_prot=sql_get_refseq
    cur1 = cur1_r
  else:
    sql_get_prot=sql_get_uniprot
    cur1 = cur1_u

  cur1.execute(sql_get_prot%(acc,))

  row = cur1.fetchone()

  if (not row):
    sys.stderr.write("*** %s *** not found\n"%(acc))
    exit(1)

  header = ">%s|%s"%(row['db'],row['acc'])
  if (row['id']):
    header += "|%s"%(row['id'])

  header += " "+row['descr']

  start = 0
  if (sub_range):
    start, stop = sub_range.split('-')
    start, stop = int(start), int(stop)

  if (start > 0):
    seq = row['seq'][start-1:stop]
    print("%s @C%d" %(header, start+1))
  else:
    seq = row['seq']

  print(header)
  print('\n'.join(textwrap.wrap(seq)))
