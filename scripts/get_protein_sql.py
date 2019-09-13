#!/usr/bin/env python

## get_protein_sql.py -- 
## get a protein sequence from a local Uniprot or NCBI/Refseq mySQL database using the accession
##

import sys
import re
import textwrap
import MySQLdb.cursors

db = MySQLdb.connect(db='uniprot', host='wrpxdb.bioch.virginia.edu', user='web_user', passwd='fasta_www',
                     cursorclass=MySQLdb.cursors.DictCursor)

cur1 = db.cursor()
cur2 = db.cursor()

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
  else:
    sql_get_prot=sql_get_uniprot


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
