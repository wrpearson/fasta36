#!/usr/bin/env python3

import sys
import re
import textwrap
import argparse
import MySQLdb.cursors

db = MySQLdb.connect(db='uniprot', host='wrpxdb.bioch.virginia.edu', user='web_user', passwd='fasta_www',
                     cursorclass=MySQLdb.cursors.DictCursor)

cur1 = db.cursor()
cur2 = db.cursor()
get_iso_acc='select acc from annot2_iso where prim_acc="%s"'
get_fasta_info='select db, acc, id, descr, seq from annot2 join protein using(acc) where acc="%s"'
get_iso_fasta_info='select db, acc, id, descr, seq from annot2_iso join protein_iso using(acc) where prim_acc="%s"'

fasta_seqs=[]

for acc in sys.argv[1:]:

  if (re.search(r':',acc)):
    (acc, sub_range) = acc.split(':')

  if (re.match(r'^(sp|tr|iso|ref)\|',acc)):
      acc=acc.split('|')[1]

  cur1.execute(get_fasta_info%(acc,))
  row = cur1.fetchone()
  if (row):
    fasta_seqs.append(row)
  else:
    sys.stderr.write("***error*** %s sequence not found\n"%(acc))
    continue

  cur2.execute(get_iso_fasta_info%(acc,))
  for row in cur2:
    fasta_seqs.append(row)

  for row in fasta_seqs:
    print(">%s|%s|%s %s"%(row['db'],row['acc'],row['id'],row['descr']))
    print('\n'.join(textwrap.wrap(row['seq'])))




