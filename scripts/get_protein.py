#!/usr/bin/python

import sys
import re
import textwrap
from urllib2 import urlopen

ncbi_url = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?" 
uniprot_url = "https://www.uniprot.org/uniprot/"

sub_range = ''
for acc in sys.argv[1:]:

  if (re.search(r':',acc)):
    (acc, sub_range) = acc.split(':')

  if (re.match(r'^(sp|tr|iso|ref)\|',acc)):
      acc=acc.split('|')[1]

  if (re.match(r'[NX]P_',acc)):
    db_type="protein"

    seq_args = "db=%s&id=" % (db_type) + ",".join(sys.argv[1:])  + "&rettype=fasta"
    seq_html = urlopen(ncbi_url + seq_args).read()
  else:
    seq_html = urlopen(uniprot_url + acc + ".fasta").read()

  header=''
  seq = ''
  for line in seq_html.split('\n'):
    if (line and line[0]=='>'):
      # print out old one if there
      if (header):
        if (sub_range):
          start, stop = sub_range.split('-')
          start, stop = int(start), int(stop)
          if (start > 0):
            start -= 1
          new_seq = seq[start:stop]
        else:
          start = 0
          new_seq = seq

        if (start > 0):
          print "%s @C%d" %(header, start+1)
        else:
          print header
        print '\n'.join(textwrap.wrap(new_seq))

      header = line;
      seq = ''
    else:
      seq += line

start=0
if (sub_range):
  start, stop = sub_range.split('-')
  start, stop = int(start), int(stop)
  if (start > 0):
    start -= 1
    new_seq = seq[start:stop]
else:
  new_seq = seq

if (start > 0):
  print "%s @C:%d" %(header, start+1)
else:
  print header

print '\n'.join(textwrap.wrap(new_seq))
