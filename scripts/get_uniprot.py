#!/usr/bin/python

import sys
from urllib import urlopen

ARGV = sys.argv[1:];

for acc in ARGV :
  url = "https://www.uniprot.org/uniprot/" + acc + ".fasta"
#  print url
  fa_seq = urlopen(url).read()
  print fa_seq
