#!/usr/bin/env python3

import sys
import re
from urllib.request import urlopen

ARGV = sys.argv[1:];

for acc_arg in ARGV :

## check for acc format
  if (re.search(r'\|',acc_arg)):
    acc_info=acc_arg.split('|')  # sp|P09488|GSTM1_HUMAN
    acc = acc_info[1]            # P09488
  else:
    acc = acc_arg

  url = "https://www.uniprot.org/uniprot/" + acc + ".fasta"
#  print url
  fa_seq = urlopen(url).read().decode('utf-8')
  print(fa_seq)
