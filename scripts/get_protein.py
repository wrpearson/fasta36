#!/usr/bin/python

import sys
import re
from urllib2 import urlopen

ncbi_url = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?" 
uniprot_url = "https://www.uniprot.org/uniprot/"

for acc in sys.argv[1:]:
  if (re.match(r'[NX]P_',acc)):
    db_type="protein"

    seq_args = "db=%s&id=" % (db_type) + ",".join(sys.argv[1:])  + "&rettype=fasta"
    seq_html = urlopen(ncbi_url + seq_args).read()
  else:
    seq_html = urlopen(uniprot_url + acc + ".fasta").read()

  print seq_html

    
