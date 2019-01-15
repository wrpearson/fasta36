#!/usr/bin/python

import sys
import re
from urllib2 import urlopen


db_type="protein"
if (re.match(r'[NX]M_',sys.argv[1])):
    db_type="nucleotide"

seq_url = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?"
seq_args = "db=%s&id=" % (db_type) + ",".join(sys.argv[1:])  + "&rettype=fasta"

seq_html = urlopen(seq_url + seq_args).read()

print seq_html
