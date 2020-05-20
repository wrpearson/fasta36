#!/usr/bin/env python3

import sys
import re
import requests

db_type="protein"
if (re.match(r'[A-Z]M_\d+',sys.argv[1])):
    db_type="nucleotide"

acc_str = ",".join(sys.argv[1:])


seq_url = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?"
seq_args = "db=%s&id=" % (db_type) + acc_str  + "&rettype=fasta"
url_string = seq_url+seq_args

try: 
    req = requests.get(url_string)
except requests.exceptions.RequestException as e:
    seq_html = ''
    sys.stderr.print(e.response.text+'\n')
else:
    seq_html=req.text
    print(seq_html,end='')
