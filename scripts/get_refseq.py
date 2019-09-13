#!/usr/bin/env python3

import sys
import re
from urllib.request import Request, urlopen
from urllib.error import URLError

db_type="protein"
if (re.match(r'[A-Z]M_\d+',sys.argv[1])):
    db_type="nucleotide"

acc_str = ",".join(sys.argv[1:])


seq_url = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?"
seq_args = "db=%s&id=" % (db_type) + acc_str  + "&rettype=fasta"
url_string = seq_url+seq_args

req = Request(url_string)

try: 
    resp = urlopen(req)
except URLError as e:
    seq_html = ''
    if hasattr(e, 'reason'):
      sys.stderr.write("Sequence [%s] not found: %s\n"%(acc_str,str(e.reason)))
    elif hasattr(e, 'code'):
      sys.stderr.write("Sequence [%s] not found: %s\n"%(acc_str,str(e.ecode)))
    else:
      sys.stderr.write("Sequence [%s] not found\n"%(acc_str))
    exit(0)

else:
    seq_html=resp.read().decode('utf-8')

print(seq_html,end='')
