#!/usr/bin/env python3

## get_protein.py -- 
## get a protein sequence from Uniprot or NCBI/Refseq using the accession
##

import sys
import re
import textwrap
import MySQLdb
import time
from urllib.request import Request, urlopen
from urllib.error import URLError

ncbi_url = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?" 
uniprot_url = "https://www.uniprot.org/uniprot/"

(host, user, passwd, db) = ("wrpxdb.its.virginia.edu","web_user","fasta_www","seqdb_demo2")

db = MySQLdb.connect(host,user,passwd,db)
cur = db.cursor()

get_refseq_sql = '''select acc, descr, seq from annot join protein using(prot_id) where acc="%s"'''
get_up_sql = '''select db, acc, id, descr, seq from uniprot.annot2 join uniprot.protein using(acc) where acc="%s"'''
get_up_iso_sql = '''select db, acc, id, descr, seq from uniprot.annot2_iso join uniprot.protein_iso using(acc) where acc="%s"'''

def get_ncbi_sql(acc, cur):
    result = cur.execute(get_refseq_sql%(acc,))
    if (result):
        acc, descr, seq = cur.fetchone()
        seq = re.sub(r'(.{60})',r'\g<1>\n',seq)
        return ">%s %s\n%s\n"%(acc, descr, seq)
    else:
        return False

def get_ncbi_www(acc):
    db_type="protein"
    seq_args = "db=%s&id=" % (db_type) + acc  + "&rettype=fasta"

##   seq_html = urlopen(ncbi_url + seq_args)

    req = Request(ncbi_url+seq_args)

    try: 
        resp = urlopen(req)
    except URLError as e:
        seq_html = ''
        if hasattr(e, 'reason'):
            sys.stderr.write("Sequence [%s] not found: %s\n"%(acc,str(e.reason)))
        elif hasattr(e, 'code'):
            sys.stderr.write("Sequence [%s] not found: %s\n"%(acc,str(e.ecode)))
        else:
            sys.stderr.write("Sequence [%s] not found\n"%(acc))
        exit(0)

    else:
        seq_html=resp.read().decode('utf-8')

    time.sleep(0.3)
    return seq_html

def get_uniprot_sql(acc, cur):
    result = cur.execute(get_up_sql%(acc,))
    if (not result):
        result = cur.execute(get_up_iso_sql%(acc,))
    if (result):
        db, acc, id, descr, seq = cur.fetchone()
# not needed because of final textwrap
#        seq = re.sub(r'(.{60})',r'\g<1>\n',seq)
        return ">%s|%s|%s %s\n%s\n"%(db, acc, id, descr, seq)
    else:
        return False

def get_uniprot_www(acc):
    seq_html = urlopen(uniprot_url + acc + ".fasta").read().decode('utf-8')
    return seq_html

sub_range = ''
for acc in sys.argv[1:]:

    if (re.search(r':',acc)):
        (acc, sub_range) = acc.split(':')

    if (re.match(r'^(sp|tr|iso|ref)\|',acc)):
        acc=acc.split('|')[1]

    if (re.match(r'[NXYW]P_',acc)):
        seq_html = get_ncbi_sql(acc,cur)
        if (not seq_html):
            seq_html = get_ncbi_www(acc)
    else:
        seq_html = get_uniprot_sql(acc,cur)
        if (not seq_html):
            seq_html = get_uniprot_www(acc)

    if not seq_html:
        continue;

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
                    print("%s @C%d" %(header, start+1))
                else:
                    print(header)

                print('\n'.join(textwrap.wrap(new_seq)))

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
    print("%s @C:%d" %(header, start+1))
else:
    print(header)

print('\n'.join(textwrap.wrap(new_seq)))
