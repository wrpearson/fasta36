#!/usr/bin/env python3

## get_protein.py -- 
## get a protein sequence from Uniprot or NCBI/Refseq using the accession
##

import sys
import re
import textwrap
import MySQLdb
import time
import requests

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

    try: 
        req = requests.get(ncbi_url+seq_args)
    except requests.exceptions.RequestException as e:
        seq_html = ''
        sys.stderr.print(e.response.text+'\n')
    else:
        seq_html=req.text

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
    up_req = requests.get(uniprot_url + acc + ".fasta")
    seq_html = up_req.text

    if (re.search(r'Error',seq_html)):
        sys.stderr.write("*** ERROR : %s not found\n"%(acc))
        seq_html=''
    return seq_html

def print_seq(seq_html, sub_range):

    if not seq_html:
        return

    lines = seq_html.split('\n');

    header = lines[0]
    seq = ''.join(lines[1:])

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

    print_seq(seq_html, sub_range)


