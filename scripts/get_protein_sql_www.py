#!/usr/bin/env python3

## get_protein.py -- 
## get a protein sequence from Uniprot or NCBI/Refseq using the accession
##

## modified to work with mysql.connector, urllib.request 7-Nov-2022
##

import sys
import re
import textwrap
import mysql.connector
import time
import urllib.request
import urllib.error

ncbi_url = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?" 
uniprot_url = "https://rest.uniprot.org/uniprotkb/"

(host, user, passwd) = ("wrpxdb.bioch.virginia.edu","web_user","fasta_www")

db_r = mysql.connector.connect(host=host,user=user,password=passwd,database="seqdb_demox")
db_u = mysql.connector.connect(host=host,user=user,password=passwd,database="uniprot")

cur_r = db_r.cursor(dictionary=True, buffered=True)
cur_u = db_u.cursor(dictionary=True, buffered=True)

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
        req = urllib.request.urlopen(ncbi_url+seq_args)
    except urllib.error.URLError as e:
        seq_html = ''
        sys.stderr.print(e.read().decode('utf-8')+'\n')
    else:
        seq_html=req.read().decode('utf-8')

    time.sleep(0.3)
    return seq_html

def get_uniprot_sql(acc, cur):

    cur.execute(get_up_sql%(acc,))

    row = cur.fetchone()

    if (row):
        return ">%s|%s|%s %s\n%s\n"%(row['db'], row['acc'], row['id'], row['descr'], row['seq'])
    else:
        cur.execute(get_up_iso_sql%(acc,))

        db, acc, id, descr, seq = cur.fetchone()

        if (row):
            return ">%s|%s|%s %s\n%s\n"%(row['db'], row['acc'], row['id'], row['descr'], row['seq'])
        else:
            return False

def get_uniprot_www(acc):

    try:
        up_req = urllib.request.urlopen(uniprot_url + acc + ".fasta")
    except urllib.error.URLError as e:
        seq_html = ''
        sys.stderr.print(e.read().decode('utf-8')+'\n')
    else:
        seq_html=up_req.read().decode('utf-8')

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
        seq_html = get_ncbi_sql(acc,cur_r)
        if (not seq_html):
            seq_html = get_ncbi_www(acc)
    else:
        seq_html = get_uniprot_sql(acc,cur_u)
        if (not seq_html):
            seq_html = get_uniprot_www(acc)

    print_seq(seq_html, sub_range)


