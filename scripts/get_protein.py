#!/usr/bin/env python3

## get_protein_www.py -- 
## get a protein sequence from the Uniprot or NCBI/Refseq web sites using the accession
##

## modified to work with urllib.request 7-Nov-2022
## modified to allow argparse arguments for identifier 20-Mar-2023

import argparse
import sys
import re
import textwrap
import time
import urllib.request
import urllib.error

def main():

  ncbi_url = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?" 
  uniprot_url = "https://rest.uniprot.org/uniprotkb/"
  sub_range = ''

  parser=argparse.ArgumentParser(description='get protein sequences from uniprot/ncbi')
  parser.add_argument('--id', help='substitute id',action='store',default='')
  parser.add_argument('accs', nargs='*', help='accessions')

  args=parser.parse_args()

  for acc in args.accs:

    if (re.search(r':',acc)):
      (acc, sub_range) = acc.split(':')

    if (re.match(r'^(sp|tr|iso|ref)\|',acc)):
        acc=acc.split('|')[1]

    if (re.match(r'[A-Z]P_\d+',acc)):   # get refseq
      db_type="protein"

      seq_args = "db=%s&id=" % (db_type) + acc  + "&rettype=fasta"

      url_string = ncbi_url + seq_args

    else:				# get uniprot
      acc_fields = acc.split('|')
      if (len(acc_fields)==1):
        url_string = uniprot_url + acc + ".fasta"
      else:
        url_string = uniprot_url + acc_fields[0] + ".fasta"

    try: 
      req = urllib.request.urlopen(url_string)
    except urllib.error.URLError as e:
      seq_html = ''
      sys.stderr.write(e.read().decode('utf-8')+'\n')
      continue

    else:
      seq_html=req.read().decode('utf-8')

    time.sleep(0.3)

    if (not sub_range):

      if (args.id):
        seq_html = re.sub('>','>%s '%(args.id),seq_html)

      print(seq_html)
    else:
      (start, stop) = sub_range.split('-')

      (start, stop) = (int(start), int(stop))

      lines = seq_html.split('\n')

      header=lines[0]
      seq = ''.join(lines[1:])

      if (start > 0):
        start -= 1

      new_seq = seq[start:stop]
      ## print the header
      if (start > 0):
        print("%s @C:%d" %(header, start+1))
      else:
        print(header)

      print('\n'.join(textwrap.wrap(new_seq)))

if __name__ == '__main__':
    main()
