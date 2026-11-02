#!/usr/bin/env python3

## get_protein_www.py -- 
## get a protein sequence from the Uniprot or NCBI/Refseq web sites using the accession
##

## modified to work with urllib.request 7-Nov-2022
## modified to allow argparse arguments for identifier 20-Mar-2023
## modified 11-Aug-2026 to de-taint accessions

import argparse
import sys
import os
import re
import textwrap
import time
import urllib.request
import urllib.error

def main():

  ncbi_url = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?" 
  uniprot_url = "https://rest.uniprot.org/uniprotkb/"
  uniparc_url = "https://rest.uniprot.org/uniparc/"
  sub_range = ''

  parser=argparse.ArgumentParser(description='get protein sequences from uniprot/ncbi')
  parser.add_argument('--id', help='substitute id',action='store',default='')
  parser.add_argument('--key', help='API_key',action='store',default='')
  parser.add_argument('--log', help='log download',action='store_true',default=False)
  parser.add_argument('accs', nargs='*', help='accessions')

  args=parser.parse_args()

  r_size = 100

  if (args.key):
    ncbi_url += f'api_key={args.key}&'
  elif os.getenv('NCBI_API_KEY'):
    ncbi_url += f'api_key={os.getenv("NCBI_API_KEY")}&'

  if (args.accs[0].startswith("@")):
    accs = []
    with open(args.accs[0][1:],'r') as f_in:
      for line in f_in:
        accs.append(line.strip('\n'))

    acc_list = []

    for ix in range(0,len(accs),r_size):
      mx = ix + r_size
      if (mx > len(accs)):
        mx = len(accs)
##      print(f"extracting: {ix}:{mx}",file=sys.stderr)
      acc_list.append(','.join(accs[ix:mx]))

  else:
    acc_list = args.accs
    
  for acc in acc_list:

    ## detaint acc -- do not allow ':' or '|'
    acc = re.sub(r'[^A-Za-z_0-9\.\,\|]','',acc)

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

      if (re.search(r'^UPI',acc)):
        uniprot_url = uniparc_url

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
        ## detaint since output is being printed:
        args.id = re.sub(r'[^A-Za-z_0-9\|]','',args.id)
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
