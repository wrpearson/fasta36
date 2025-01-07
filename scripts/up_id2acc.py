#!/usr/bin/env python3

## get_protein_www.py -- 
## get a protein sequence from the Uniprot or NCBI/Refseq web sites using the accession
##

## modified to work with urllib.request 7-Nov-2022
## modified to allow argparse arguments for identifier 20-Mar-2023

import argparse
import time
import sys
import re
import json
from urllib import request, parse, error

def main():

  uniprot_url = "https://rest.uniprot.org/idmapping/"
  uniprot_run_url = uniprot_url+"run"
  uniprot_status_url = uniprot_url+"status/"
  uniprot_results_url = uniprot_url+"results/"

  parser=argparse.ArgumentParser(description='get protein sequences from uniprot/ncbi')
  parser.add_argument('--fail', action='store_true', help='show failures', default=False)
  parser.add_argument('ids', nargs='*', help='Uniprot IDs')

  args=parser.parse_args()

  id_str = ','.join(args.ids)

  url_data = {'ids':id_str, 'from':'UniProtKB_AC-ID', 'to':'UniProtKB'}
  enc_data = parse.urlencode(url_data).encode()

  try: 
    req = request.urlopen(uniprot_run_url, enc_data)
  except error.URLError as e:
    sys.stderr.write(e.read().decode('utf-8')+'\n')
    exit(1)

  else:
    jobId_html=req.read().decode('utf-8')

  jobId_json = json.loads(jobId_html)

  if ('jobId' not in jobId_json):
    sys.stderr.write("Bad jobId:\n%s\n"%(str(jobId_json)))
    exit(1)
  else:
    jobId=jobId_json['jobId']
    
  req_cnt = 0
  while (req_cnt < 10):

    req_cnt += 1
    try: 
      req = request.urlopen(uniprot_status_url+jobId)
    except error.URLError as e:
      sys.stderr.write(e.read().decode('utf-8')+'\n')
      exit(1)

    else:
      status_html=req.read().decode('utf-8')

    status_json = json.loads(status_html)

    if ('results' in status_json):
      break

    if ('jobStatus' in status_json and status_json['jobStatus'] == "FINISHED"):
      break

    time.sleep(1)

  ## have results or finished, get result
  
  if (args.fail):
    if ('failedIds' in status_json):
      for fail_id in status_json['failedIds']:
        print('\t'.join([fail_id,'notFound']))

  if ('results' in status_json):
    for result in status_json['results']:
      print('\t'.join([result['from'], result['to']['primaryAccession']]))
  else:
    try: 
      req = request.urlopen(uniprot_results_url+jobId)
    except error.URLError as e:
      sys.stderr.write(e.read().decode('utf-8')+'\n')
      exit(1)

    else:
      results_html=req.read().decode('utf-8')

    results_json = json.loads(results_html)
      
    for result in status_json['results']:
      print('\t'.join([result['from'], result['to']]))

if __name__ == '__main__':
    main()
