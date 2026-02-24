#!/usr/bin/env python3

################################################################
# copyright (c) 2022 by William R. Pearson and The Rector &
# Visitors of the University of Virginia */
################################################################
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
# http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing,
# software distributed under this License is distributed on an "AS
# IS" BASIS, WITHOUT WRRANTIES OR CONDITIONS OF ANY KIND, either
# express or implied.  See the License for the specific language
# governing permissions and limitations under the License. 
################################################################

## updated 14-April-2024 to work on current InterPro API

# ann_pfam_www0.py takes an annotation file from fasta36 -V with a line of the form:

# sp|P0810|GSTM1_RAT [tab] seqlen
#
# and returns:
# >P08010|GSTM2_RAT
# 3	-	81	GST_N~1
# 105	-	189	GST_C~2

# This version has been re-written in python to use the EBI/InterPro API parsing json
# to get Pfam protein information:  https://www.ebi.ac.uk/interpro/api/entry/pfam/protein/uniprot/P09488
# to get Pfam domain information: https://www.ebi.ac.uk/interpro/api/entry/pfam/PF02798
#
# currently, it does not provide clan information, because the
# EBI/Interpro/Pfam API does not provide clan information
#
# to get clan (set) information:
# https://www.ebi.ac.uk/interpro/api/set/pfam/entry/pfam/PF02798
#
# which provides (among other things):
# "results": {"metadata": {
#                "accession": "CL0172",
#                "name": "Thioredoxin",
#                "source_database": "pfam"}, ... }

import fileinput
import sys
import re
import json
import argparse
import urllib.request
import urllib.error

interpro_pfam_prot_url = "https://www.ebi.ac.uk/interpro/api/entry/pfam/protein/uniprot/"
interpro_domain_url = "https://www.ebi.ac.uk/interpro/api/entry/pfam/"
interpro_clan_url = "https://www.ebi.ac.uk/interpro/api/set/pfam/entry/pfam/"

def get_pfam_id_www( acc):

    if (acc == 'NODOM'):
        return 'NODOM'

    try:
        req = urllib.request.urlopen(interpro_domain_url + acc)

    except urllib.error.URLError as e:
        prot_info = ''
        sys.stderr.write(e.read().decode('utf-8')+'\n')

    else:
        prot_info = req.read().decode('utf-8')

    json_info= json.loads(prot_info)

## json looks like:
##  'results': [
##      {'metadata': {'accession': 'PF02798', 'name': 'Glutathione S-transferase, N-terminal domain', 'source_database': 'pfam', 'type': 'domain', 'integrated': 'IPR004045', 'member_databases': None, 'go_terms': None}, 
##       'proteins': [
##           {'accession': 'p09488', 
##            'protein_length': 218,
##            'source_database': 'reviewed',
##            'organism': '9606',
##            'entry_protein_locations': [
##                {'fragments': [{'start': 5,'end': 82,'dc-status': 'CONTINUOUS'}],
##                 'model': 'PF02798',
##                 'score': 5.2e-21}]
##            }
##       ]

    pf_dom_id  = json_info['metadata']['name']['short']

    return(pf_dom_id)

def get_clan_info_www(pf_acc):

    if (pf_acc=='NODOM'):
        return {'name':'NODOM','accession':'NODOM'}

    try:
        req = urllib.request.urlopen(interpro_clan_url + pf_acc)

    except urllib.error.URLError as e:
        prot_info = ''
        sys.stderr.write(e.read().decode('utf-8')+'\n')

    else:
        clan_info = req.read().decode('utf-8')


    json_info= json.loads(clan_info)

## results look like:
## "results": [
##      {
##          "metadata": {
##              "accession": "CL0172",
##              "name": "Thioredoxin",
##              "source_database": "pfam"
##          },
## ...
##      }
## ]
##

    if ('results') in json_info:
        return(json_info['results'][0]['metadata'])
    else:
        return None

def get_seq_acc(seq_id):

    if (re.search(r'^gi\|',seq_id)):
        fields = seq_id.split('|')
        acc = fields[3]

    elif (re.search(r'^(sp|tr|up)\|', seq_id)):
        fields = seq_id.split('|')
        (sdb, acc) = fields[0:2]
    else:
      acc = re.split(r'\s',seq_id)[0]

    acc = re.sub(r'\.\d+$','',acc)

    return acc

def get_pfam_www(acc):

    try:
        req = urllib.request.urlopen(interpro_pfam_prot_url + acc)

    except urllib.error.URLError as e:
        prot_info = ''
        sys.stderr.write(e.read().decode('utf-8')+'\n')
        return ([], 0)

    else:
        prot_info = req.read().decode('utf-8')

        if (len(prot_info) == 0):
            return ([],0)

    json_info= json.loads(prot_info)

    pf_dom_list = []

    prot_len = json_info['results'][0]['proteins'][0]['protein_length']

    for result in json_info['results']:

        pfam_info = result['metadata']
        pf_acc = pfam_info['accession']
        pf_name = pfam_info['name']

        for protein in result['proteins']:
            for entry in protein['entry_protein_locations']:
                for frag in entry['fragments']:
                    pf_dom_list.append({'pf_acc':pf_acc, 'start':frag['start'], 'end':frag['end'], 'name':pf_name})
                
    pf_dom_list.sort(key = lambda x: x['start'])

    return(pf_dom_list, prot_len)

def add_nodoms(pf_dom_list, seq_len, min_nodom=10):

    if (len(pf_dom_list) == 0):
        return pf_dom_list

    prev_dom = {'end':0}

    npf_domains = []

    for curr_dom in pf_dom_list:
        if(curr_dom['start'] - prev_dom['end'] > min_nodom):
            new_dom = {'start':prev_dom['end']+1, 'end':curr_dom['start']-1, 'pf_acc':'NODOM'}
            npf_domains.append(new_dom)
        npf_domains.append(curr_dom)
        prev_dom = {'end':curr_dom['end']}

    if (seq_len - pf_dom_list[-1]['end'] > min_nodom):
        new_dom = {'start':pf_dom_list[-1]['end']+1, 'end':seq_len, 'pf_acc':'NODOM'}
        npf_domains.append(new_dom)

    return npf_domains

def print_doms(seq_id, color_ix, args, dom_colors, dom_names, clan_info):
               
    this_acc = get_seq_acc(seq_id)

    if (args.lav):
        print_fmt = "%d\t%d\t%s~%s"
    else:
        print_fmt = "%d\t-\t%d\t%s~%s"

    (pf_dom_list, prot_len)  = get_pfam_www(this_acc)

    ## add no-doms if requested
    if (args.neg_doms):
        pf_dom_list = add_nodoms(pf_dom_list, prot_len, args.min_nodom)

    for dom in pf_dom_list:

        pf_acc = dom['pf_acc']

        ## check if domain has color number
        if (pf_acc in dom_colors):
            dom_color = dom_colors[pf_acc]
        else:
            dom_color = dom_colors[pf_acc] = str(color_ix)
            color_ix += 1

        ## check if domain_acc has short name
        if (pf_acc not in dom_names):
            pf_id = dom_names[pf_acc] = get_pfam_id_www(pf_acc)
            pf_clan_info = clan_info[pf_acc] = get_clan_info_www(pf_acc)
        else:
            pf_id = dom_names[pf_acc]
            pf_clan_info = clan_info[pf_acc]


        this_clan_info = clan_info[pf_acc]

        ## display id or acc?
        if (this_clan_info and not args.no_clans):
            pf_info = 'C.'+ this_clan_info['name']
        else:
            pf_info = pf_id

        if (args.pfam_acc):
            if (pf_info and not args.no_clans):
                pf_info = this_clan_info['accession']
            else:
                pf_info = pf_acc

        if (args.acc_comment):
            pf_info = "%s{%s}"%(pf_info,pf_acc)

        if (args.bound_comment):
            dom_color = "%d:%d"%(dom['start'],dom['end'])

        print(print_fmt%(dom['start'],dom['end'],pf_info,dom_color))

    return color_ix

def read_print_fd(fd, args):

    dom_colors = {'NODOM':'0'}
    dom_names = {}
    clan_info =  {'NODOM':{'name':'NODOM','accession':'NODOM'}}

    color_ix = 1

    for line in fd:
        line = line.strip('\n')
        seq_id = line.split('\t')[0]

        print(">%s"%(seq_id))
        color_ix = print_doms(seq_id, color_ix, args, dom_colors, dom_names, clan_info)

def main() :

    parser=argparse.ArgumentParser(description='ann_pfam_www.py P12345')
    parser.add_argument('--lav',dest='lav',action='store_true',default=False)
    parser.add_argument('--acc_comment',dest='acc_comment',action='store_true',default=False)
    parser.add_argument('--bound_comment',dest='bound_comment',action='store_true',default=False)
    parser.add_argument('--neg_doms',dest='neg_doms',action='store_true',default=False)
    parser.add_argument('--neg-doms',dest='neg_doms',action='store_true',default=False)
    parser.add_argument('--no-clans',dest='no_clans',action='store_true',default=False)
    parser.add_argument('--no_clans',dest='no_clans',action='store_true',default=False)
    parser.add_argument('--min_nodom',dest='min_nodom',action='store',default=10)
    parser.add_argument('--no_over',dest='no_over',action='store_true',default=False)
    parser.add_argument('--no-over',dest='no_over',action='store_true',default=False)
    parser.add_argument('--pfacc',dest='pfam_acc',action='store_true',default=False)
    parser.add_argument('--pfam_acc',dest='pfam_acc',action='store_true',default=False)

    parser.add_argument('files', metavar='FILE', help='files to read, stdin if empty', nargs='*')
    args=parser.parse_args()

    ## check for stdin input
    if (len(args.files) == 0):
        read_print_fd(sys.stdin, args)
        return

    try:
        fd = open(args.files[0],'r')
        read_print_fd(fd, args)
        fd.close()

    except:
        color_ix = 1
        dom_colors = {'NODOM':'0'}
        dom_names = {}
        clan_info = {'NODOM':{'name':'NODOM','accession':'NODOM'}}

        for seq_id in args.files:
            print(">%s"%(seq_id))
            print_doms(seq_id, color_ix, args, dom_colors, dom_names, clan_info)

if __name__ == '__main__':
    main()
