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

# ann_pfam_www0.py takes an annotation file from fasta36 -V with a line of the form:

# sp|P0810|GSTM1_RAT [tab] seqlen
#
# and returns:
# >P08010|GSTM2_RAT
# 3	-	81	GST_N~1
# 105	-	189	GST_C~2

# This version has been re-written in python from ann_pfam_sql.pl
# by default, it shows clan assignments (not possible with ann_pfam_www.py)
# currently it does not implement virtual domains

import fileinput
import sys
import re
import json
import argparse
import mysql.connector

## import urllib.request
## import urllib.error

sql_get_pfam_acc = '''
SELECT seq_start, seq_end, model_start, model_end, model_length, pfamA_acc, pfamA_id, auto_pfamA_reg_full, domain_evalue_score as evalue, length
FROM pfamseq
JOIN pfamA_reg_full_significant using(pfamseq_acc)
JOIN pfamA USING (pfamA_acc)
WHERE in_full = 1
AND  pfamseq_acc="%s"
ORDER BY seq_start
'''

sql_get_upfam_acc = '''
SELECT seq_start, seq_end, model_start, model_end, model_length, pfamA_acc, pfamA_id, auto_uniprot_reg_full as auto_pfamA_reg_full, domain_evalue_score as evalue, length
FROM uniprot
JOIN uniprot_reg_full using(uniprot_acc)
JOIN pfamA USING (pfamA_acc)
WHERE in_full = 1
AND  uniprot_acc="%s"
ORDER BY seq_start
'''

sql_get_pfam_clan = '''
SELECT clan_acc, clan_id
FROM clan
JOIN clan_membership using(clan_acc)
WHERE pfamA_acc="%s"
'''

def get_seq_acc(seq_id):

    if (re.search(r'^gi\|',seq_id)):
        (tmp, gi, sdb, acc, id) = seq_id.split('|')

    elif (re.search(r'^(sp|tr|up)\|', seq_id)):
        (sdb, acc, id) = seq_id.split('|')
    else:
      acc = re.split(r'\s',seq_id)[0]

    acc = re.sub(r'\.\d+$','',acc)

    return acc

def get_pfam_clan(pf_acc, pf_id, db_cursor):

    db_cursor.execute(sql_get_pfam_clan%(pf_acc,))

    row = db_cursor.fetchone()

    if (not row):
        return (pf_acc, pf_id)
    else:
        return ("C."+row['clan_acc'],"C."+row['clan_id'])
    
def get_pfam_sql(acc, db_cursor):

    pf_dom_list = []
    prot_len = 0

    db_cursor.execute(sql_get_pfam_acc%(acc,))

    for row in db_cursor:
        prot_len = row['length']
        pf_dom_list.append(row)
        
    if (len(pf_dom_list) == 0):
        db_cursor.execute(sql_get_upfam_acc%(acc,))

        for row in db_cursor:
            prot_len = row['length']
            pf_dom_list.append(row)
            
    return(pf_dom_list, prot_len)

def add_nodoms(pf_dom_list, seq_len, min_nodom=10):

    if (len(pf_dom_list) == 0):
        return pf_dom_list

    prev_dom = {'seq_end':0}

    npf_domains = []

    for curr_dom in pf_dom_list:
        if(curr_dom['seq_start'] - prev_dom['seq_end'] > min_nodom):
            new_dom = {'seq_start':prev_dom['seq_end']+1, 'seq_end':curr_dom['seq_start']-1, 
                       'pfamA_acc':'NODOM', 'pfamA_id':'NODOM','clan_acc':'NODOM','clan_id':'NODOM'}
            npf_domains.append(new_dom)
        npf_domains.append(curr_dom)
        prev_dom = {'seq_end':curr_dom['seq_end']}

    if (seq_len - pf_dom_list[-1]['seq_end'] > min_nodom):
        new_dom = {'seq_start':pf_dom_list[-1]['seq_end']+1, 'seq_end':seq_len, 
                   'pfamA_acc':'NODOM', 'pfamA_id':'NODOM','clan_acc':'NODOM','clan_id':'NODOM'}
        npf_domains.append(new_dom)

    return npf_domains

def print_doms(seq_id, color_ix, args, dom_colors, dom_names, db_cursor):
               
    this_acc = get_seq_acc(seq_id)

    if (args.lav):
        ## --lav is missing a '\t-\t' field
        print_fmt = "%d\t%d\t%s~%s"
    else:
        ## default format
        print_fmt = "%d\t-\t%d\t%s~%s"

    (pf_dom_list, prot_len)  = get_pfam_sql(this_acc, db_cursor)

    if (not args.no_clans):
        for dom in pf_dom_list:
            (dom['clan_acc'], dom['clan_id']) = get_pfam_clan(dom['pfamA_acc'],dom['pfamA_id'],db_cursor)

    ## add no-doms if requested
    if (args.neg_doms):
        pf_dom_list = add_nodoms(pf_dom_list, prot_len, args.min_nodom)

    for dom in pf_dom_list:

        pf_acc = dom['pfamA_acc']

        ## check if domain has color number
        if (pf_acc in dom_colors):
            dom_color = dom_colors[pf_acc]
        else:
            dom_color = dom_colors[pf_acc] = str(color_ix)
            color_ix += 1

        ## display id or acc?
        pf_info = dom['clan_id']
        if (args.no_clans):
            pf_info = dom['pfamA_id']

        if (args.pfam_acc):
            pf_info = dom['clan_acc']
            if (args.no_clans):
                pf_info = dom['pfamA_acc']
            
        if (args.acc_comment):
            pf_info = "%s{%s}"%(pf_info,pf_acc)

        if (args.bound_comment):
            dom_color = "%d:%d"%(dom['seq_start'],dom['seq_end'])

        print(print_fmt%(dom['seq_start'],dom['seq_end'],pf_info,dom_color))

    return color_ix

def read_print_fd(fd, args, db_cursor):

    dom_colors = {'NODOM':'0'}
    dom_names = {}

    color_ix = 1

    for line in fd:
        line = line.strip('\n')
        seq_id = line.split('\t')[0]

        print(">%s"%(seq_id))
        color_ix = print_doms(seq_id, color_ix, args, dom_colors, dom_names,db_cursor)

def main() :

    parser=argparse.ArgumentParser(description='ann_pfam_www.py P12345')
    ## db parameters
    parser.add_argument('--host',dest='host',action='store',default='wrpxdb.bioch.virginia.edu')
    parser.add_argument('--user',dest='user',action='store',default='web_user')
    parser.add_argument('--passwd',dest='passwd',action='store',default='fasta_www')
    parser.add_argument('--db',dest='db',action='store',default='pfam35')

    ## domain presentation parameters
    parser.add_argument('--lav',dest='lav',action='store_true',default=False)
    parser.add_argument('--acc_comment',dest='acc_comment',action='store_true',default=False)
    parser.add_argument('--bound_comment',dest='bound_comment',action='store_true',default=False)
    parser.add_argument('--no_clans',dest='no_clans',action='store_true',default=False)
    parser.add_argument('--no-clans',dest='no_clans',action='store_true',default=False)
    parser.add_argument('--neg_doms',dest='neg_doms',action='store_true',default=False)
    parser.add_argument('--neg-doms',dest='neg_doms',action='store_true',default=False)
    parser.add_argument('--min_nodom',dest='min_nodom',action='store',default=10)
    parser.add_argument('--no_over',dest='no_over',action='store_true',default=False)
    parser.add_argument('--no-over',dest='no_over',action='store_true',default=False)
    parser.add_argument('--pfacc',dest='pfam_acc',action='store_true',default=False)
    parser.add_argument('--pfam_acc',dest='pfam_acc',action='store_true',default=False)

    parser.add_argument('files', metavar='FILE', help='files to read, stdin if empty', nargs='*')
    args=parser.parse_args()

    pf_db = mysql.connector.connect(db=args.db, host=args.host, user=args.user, passwd=args.passwd)

    pf_cur = pf_db.cursor(dictionary=True, buffered=True)

    ## check for stdin input
    if (len(args.files) == 0):
        read_print_fd(sys.stdin, args, pf_cur)
        return

    try:
        fd = open(args.files[0],'r')
        read_print_fd(fd, args, pf_cur)
        fd.close()

    except:
        color_ix = 1
        dom_colors = {'NODOM':'0'}
        dom_names = {}

        for seq_id in args.files:
            print(">%s"%(seq_id))
            print_doms(seq_id, color_ix, args, dom_colors, dom_names,pf_cur)

if __name__ == '__main__':
    main()
