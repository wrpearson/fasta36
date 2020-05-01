#!/usr/bin/env python3

################################################################
# copyright (c) 2017,2018 by William R. Pearson and The Rector &
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

################################################################
# annot_blast_btop4.py --query query.file --ann_script ann_pfam_www.pl --include_doms blast_tab_btop_file
################################################################
# annot_blast_btop4.py associates domain annotation information and
# subalignment scores with a blast tabular (-outfmt 6 or -outfmt 7)
# file that contains the raw score and the BTOP alignment encoding
# This file can be generated from "blastp/n" or "blast_formatter"
# using the command:
#   blast_formatter -archive blast_output.asn -outfmt '7 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore score btop'  > blast_output.tab_annot
#
# If the BTOP field or query_file is not available, the script
# produces domain content without sub-alignment scores.
################################################################
## 2-Dec-2019
# added --have_qslen, --raw_score/--no_raw_score
# made more robust to multiple HSPs when using --ann_file
#
################################################################
## 4-Nov-2018
# add --include_doms, which adds a new field with the coordinates of
# the domains in the protein (independent of alignment)
#
################################################################
## 21-July-2018
# include sequence length (actually alignment end) to produce NODOM's (no NODOM's without length).
#
################################################################
## 13-Jan-2017
# modified to provide query/subject coordinates and identities if no
# query sequence -- does not decrement for reverse-complement fastx/blastx DNA
################################################################
## 16-Nov-2015
# modify to allow multi-query blast searches
################################################################
## 19-Dec-2015
# add -q_annot_script to annotate query sequence
#

import argparse
import fileinput
import sys
import re
import shutil
import subprocess
from math import log

# read lines of the form:
# gi|121694|sp|P20432.1|GSTT1_DROME	gi|121694|sp|P20432|GSTT1_DROME	100.00	209	0	0	1	209	1	209	6e-156	433	1113	209
# gi|121694|sp|P20432.1|GSTT1_DROME	gi|1170090|sp|P04907|GSTF3_MAIZE	26.77	198	123	7	4	185	6	197	2e-08	51.2	121	FL1YG ... 1NKRA1YW1
# gi|121694|sp|P20432.1|GSTT1_DROME	gi|81174731|sp|P0ACA5|SSPA_ECO57	39.66	58	32	2	43	100	49	103	8e-06	43.9	102	EDFLLI ... V-I-NEQS3FM
# gi|121694|sp|P20432.1|GSTT1_DROME	gi|121695|sp|P12653|GSTF1_MAIZE	27.62	181	107	7	32	203	34	199	9e-05	40.8	94	LI1LF ... N-1AS1CLLM1

# and report the domain content ala -m 8CC

def init_blosum62():

    # ncbi_blaa -- list of amino acids
    ncbi_blaa = "A R N D C Q E G H I L K M F P S T W Y V B Z X *".split(' ')

    # blosum62: 2D dict of scoring matrix values

    blosum62 = {}
    blosum62['A'] = dict(zip(ncbi_blaa,[ 4,-1,-2,-2, 0,-1,-1, 0,-2,-1,-1,-1,-1,-2,-1, 1, 0,-3,-2, 0,-2,-1, 0,-4]))
    blosum62['R'] = dict(zip(ncbi_blaa,[-1, 5, 0,-2,-3, 1, 0,-2, 0,-3,-2, 2,-1,-3,-2,-1,-1,-3,-2,-3,-1, 0,-1,-4]))
    blosum62['N'] = dict(zip(ncbi_blaa,[-2, 0, 6, 1,-3, 0, 0, 0, 1,-3,-3, 0,-2,-3,-2, 1, 0,-4,-2,-3, 3, 0,-1,-4]))
    blosum62['D'] = dict(zip(ncbi_blaa,[-2,-2, 1, 6,-3, 0, 2,-1,-1,-3,-4,-1,-3,-3,-1,0,-1,-4,-3,-3,4,1,-1,-4]))
    blosum62['C'] = dict(zip(ncbi_blaa,[ 0,-3,-3,-3, 9,-3,-4,-3,-3,-1,-1,-3,-1,-2,-3,-1,-1,-2,-2,-1,-3,-3,-2,-4]))
    blosum62['Q'] = dict(zip(ncbi_blaa,[-1, 1, 0, 0,-3, 5, 2,-2, 0,-3,-2,1,0,-3,-1,0,-1,-2,-1,-2,0,3,-1,-4]))
    blosum62['E'] = dict(zip(ncbi_blaa,[-1, 0, 0, 2,-4, 2, 5,-2, 0,-3,-3,1,-2,-3,-1,0,-1,-3,-2,-2,1,4,-1,-4]))
    blosum62['G'] = dict(zip(ncbi_blaa,[ 0,-2, 0,-1,-3,-2,-2, 6,-2,-4,-4,-2,-3,-3,-2,0,-2,-2,-3,-3,-1,-2,-1,-4]))
    blosum62['H'] = dict(zip(ncbi_blaa,[-2, 0, 1,-1,-3, 0, 0,-2, 8,-3,-3,-1,-2,-1,-2,-1,-2,-2,2,-3,0,0,-1,-4]))
    blosum62['I'] = dict(zip(ncbi_blaa,[-1,-3,-3,-3,-1,-3,-3,-4,-3,4,2,-3,1,0,-3,-2,-1,-3,-1,3,-3,-3,-1,-4]))
    blosum62['L'] = dict(zip(ncbi_blaa,[-1,-2,-3,-4,-1,-2,-3,-4,-3,2,4,-2,2,0,-3,-2,-1,-2,-1,1,-4,-3,-1,-4]))
    blosum62['K'] = dict(zip(ncbi_blaa,[-1, 2, 0,-1,-3, 1, 1,-2,-1,-3,-2,5,-1,-3,-1,0,-1,-3,-2,-2,0,1,-1,-4]))
    blosum62['M'] = dict(zip(ncbi_blaa,[-1,-1,-2,-3,-1, 0,-2,-3,-2,1,2,-1,5,0,-2,-1,-1,-1,-1,1,-3,-1,-1,-4]))
    blosum62['F'] = dict(zip(ncbi_blaa,[-2,-3,-3,-3,-2,-3,-3,-3,-1,0,0,-3,0,6,-4,-2,-2,1,3,-1,-3,-3,-1,-4]))
    blosum62['P'] = dict(zip(ncbi_blaa,[-1,-2,-2,-1,-3,-1,-1,-2,-2,-3,-3,-1,-2,-4,7,-1,-1,-4,-3,-2,-2,-1,-2,-4]))
    blosum62['S'] = dict(zip(ncbi_blaa,[ 1,-1, 1, 0,-1, 0, 0, 0,-1,-2,-2,0,-1,-2,-1,4,1,-3,-2,-2,0,0,0,-4]))
    blosum62['T'] = dict(zip(ncbi_blaa,[ 0,-1, 0,-1,-1,-1,-1,-2,-2,-1,-1,-1,-1,-2,-1,1,5,-2,-2,0,-1,-1,0,-4]))
    blosum62['W'] = dict(zip(ncbi_blaa,[-3 -3,-4,-4,-2,-2,-3,-2,-2,-3,-2,-3,-1,1,-4,-3,-2,11,2,-3,-4,-3,-2,-4]))
    blosum62['Y'] = dict(zip(ncbi_blaa,[-2,-2,-2,-3,-2,-1,-2,-3, 2,-1,-1,-2,-1,3,-3,-2,-2,2,7,-1,-3,-2,-1,-4]))
    blosum62['V'] = dict(zip(ncbi_blaa,[ 0,-3,-3,-3,-1,-2,-2,-3,-3,3,1,-2,1,-1,-2,-2,0,-3,-1,4,-3,-2,-1,-4]))
    blosum62['B'] = dict(zip(ncbi_blaa,[-2,-1, 3, 4,-3, 0, 1,-1, 0,-3,-4,0,-3,-3,-2,0,-1,-4,-3,-3, 4, 1,-1,-4]))
    blosum62['Z'] = dict(zip(ncbi_blaa,[-1, 0, 0, 1,-3, 3, 4,-2, 0,-3,-3,1,-1,-3,-1,0,-1,-3,-2,-2,1,4,-1,-4]))
    blosum62['X'] = dict(zip(ncbi_blaa,[ 0,-1,-1,-1,-2,-1,-1,-1,-1,-1,-1,-1,-1,-1,-2,0,0,-2,-1,-1,-1,-1,-1,-4]))
    blosum62['*'] = dict(zip(ncbi_blaa,[-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,1]))

    if (len(blosum62.keys()) != len(ncbi_blaa)):
        sys.stderr.write(" blosum62 length mismatch %d != %d\n" %(len(blosum62), len(ncbi_blaa)))
        print(' '.join(ncbi_blaa),file=sys.stderr)
        print(' '.join(blosum62.keys()),file=sys.stderr)

        exit(1)
        
    blosum62_diag = {x:blosum62[x][x] for x in ncbi_blaa}

    return (blosum62, blosum62_diag, -11, -1)

################
# read_annots (\@hit_list)
# input: hit_entry['s_seq_id, etc'], target
# output: modified $hit_entry['domains']
#         modified $hit_entry['sites']
#
# extend to make robust to multiple hits on the same subject

def read_annots(Reader):

    target_set = {}

    current_domain = ""
    hit_ix = 0
    seq_domains = []
    seq_sites = []

    subj_domains = {}

    for line in Reader:
        if (line[0]=='='):
            continue

        line = line.strip("\n")

        # check for header
        if (line[0] == '>'):
            if (current_domain):  # previous domains/sites have already been found and parsed
                if (current_domain not in  target_set):
                    target_set[current_domain] = {}
                    target_set[current_domain]['domains'] = [ d for d in seq_domains ]  # previous domains
                    target_set[current_domain]['sites'] = [ s for s in seq_sites ]      # previous sites
                else:
                    sys.stderr.write("*** phase error: %s duplicate\n"%(current_domain))

            seq_domains = [];   # current domains
            seq_sites = [];     # current sites
            current_domain = line.split(' ')[0][1:]

        else:			# check for data
            a_fields = line.split('\t')
            a_fields[0]=int(a_fields[0])
            if (a_fields[1] == '-'):
                a_fields[2]=int(a_fields[2])
                annot_info = dict(zip(('d_pos','type','d_end','descr'), a_fields))
                re_df=re.compile(r' :(\d+)$')
                annot_info['descr'] =  re_df.sub(r'~\1',annot_info['descr'])
                seq_domains.append(annot_info)
            else:
                annot_info = dict(zip(('d_pos','type', 'd_val', 'descr', a_fields)))
                annot_info['d_end'] = annot_info['d_pos']
                seq_sites.append(annot_info)
    
    Reader.close()

    # get the last one
    if (current_domain):  # previous domains/sites have already been found and parsed
        if (current_domain not in  target_set):
            target_set[current_domain] = {}
            target_set[current_domain]['domains'] = [ d for d in seq_domains ]  # previous domains
            target_set[current_domain]['sites'] = [ s for s in seq_sites ]      # previous sites
        # else:
            # sys.stderr.write("*** phase error: %s duplicate\n"%(current_domain))

    return target_set

################
# merge_annots(hit_r):
#
# take different annotations in hit_r and put them in one list
#
def merge_annots(hit_r):
    merged_annots = []
    if ('q_aligned_domains' in hit_r):
        for annot in hit_r['q_aligned_domains']:
            annot['target']=1
            merged_annots.append(annot)

    if ('aligned_domains' in hit_r):
        for annot in hit_r['aligned_domains']:
            annot['target']=0
            merged_annots.append(annot)

    merged_annots = sorted(merged_annots, key=lambda x: x['qa_start'])

    return(merged_annots)


################
# get_file_annots(file_name)
#
def get_file_annots(file_name, hit_list):

    with open(file_name,'r') as Reader:
        ann_set = read_annots(Reader)

    return ann_set

################
# get_script_annots(script_name, hit_list)
#
# set up stdin/stdout pipe to send in hit list info and read results
#
def get_script_annots(script_name, hit_list, key_list):

    seq_set = {}

    proc = subprocess.Popen(script_name, stdin=subprocess.PIPE, stdout=subprocess.PIPE, shell=True, encoding='utf-8')

    for hit in hit_list:
        (seq_id, seq_len) = (hit[key_list[0]],hit[key_list[1]])
        if (seq_id not in seq_set):
            proc.stdin.write(("%s\t%s\n"%(seq_id,seq_len)))

    proc.stdin.close()

    while (proc.returncode is None):
        proc.poll()

    return read_annots(proc.stdout)

################
#
# link_annots(hit_list, annot_set)
#
# put 'domains' and 'sites' into each hit in the hit list
#
def link_annots(hit_list, annot_set):

    for hit in hit_list:
        seqid = hit['s_seq_id']
        if (seqid in annot_set):
            if ('domains' in annot_set[seqid]):
                hit['domains']=annot_set[seqid]['domains']
            if ('sites' in annot_set[seqid]):
                hit['sites']=annot_set[seqid]['sites']

# input: a blast BTOP string of the form: "1VA160TS7KG10RK27"
# returns a list_ref of tokens: (1, "VA", 60, "TS", 7, "KG, 10, "RK", 27)
#
def decode_btop(btop_str):

    tokens = re.split(r'(\d+)',btop_str)  # split with capture returns both strings between and separator (\d+)

    if not tokens[0]:
        tokens = tokens[1:]

    out_tokens = []

    for token in tokens:
        if re.match(r'\d+',token):
            out_tokens.append(token)
        else:
            mis_tokens = re.split(r'(..)',token)  # split with capture
            for mis in mis_tokens:
                if (mis):
                    out_tokens.append(mis)

    return out_tokens

def parse_query_lib(query_file):

    query_seqs = {}

    with open(query_file,"r") as qfd:
        header=''
        seq_data=''

        for line in qfd:
            line = line.strip("\n")
            
            if (line[0]=='>'):
                if (header):
                    # save existing sequence
                    seq_data = '' + seq_data
                    query_seqs[header]=seq_data

                header = line[1:].split(' ')[0]
            else:
                line = re.sub(r'[^A-Za-z]','',line)
                seq_data += line.upper()

    # save last entry
    if (header):
        query_seqs[header]=seq_data

    return query_seqs

# given: (1) a query sequence; (2) an encoded alignment; (3) a scoring matrix
# calculate:
# (1) the overall score
# (2) a per residue dictionary of scores and mappings from query -> subject
# (2) a per residue dictionary of scores and mappings from subject -> query
#
def alignment_score(query_r, hit, matrix_2d, matrix_diag, g_open, g_ext):

    query_start, subj_start = (int(hit['q_start']),int(hit['s_start']))
    btop_align_r = decode_btop(hit['BTOP'])
    hit['btop_align'] = btop_align_r

    q_map = []
    s_map = []

    gap0, gap1 = (0 ,0)

    q_ix = query_start - 1 # start from zero
    s_ix = subj_start - 1

    score, m_score = (0, 0)
    seq0, seq1 = ("","")

    for btop in btop_align_r:
        if (re.search(r'^\d+$',btop)):   # matching query sequence, add it up
            for i in  range(0,int(btop)):
                res = query_r[q_ix]
                score += matrix_diag[res]
                q_map.append({'s':score, 'y_ix':s_ix, 'res':res})
                s_map.append({'s':score, 'y_ix':q_ix, 'res':res})
                q_ix += 1
                s_ix += 1
        else:
            seq0, seq1 = (btop[0],btop[1])

            if (re.search(r'\-',btop)):  # is there a gap?
                if (seq0 == '-'):      # is it in query?
                    if gap0:             # are we in a gap?
                        score += g_ext
                    else:
                        score += g_open+g_ext
                        gap0 = True
                    # 'y_ix':-1 indicates alignment to gap
                    s_map.append({'s':score, 'y_ix':-1, 'res':seq1})
                    s_ix += 1
                else:			# gap is in subject
                    if gap1:
                        score += g_ext
                    else:
                        score += g_open+g_ext
                        gap1 = True
                    # 'y_ix':-1 indicates alignment to gap
                    q_map.append({'s':score, 'y_ix':-1, 'res':seq0})
                    q_ix += 1
            else:			# mismatch, not gap
                score += matrix_2d[seq0][seq1]
                gap1=gap0 = False
            
                q_map.append({'s':score, 'y_ix':s_ix, 'res':seq0})
                s_map.append({'s':score, 'y_ix':q_ix, 'res':seq1})

                q_ix += 1
                s_ix += 1

    return score, q_map, s_map

################################################################
# sub_alignment_stats -- calculate stats for ONE domain entry
# given x_map, xa_start, xa_end  where x=q/s depending on target
# domain_r : domain boundaries
# 
# calculate a score, identity and boundaries in both sequences and return values
#

def one_sub_alignment_stats(domain_r, x_map, y_map, xa_start, xa_end, ya_start, ya_end):

    td_start, td_end = (domain_r['d_pos'],domain_r['d_end'])

    if (td_end < xa_start or td_start > xa_end):
        return 0

    if (td_start < xa_start):
        td_start = xa_start

    if (td_end > xa_end):
        td_end = xa_end

    td_start -= xa_start
    td_end -= xa_start

    left_score = 0
    if (td_start>0) :
        left_score = x_map[td_start-1]['s']

    score = x_map[td_end]['s'] - left_score

    # map[] coordinates are 0-based
    # ya_start = x_map[td_start]['y_ix']+1
    # ya_end = x_map[td_end]['y_ix']+1

    #### identity calculation:
    n_len = 0
    n_id = 0
    for xi in range(td_start, td_end+1):
        this_x = x_map[xi]

        x_res=this_x['res']
        if (this_x['y_ix'] >= 0):
            n_len += 1
            y_res = y_map[this_x['y_ix']-ya_start+1]['res']
            if (x_res.upper() == y_res.upper()):
                n_id += 1

    ident = float(n_id)/float(n_len)

    return score, ident, td_start+xa_start-1, td_end+xa_start-1, x_map[td_start]['y_ix'], x_map[td_end]['y_ix']

################
# get domain scores, idents, boundaries for list of domains
#
def do_sub_alignment_stats(domain_list, x_map, y_map, xa_start, xa_end, ya_start, ya_end, keys_str):

    aligned_doms = []

    for domain in domain_list:
        subalign_data = one_sub_alignment_stats(domain, x_map, y_map, xa_start, xa_end, ya_start, ya_end)
        if (subalign_data and len(subalign_data)==6):
            sub_data = dict(zip(keys_str,subalign_data))
            for k in ('type','descr'):
                sub_data[k] = domain[k]
            aligned_doms.append(sub_data)

    return(aligned_doms)

####
# print raw domain info:
# |DX:%d-%d;C=dom_info|XD:%d-%d:C=dom_info
#
def format_dom_info(q_dom_r, dom_r):

    dom_str = ""
    for dom  in q_dom_r:
        dom_str += "|DX:%d-%d;C=%s"%(dom['d_pos'],dom['d_end'], dom['descr'])

    for dom  in dom_r:
        dom_str += "|XD:%d-%d;C=%s"%(dom['d_pos'],dom['d_end'], dom['descr'])

    return dom_str


def format_annot_info(annot_list_r, hit):

    annot_str = "";

    # two types of annotations, domains and sites.

    score_scale = hit['score']/hit['raw_score']

    for annot_r in (annot_list_r ):

        if (annot_r['type'] == '-'):
            fsub_score = annot_r['score']/hit['raw_score']

            ns_score, s_bit = (int(annot_r['score'] * score_scale + 0.5),
                               fsub_score * hit['bits'])

            qval = 0.0
            if (hit['evalue'] == 0.0):
                if (s_bit > 50.0):
                    qval = 3000.0
                else:
                    qval = -10.0 * (2.0*log(400.0) + s_bit)/log(10.0)
            else:
                qval = -10.0*log(hit['evalue'])*fsub_score/log(10.0)

            if qval < 0.0:
                qval = 0.0

            rx_str = 'XR'
            if (annot_r['target']):
                rx_str = "RX"
            annot_str += ';'.join(("|%s:%d-%d:%d-%d:s=%d"%(rx_str,
                                      annot_r['qa_start']+1,annot_r['qa_end']+1,
                                      annot_r['sa_start']+1,annot_r['sa_end']+1,ns_score),
                                   "b=%.1f"%(s_bit),"I=%.3f"%(annot_r['ident']),
                                   "Q=%.1f"%(qval),"C=%s"%(annot_r['descr'])))

        else:        # site annotation
            ann_type = annot_r['type'];
            site_str = "|%cX"%(ann_type)
            if (annot_r['target'] == 1):
                site_str = "|X%c"%(ann_type)
            elif (annot_r['target'] == 2):
                site_str = "|%c%c"%(ann_type, ann_type)

            annot_str += "%s:"%(site_str)
            annot_str += "%d%s%s%d%s"%(annot_r['qa_pos'], annot_r['q_res'], annot_r['m_symb'],
                                        annot_r['sa_pos'], annot_r['s_res'])

    return annot_str

def main(args):

    blosum62, blosum62_diag, g_open, g_ext = init_blosum62()

    if (args.query_file):
        # query_lib_r has a set of query sequences
        query_lib_r = parse_query_lib(args.query_file)
    else:
        sys.stderr.write("--query required\n")
        exit(1)

    tab_fields = "q_seqid s_seqid percid alen mismatch gopen q_start q_end s_start s_end evalue bits BTOP".split(' ')
    int_fields = "alen mismatch gopen q_start q_end s_start s_end".split(' ')
    float_fields = "percid evalue bits score".split(' ')

    if (args.have_qslen):
        tab_fields = "q_seqid q_len s_seqid s_len percid alen mismatch gopen q_start q_end s_start s_end evalue bits BTOP".split(' ')
        int_fields = "q_len s_len alen mismatch gopen q_start q_end s_start s_end".split(' ')

    # the fields that are displayed are listed here.  By default, all fields except score and BTOP are displayed.
    out_tab_fields = tab_fields[0:-1]
    in_tab_fields = tab_fields[0:-1]

    if (args.raw_out):
        out_tab_fields.append("raw_score")

    if (args.raw_in):
        in_tab_fields.append("score")

    ## always add BTOP
    in_tab_fields.append("BTOP")
    tab_fields = in_tab_fields

    if (args.out_fields):
        out_tab_fields = out_fields.split(" ")

    header_lines = []

    next_line = ""
    have_data = False

    hit_list = []
    q_hit_list = []

    for line in fileinput.input(args.files):
        if (line[0] == '#'):
            if (have_data):
                next_line = line
                have_data = False
                break
            else:
                header_lines.append(line)
            continue

        have_data = True
        line = line.strip('\n')
        if (line):
            this_data = dict(zip(tab_fields, line.split("\t")))
            for k in this_data.keys():
                if (k in int_fields):
                    this_data[k] = int(this_data[k])
                if (k in float_fields):
                    this_data[k] = float(this_data[k])
            hit_list.append(this_data)

    # get the query annotations
    q_hit_list = []
    if (args.q_ann_file):
        q_seqid = hit_list[0]['q_seqid']
        q_hit_list.append({'s_seq_id':q_seqid, 's_end':len(query_lib_r[q_seqid])})
        q_annots = get_file_annots(args.q_ann_file, q_hit_list)
        link_annots(q_hit_list, q_annots)

    elif (args.q_ann_script):
        args.q_ann_script = re.sub(r'\+',' ',args.q_ann_script)

        if (args.q_ann_script and shutil.which(args.q_ann_script.split(" ")[0])):
            q_seqid = hit_list[0]['q_seqid']
            q_hit_list.append({'s_seq_id':q_seqid, 's_end':len(query_lib_r[q_seqid])})
            q_annots = get_script_annots(args.q_ann_script, q_hit_list, ['s_seq_id','s_end'])
            link_annots(q_hit_list, q_annots)

    # get the subject annotations
    # first set up the list with sequence lengths
    if (args.ann_file or args.ann_script):
        s_len = 100000
        for hit in hit_list:
            hit['s_seq_id']=hit['s_seqid']
            if (not args.have_qslen):
                hit['s_end']=s_len
        
        if (args.ann_file):
            s_annots = get_file_annots(args.ann_file, hit_list)
            link_annots(hit_list, s_annots)
        elif (args.ann_script):
            args.ann_script = re.sub(r'\+',' ',args.ann_script)
            if (shutil.which(args.ann_script.split(" ")[0])):
                s_annots = get_script_annots(args.ann_script, hit_list,['s_seq_id','s_end'])
                link_annots(hit_list, s_annots)

    for line in header_lines:
        print(line, end='')

    header_lines = [next_line]

    # now get query annotation if available

    for hit in hit_list:
        list_covered = []

        # If I have an encoded aligment {BTOP} and a query sequence query_lib_r && query_lib_r[hit['q_seqid']]
        # then I can calculate sub-alignment scores
        if ('BTOP' in hit and query_lib_r and hit['q_seqid'] in query_lib_r):

            # calculate raw_score and mappings
            hit['raw_score'], q_map, s_map = alignment_score(query_lib_r[hit['q_seqid']],
                                                             hit,blosum62, blosum62_diag, g_open, g_ext)
            if ('score' not in hit):
                hit['score'] = hit['raw_score']

            # calculate sub-alignment scores in subject/library coordinates

            if ('domains' in hit and len(hit['domains'])>0):
                hit['aligned_domains'] = do_sub_alignment_stats(hit['domains'], s_map, q_map, hit['s_start'],hit['s_end'],hit['q_start'],hit['q_end'],
                                                                ('score','ident','sa_start', 'sa_end', 'qa_start', 'qa_end'))

            # calculate sub-alignment scores in query coordinates
            if (len(q_hit_list) > 0 and 'domains' in q_hit_list[0] and len(q_hit_list[0]['domains'])>0):
                hit['q_aligned_domains'] = do_sub_alignment_stats(q_hit_list[0]['domains'], q_map, s_map, hit['q_start'],hit['q_end'],hit['s_start'],hit['s_end'],
                                                                  ('score','ident','qa_start', 'qa_end', 'sa_start', 'sa_end'))
        ################
        ## final output display

        print("\t".join([str(hit[x]) for x in out_tab_fields]),end='') # show fields from original blast tabular file

        merged_annots_r = merge_annots(hit)                # merge the four possible annotation lists into one.

        if (len(merged_annots_r)>0):
            print("\t"+format_annot_info(merged_annots_r, hit),end='')
            if (args.dom_info):
                if (len(q_hit_list) > 0 and 'domains' in q_hit_list[0]):
                    print("\t"+format_dom_info(q_hit_list[0]['domains'], hit['domains']),end='')
                else:
                    print("\t"+format_dom_info([], hit['domains']),end='')
        elif (len(list_covered)>0):
            print("\t" + ";".join(list_covered))
            if (args.dom_info):
                print("\t"+format_dom_info(q_hit_list[0]['domains'], hit['domains']),end='')
        print()

    for line in header_lines:
        print(line,end="")


if __name__ == '__main__':

    print('# ' + ' '.join(sys.argv))

    parser=argparse.ArgumentParser(description='annot_blast_btop4.py : annotate blast tabular format with BTOP ')
    # not implemented
    # parser.add_argument('--matrix', help='scoring matrix',dest='matrix',action='store',default='BL62')
    parser.add_argument('--ann_script', help='script for subject annotations',dest='ann_script',action='store')
    parser.add_argument('--q_ann_script', help='script for query annotations',dest='q_ann_script',action='store')
    parser.add_argument('--ann_file', help='subject annotation file',dest='ann_file',action='store')
    parser.add_argument('--q_ann_file', help='query annotation file',dest='q_ann_file',action='store')
    parser.add_argument('--have_qslen', help='query/subject lenghts in tab file',dest='have_qslen',action='store_true',default=False)
    parser.add_argument('--dom_info', help='show unaligned domain coordinates',dest='dom_info',action='store_true',default=False)
    parser.add_argument('--sub2query', help='get query annots from self-subject',dest='sub_query',action='store_true',default=False)
    parser.add_argument('--query', help='file of query sequences',dest='query_file',action='store')
    parser.add_argument('--out_fields', help='names/order of output fields',dest='out_fields',action='store')
    parser.add_argument('--raw_score', help='raw score after bit score',dest='raw_in',action='store_true',default=True)
    parser.add_argument('--no_raw_score', help='raw score after bit score',dest='raw_in',action='store_false', default=True)
    parser.add_argument('--no-raw_score', help='raw score after bit score',dest='raw_in',action='store_false', default=True)
    parser.add_argument('--raw_score_out', help='display raw score',dest='raw_out',action='store_true',default=False)
    parser.add_argument('files', metavar='FILE', help='Blast tabular BTOP files to read', nargs='*')

    args=parser.parse_args()

    main(args)
