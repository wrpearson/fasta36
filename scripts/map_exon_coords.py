#!/usr/bin/env python3
# 
# given a -m8CB file with exon annotations for the query and subject,
# provide a function that maps subject coordinates to query, or vice versa
#
# appropriate exon annotations can be produced by:
#
# fasta36 -q -m 8CBL -V \!ann_exons_up_sql.pl+--gen_coord+--exon_label -V q\!ann_exons_up_sql.pl+--gen_coord+--exon_label \!get_protein.py+P30711 \!get_protein.py+Q2NL00 > hum_v_bov_gstt1.m8CBL
# ann_exons_up_sql.pl --gen_coord --exon_label is required to produce the appropriate exons
#
#
################################################################
# copyright (c) 2018 by William R. Pearson and The Rector &
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

import fileinput
import sys
import re
import argparse

################
# "domain" class that describes a domain/exon alignment annotation
#
class exonInfo:
    def __init__(self, name, q_target, p_start, p_end, chrom, d_start, d_end, full_text):
        self.name = name
        self.q_target = q_target
        self.p_start = p_start
        self.p_end = p_end
        self.chrom = chrom
        self.d_start = d_start
        self.d_end = d_end
        self.text = full_text
        self.plus_strand = True
        if (d_start > d_end):
            self.plus_strand = False

    def __str__(self):
        rxr_str = "XD"
        if (self.q_target):
            rxr_str="DX"
        return '|%s:%i-%i:%s{%s:%i-%i}' % (rxr_str, self.p_start, self.p_end, self.name, self.chrom, self.d_start, self.d_end)

class exonAlign:
    def __init__(self, name, q_target, qp_start, qp_end, sp_start, sp_end, full_text):
        self.exon = None

        self.name = name
        self.q_target = q_target

        self.q_start = qp_start
        self.q_end = qp_end
        self.s_start = sp_start
        self.s_end = sp_end

        self.text = full_text
        self.out_str = ''

    def __str__(self):
        rxr_str = "RX"
        if (self.q_target):
            rxr_str="XR"
        return "[%s:%i-%i:%i-%i::%s" % (rxr_str,self.q_start, self.q_end, self.s_start, self.s_end, self.name)

    def print_bar_str(self):  # checking for 'NADA' 
        if (not self.out_str):
            self.out_str = self.text
        return str("|%s"%(self.out_str))

# Parses domain annotations after split at '|'

#
def parse_exon_align(text):
    # takes a domain in string form, turns it into a domain object
    # looks like: RX:5-82:5-82:s=397;b=163.1;I=1.000;Q=453.6;C=C.Thioredoxin~1
    # could also look like: RX:5-82:5-82:s=397;b=163.1;I=1.000;Q=453.6;C=C.Thioredoxin{PF012445}~1

    # get RX/XR and qstart/qstop sstart/sstop as strings
    m = re.search(r'^(\w+):(\d+)-(\d+):(\d+)-(\d+):',text)
    if (m):
        (RXRState, qstart_s, qend_s, sstart_s, send_s) = m.groups()
    else:
        sys.stderr.write("could not parse exon location: %s\n"%(text))

    # get domain name/color (and possibly {info})

    (name, color_s) = re.search(r';C=([^~]+)~(.+)$',text).groups()
    info_s=""

    if (re.search(r'\}$',name)):
        (name, info_s) = re.search(r'([^\{]+)(\{[^\}]+\})$',name).groups()

    q_target = True
    if (RXRState=='XR'):
        q_target = False

    exon_align = exonAlign(name, q_target, int(qstart_s), int(qend_s), int(sstart_s), int(send_s), 
                           text)

    return exon_align

################
# exon_info is like domain, but no scores
#
def parse_exon_info(text):
    # takes a domain in string form, turns it into a domain object
    # looks like: DX:1-100;C=C.Thioredoxin~1

    (RXRState, start_s, end_s,name, color) = re.search(r'^(\w+):(\d+)-(\d+);C=([^~]+)~(.*)$',text).groups()
    info = ""
    if (re.search(r'\}$',name)):
        (name, info) = re.search(r'([^\{]+)(\{[^\}]+\})$',name).groups()

    gene_re = re.search(r'^\{(\w+):(\d+)\-(\d+)\}',info)
    if (gene_re):
        (chrom, d_start, d_end) = gene_re.groups()
    else:
        sys.stderr.write("genome info not found: %s\n" % (text))

    q_target = True;
    if (RXRState == 'XD'):
        q_target = False

    exon_info = exonInfo(name, q_target, int(start_s), int(end_s), chrom, int(d_start), int(d_end), text)

    return exon_info

####
# parse_protein(result_line)
# takes a protein in string format, turns it into a dictionary properly
# looks like:   sp|P30711|GSTT1_HUMAN   up|Q2NL00|GSTT1_BOVIN   86.67   240     32      0       1       240     1       240     1.4e-123        444.0   16VI7DR6IT3IR15KQ3AI6TI11TA7YH8RC12TA3SN10FL10QETM2AT6VMTA2LV2DG4ND6PS24EK6TA11DV14FSPQ5IL3LMML1WK5RQ   |XR:4-76:4-76:s=327;b=134.6;I=0.895;Q=367.8;C=C.Thioredoxin~1|RX:5-82:5-82:s=356;b=146.5;I=0.902;Q=403.3;C=C.Thioredoxin~1|RX:83-93:83-93:s=52;b=21.4;I=0.818;Q=30.9;C=NODOM~0|XR:77-93:77-93:s=86;b=35.4;I=0.882;Q=72.6;C=NODOM~0|RX:94-110:94-110:s=88;b=36.2;I=0.882;Q=75.0;C=vC.GST_C~2v|XR:94-110:94-110:s=88;b=36.2;I=0.882;Q=75.0;C=vC.GST_C~2v|RX:111-201:111-201:s=409;b=168.3;I=0.868;Q=468.3;C=C.GST_C~2|XR:111-201:111-201:s=409;b=168.3;I=0.868;Q=468.3;C=C.GST_C~2|RX:202-240:202-240:s=154;b=63.4;I=0.795;Q=155.9;C=NODOM~0|XR:202-240:202-240:s=154;b=63.4;I=0.795;Q=155.9;C=NODOM~0
#
def parse_protein(line_data,fields, req_name):
    # last part (domain annotions) split('|') and parsed by parse_domain()

    data = {}
    data = dict(list(zip(fields, line_data)))
    if (re.search(r'\|',data['qseqid'])):
        data['qseq_acc'] = data['qseqid'].split('|')[1]
    else:
        data['qseq_acc'] = data['qseqid']

    if (re.search(r'\|',data['sseqid'])):
        data['sseq_acc'] = data['sseqid'].split('|')[1]
    else:
        data['sseq_acc'] = data['sseqid']

    Qexon_list = []
    Sexon_list = []

    Qinfo_list = []
    Sinfo_list = []

    counter = 0

    if ('align_annot' in data and len(data['align_annot']) > 0):
        for exon_str in data['align_annot'].split('|')[1:]:
            if (req_name and not re.search(req_name, exon_str)):
                continue

            counter += 1
            exon = parse_exon_align(exon_str)
            if (exon.q_target):
                Qexon_list.append(exon)
            else:
                Sexon_list.append(exon)

    data['q_exalign_list'] = Qexon_list
    data['s_exalign_list'] = Sexon_list

    if ('exon_info' in data and len(data['exon_info']) > 0):
        for info_str in data['exon_info'].split('|')[1:]:
            if (not re.search(r'^[DX][XD]',info_str)):
                continue

            dinfo = parse_exon_info(info_str)

            if (dinfo.q_target):
                Qinfo_list.append(dinfo)
            else:
                Sinfo_list.append(dinfo)


        # put links to info_list into exon_list so info_list names can
        # be changed  -- give S/Qinfo's the S/Qdom ids of the overlapping domain

        # find_info_overlaps(Qinfo_list, Qexon_list)
        # find_info_overlaps(Sinfo_list, Sexon_list)
        
    data['q_exinfo_list'] = Qinfo_list
    data['s_exinfo_list'] = Sinfo_list

    return data

################
#
# decode_btop() - 
#   input: a blast BTOP string of the form: "1VA160TS7KG10RK27"
#   returns a list_ref of tokens: (1, "VA", 60, "TS", 7, "KG, 10, "RK", 27)
def decode_btop(btop_str):
    out_tokens = []
    for token in re.split(r'(\d+)',btop_str):
        if (not token): continue
        if re.match(r'\d+',token):
            out_tokens.append(token)
        else:
            for mismat in re.split(r'(..)',token):
                if (mismat): out_tokens.append(mismat)

    return out_tokens

################
#
#  map_align(btop, q_start, s_start)
#    input: btop
#    output: q_pos_arr, s_pos_arr
#
def map_align(btop_str, q_start, s_start):

    q_pos = q_start
    s_pos = s_start

    q_pos_arr = []
    s_pos_arr = []

    btop_tokens = decode_btop(btop_str)

    for t in btop_tokens:
        if (re.match(r'\d+',t)):
            for i in range(int(t)) :
                q_pos_arr.append(q_pos)
                q_pos += 1
                s_pos_arr.append(s_pos)
                s_pos += 1
        elif (re.match(r'\-\w',t)):
            q_pos_arr.append(q_pos)
            s_pos_arr.append(s_pos)
            s_pos += 1
        elif (re.match(r'\w\-',t)):
            q_pos_arr.append(q_pos)
            q_pos += 1
            s_pos_arr.append(s_pos)
        else:
            q_pos_arr.append(q_pos)
            q_pos += 1
            s_pos_arr.append(s_pos)
            s_pos += 1

    return q_pos_arr, s_pos_arr

################
#
# map_coords(from_coords, to_coords, coord_list)
#
def map_coords(from_coords, to_coords, coord_list):

    mapped_coords = []

    fx = 0
    mx = 0
    while mx < len(coord_list):
        this_from_coord = coord_list[mx]
        while (from_coords[fx] < this_from_coord):
            fx += 1
            continue

        mapped_coords.append(to_coords[fx])
        mx += 1

    return mapped_coords
    
################
#
# map_align_coords()  given a BTOP, q_start, s_start, and s_target, generate s_coords for list of q_coords
#
def map_align_coords(btop_str, q_start, s_start, s_target, coord_list):
    
    (q_coords, s_coords) = map_align(btop_str, q_start, s_start)

    sorted_coord_list = sorted(coord_list)

    if (s_target):
        s_mapped_coords = map_coords(q_coords, s_coords, sorted_coord_list)
    else:
        s_mapped_coords = map_coords(s_coords, q_coords, sorted_coord_list)

    coord_dict={}
    for ix, s_coord in enumerate(sorted_coord_list):
        coord_dict[s_coord]=s_mapped_coords[ix]
        
    return [ coord_dict[c] for c in coord_list ]


################
#
# aa_to_exon()  --- given a coordinate and the corresponding exon map, return the exon coordinate
# (can only be done for aligned exons)
#
# this version of the function must use an info_list, not an
# align_list, because it uses p_start/p_end rather than qp_start/sp_start, etc.
# a version using qp_start/sp_start would also need a target argument
#
def aa_to_exon(aa_coords, exon_info_list):

    if len(exon_info_list) < 1:
        return []

    sorted_aa_coords = sorted(aa_coords)

    pos_strand = True
    if (exon_info_list[0].d_start > exon_info_list[0].d_end):
        pos_strand = False

    ex_x = 0
    exon_coords = []

    aap_x = 0
    this_aap = sorted_aa_coords[aap_x]
    while (ex_x < len(exon_info_list)):
        this_exon = exon_info_list[ex_x]
        if (this_aap <= this_exon.p_end and this_aap >= this_exon.p_start):
            aa_dna_offset = (this_aap - this_exon.p_start) * 3 

            if (pos_strand):
                aa_dna_pos = this_exon.d_start + aa_dna_offset
            else:
                aa_dna_pos = this_exon.d_start - aa_dna_offset

            exon_coords.append({'chrom':this_exon.chrom, 'dpos':aa_dna_pos})
            aap_x += 1
            if (aap_x < len(sorted_aa_coords)):
                this_aap = sorted_aa_coords[aap_x]
            else:
                break
        else:
            ex_x += 1

    aa_coord_dict = {}
    for aap_x, aap in enumerate(sorted_aa_coords):
        aa_coord_dict[aap] = exon_coords[aap_x]

    return [aa_coord_dict[ax] for ax in aa_coords]

################
# set_data_fields() -- initialize field[] used to generate data[] dict
#
def set_data_fields(args, line_data) :

    field_str = 'qseqid sseqid pident length mismatch gapopen q_start q_end s_start s_end evalue bitscore BTOP align_annot'
    field_qs_str = 'qseqid q_len sseqid s_len pident length mismatch gapopen q_start q_end s_start s_end evalue bitscore BTOP align_annot'

    if (len(line_data) > 1) :
        if ((not args.have_qslen) and  re.search(r'\d+',line_data[1])):
            args.have_qslen=True

        if ((not args.exon_info) and re.search(r'^\|[DX][XD]\:',line_data[-1])):
            args.exon_info = True

    end_field = -1
    fields = field_str.split(' ')

    if (args.have_qslen):
        fields = field_qs_str.split(' ')

    if (args.exon_info):
        fields.append('exon_info')
        end_field = -2

    return (fields, end_field)

################################################################
#
# main program 

def main():

    print("#"+" ".join(sys.argv))

    data_fields_reset=False

    parser=argparse.ArgumentParser(description='map_exon_coords.py result_file.m8CB saa:coord : map subject coordinate to query genomic coordinate')
    parser.add_argument('--have_qslen', help='bl_tab fields include query/subject lengths',dest='have_qslen',action='store_true',default=False)
    parser.add_argument('--exon_info', help='raw domain coordinates included',action='store_true',default=True)
    parser.add_argument('--subj_aa',help='subject aa coordinate to map',action='store',type=int,dest='subj_aa_coord',default=1)
    parser.add_argument('files', metavar='FILE', nargs='*', help='files to read, if empty, stdin is used')
    args=parser.parse_args()

    end_field = -1
    data_fields_reset=False

    (fields, end_field) = set_data_fields(args, [])

    if (args.have_qslen and args.exon_info):
        data_fields_reset=True

    saved_qexon_list = []
    qexon_list = []

    for line in fileinput.input(args.files):
    # pass through comments
        if (line[0] == '#'):
            print(line, end='')	# ',' because have not stripped
            continue

        ################
        # break up tab fields, check for extra fields
        line = line.strip('\n')
        line_data = line.split('\t')
        if (not data_fields_reset):     # look for --have_qslen number, --exon_info data, even if not set
            (fields, end_field) = set_data_fields(args, line_data)
            data_fields_reset = True

        ################
        # get exon annotations
        # produces:  data['q_exalign_list'], data['s_exalign_list']
        #            data['q_exinfo_list'],  data['s_exinfo_list']
        data = parse_protein(line_data,fields,"exon")	# get score/alignment/domain data

        # extract aligned query_coordinates
        q_coords = []
        sa_from_qa = []
        for q_ex in data['q_exalign_list']:
            q_coords.append(q_ex.q_start)
            q_coords.append(q_ex.q_end)
            sa_from_qa.append(q_ex.s_start)
            sa_from_qa.append(q_ex.s_end)

        s_coords = []
        qa_from_sa = []
        for s_ex in data['s_exalign_list']:
            s_coords.append(s_ex.s_start)
            s_coords.append(s_ex.s_end)
            qa_from_sa.append(s_ex.q_start)
            qa_from_sa.append(s_ex.q_end)

        ################
        # map aligned coordinates in query to subject exons
        # -- this is not necessary -- it already in data['q_exalign_list'].s_start/s_end
        # s_target=True
        # sa_from_qa = map_align_coords(data['BTOP'], int(data['q_start']), int(data['s_start']),
        #                                    s_target, qa_coords)
        sex_from_qa2sa = aa_to_exon(sa_from_qa, data['s_exinfo_list'])
        qex_from_sa2qa = aa_to_exon(qa_from_sa, data['q_exinfo_list'])


        ################
        # print out non-exon info

        print('\t'.join([str(data[x]) for x in fields[:end_field]]), end='')

        ################
        # edit the full text to insert the other aligned coordinates
        # (also re-order the regions query-first, then subject
        # for 'q_exalign_list', I need to add the subj_genome_coords sex_from_qa2sa
        #    and they need to be second
        # for 's_exalign_list', I need to add the query_genome_coords from qex_from_sa2qa
        #    and they need to be first

        q_exalign_out=[]
        for qx, q_exon in enumerate(data['q_exalign_list']):
            sg_start = sex_from_qa2sa[2*qx]
            sg_end = sex_from_qa2sa[2*qx+1]
            sg_replace="::%s:%d-%d}"%(sg_start['chrom'],sg_start['dpos'],sg_end['dpos'])

            this_outstr=re.sub(r'\}',sg_replace,q_exon.text)
            q_exalign_out.append(this_outstr)

        s_exalign_out=[]
        for sx, s_exon in enumerate(data['s_exalign_list']):
            qg_start = qex_from_sa2qa[2*sx]
            qg_end = qex_from_sa2qa[2*sx+1]
            qg_replace="{%s:%d-%d::"%(qg_start['chrom'],qg_start['dpos'],qg_end['dpos'])

            this_outstr=re.sub(r'\{',qg_replace,s_exon.text)
            s_exalign_out.append(this_outstr)

        print("\t|"+"|".join(q_exalign_out+s_exalign_out)+"\t"+line_data[-1])

################
# run the program ...

if __name__ == '__main__':
    main()

