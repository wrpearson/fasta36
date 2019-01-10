#!/usr/bin/env python
# 
# given a -m8CB file with exon annotations for the query and subject,
# provide a function that maps subject coordinates to query, or vice versa

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
import copy

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
    data = dict(zip(fields, line_data))
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

################################################################
# map_subj_aa(q_exons, s_exons, s_aa_pos)
#
# return s_dna_pos, q_aa_pos, q_dna_pos for given s_aa_pos coordinate
# 
# def map_subj_aa(q_exons, s_exons, s_aa_pos):
    
# go through the q_exons to find the appropriate segment
# go through the s_exons to find the appropriate segment


################
#

def set_data_fields(args, line_data) :

    field_str = 'qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore BTOP align_annot'
    field_qs_str = 'qseqid qlen sseqid slen pident length mismatch gapopen qstart qend sstart send evalue bitscore BTOP align_annot'

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
# print "#"," ".join(sys.argv)

def main():

    data_fields_reset=False

    parser=argparse.ArgumentParser(description='map_exon_coords.py result_file.m8CB saa:coord : map subject coordinate to query genomic coordinate')
    parser.add_argument('--have_qslen', help='bl_tab fields include query/subject lengths',dest='have_qslen',action='store_true',default=False)
    parser.add_argument('--exon_info', help='raw domain coordinates included',action='store_true',default=True)
    parser.add_argument('--subj_aa',help='subject aa coordinate to map',action='store',type=int,dest='subj_aa_coord',default=1)
    parser.add_argument('files', metavar='FILE', nargs='*', help='files to read, if empty, stdin is used')
    args=parser.parse_args()

    if (args.have_qslen and args.exon_info):
        data_fields_reset=True

    saved_qexon_list = []
    qexon_list = []

    for line in fileinput.input(args.files):
    # pass through comments
        if (line[0] == '#'):
            print line,	# ',' because have not stripped
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
        data = parse_protein(line_data,fields,"exon")	# get score/alignment/domain data

        # look at results:

        for ea in data['q_exalign_list']: print ea
        for ea in data['s_exalign_list']: print ea
        for ea in data['q_exinfo_list']: print ea
        for ea in data['s_exinfo_list']: print ea


################
# run the program ...

if __name__ == '__main__':
    main()

