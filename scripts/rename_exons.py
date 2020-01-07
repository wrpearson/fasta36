#!/usr/bin/env python3
# 
# given a -m8CB file with exon annotations for the query and subject,
# adjust the subject exon names to match the query exon names
#
# see test_py.sh for sample use
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
import copy

################
# "domain" class that describes a domain/exon alignment annotation
#
class DomAlign:
    def __init__(self, name, info, color, qstart, qend, sstart, send, raw_score, bit_score, ident, qscore, RXRState, fulltext):
        self.name = name
        self.info = info
        self.color_type = ''
        if (not re.search(r'^\d+$',color)):
            m=re.search(r'^(\d+)([a-z]?\w*)$',color)
            if (m):
                (self.color, self.color_type) = m.groups()
                self.color = int(self.color)
        else:
            self.color = int(color)

        self.q_start = qstart
        self.q_end = qend
        self.s_start = sstart
        self.s_end = send
        self.raw_score = raw_score
        self.bit_score = bit_score
        self.percid = ident
        self.q_score = qscore
        self.rxr = RXRState
        self.idnum = 0
        self.overlap_list = []
        self.info_dom = None
        self.text = fulltext
        self.out_str = ''
        self.over_cnt = 0

    def append_overlap(self, overlap_dict):
        self.overlap_list.append(overlap_dict)

    def __str__(self):
#        return "[%d]name: %s : %i-%i : %i-%i I=%.1f Q=%.1f %s" % (self.idnum, self.name, self.q_start, self.q_end, self.s_start, self.s_end, self.percid, self.q_score, self.rxr)
        return "[%d:%s] %i-%i:%i-%i::%s [over:%d]" % (self.idnum, self.rxr, self.q_start, self.q_end, self.s_start, self.s_end, self.name,len(self.overlap_list))

    def print_bar_str(self):  # checking for 'NADA' 
        if (not self.out_str):
            self.out_str = self.text
        return str("|%s"%(self.out_str))

    def make_bar_str(self):  # create original string from values
        bar_str = "|%s:%d-%d:%d-%d:s=%d;b=%.1f;I=%.3f;Q=%.1f;C=%s%s~%d" % (
            self.rxr, self.q_start, self.q_end, self.s_start, self.s_end,
            self.raw_score, self.bit_score, self.percid, self.q_score, self.name, self.info, self.color)

        if (self.color_type):
            bar_str += self.color_type
        return bar_str

################
# "exonInfo" class describes raw (un-aligned) exons with genome coordinates
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


# Parses domain annotations after split at '|'
#|RX:1-38:3-40:s=37;b=17.0;I=0.289;Q=15.9;C=exon_1~1
#|RX:39-67:41-69:s=78;b=35.8;I=0.483;Q=68.7;C=exon_2~2
#|XR:1-67:3-69:s=115;b=52.8;I=0.373;Q=116.3;C=exon_1~1
#|RX:68-117:72-113:s=14;b=6.4;I=0.385;Q=0.0;C=exon_3~3
#|XR:68-124:70-119:s=-11;b=0.0;I=0.378;Q=0.0;C=exon_2~2
#|XR:125-167:120-165:s=39;b=17.9;I=0.429;Q=18.5;C=exon_3~3
#|RX:118-176:114-175:s=24;b=11.0;I=0.411;Q=1.5;C=exon_4~4
#|RX:177-200:176-198:s=27;b=12.4;I=0.435;Q=4.0;C=exon_5~5
#|XR:168-200:166-198:s=12;b=5.5;I=0.419;Q=0.0;C=exon_4~4
#
def parse_domain(text):
    # takes a domain in string form, turns it into a domain object
    # looks like: RX:5-82:5-82:s=397;b=163.1;I=1.000;Q=453.6;C=C.Thioredoxin~1
    # could also look like: RX:5-82:5-82:s=397;b=163.1;I=1.000;Q=453.6;C=C.Thioredoxin{PF012445}~1

    # get RX/XR and qstart/qstop sstart/sstop as strings
    m = re.search(r'^(\w+):(\d+)-(\d+):(\d+)-(\d+):',text)
    if (m):
        (RXRState, qstart_s, qend_s, sstart_s, send_s) = m.groups()
    else:
        sys.stderr.write("could not parse location: %s\n"%(text))

    # get score, bits, identity, Q info
    m = re.search(r's=(\-?\d+);b=(\-?[\d\.]+);I=([\d\.]+);Q=(\-?\d+\.\d*);',text)
    if (m):
        (r_score_s, b_score_s, ident_s, qscore_s) = m.groups()
    else:
        sys.stderr.write("Error: no scores: %s\n" %(text))
        r_score_s = b_score_s = qscore_s = "-1.0"

    # get domain name/color (and possibly {info})

    (name, color_s) = re.search(r';C=([^~]+)~(.+)$',text).groups()
    info_s=""

    if (re.search(r'\}$',name)):
        (name, info_s) = re.search(r'([^\{]+)(\{[^\}]+\})$',name).groups()

    dom_align = DomAlign(name, info_s, color_s, int(qstart_s), int(qend_s), int(sstart_s), int(send_s), 
                         int(r_score_s), float(b_score_s), float(ident_s),float(qscore_s), RXRState, text)

    return dom_align

# dom_info is like domain, but no scores
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

    gene_re = re.search(r'^\{([\w\.]+):(\d+)\-(\d+)\}',info)
    if (gene_re):
        (chrom, d_start, d_end) = gene_re.groups()
    else:
        (chrom, d_start, d_end) = ('',-1,-1)
#        sys.stderr.write("genome info not found: %s\n" % (text))

    q_target = True;
    if (RXRState == 'XD'):
        q_target = False

    exon_info = exonInfo(name, q_target, int(start_s), int(end_s), chrom, int(d_start), int(d_end), text)

    return exon_info

def overlap_fract(qdom, sdom):
    # checks if a query and subject domain overlap
    # if they do, return the amount of overlap with respect to each domain
    # how much of query is covered by subject, how much of subject is covered by query

    q_overlap = 0.0
    s_overlap = 0.0

    qq_len = qdom.q_end-qdom.q_start+1   # query alignment length in query coordinates
    qs_len = qdom.s_end-qdom.s_start+1   # query alignment length in subj coordinates
    sq_len = sdom.q_end-sdom.q_start+1   # subj alignment length in query coordinates
    ss_len = sdom.s_end-sdom.s_start+1   # subj alignment length in subject coordinates

    case = -1

    # case (0) no overlap at all
    if (qdom.q_end < sdom.q_start or sdom.s_end < qdom.s_start or qdom.q_start > sdom.q_end or sdom.q_start > qdom.q_end) :
        case = 0
        q_overlap = s_overlap = 0.0
    # case (1) query surrounds subject
    elif (qdom.q_start <= sdom.q_start and qdom.q_end >= sdom.q_end):
        case = 1
        s_overlap = 1.0
        q_overlap = float(sq_len)/qq_len
    # case (2) subject surrounds query
    elif (sdom.s_start <= qdom.s_start and sdom.s_end >= qdom.s_end):
        case = 2
        q_overlap = 1.0
        s_overlap = float(qs_len)/ss_len
    # case (3) query left of subject
    elif (qdom.q_start <= sdom.q_start and qdom.q_end <= sdom.q_end):
        case = 3
        q_overlap = float(qdom.q_end-sdom.q_start+1)/qq_len
        s_overlap = float(qdom.s_end-sdom.s_start+1)/ss_len
    # case (4) subject of left of query
    elif (sdom.s_start <= qdom.s_start and sdom.s_end <= qdom.s_end):
        case = 4
        q_overlap = float(sdom.q_end-qdom.q_start+1)/qq_len
        s_overlap = float(sdom.s_end-qdom.s_start+1)/ss_len

    if (q_overlap > 1.0 or s_overlap > 1.0):
        if (1):
            sys.stderr.write("***%i: qdom: %s sdom: %s\n"% (case,str(qdom),str(sdom)))
            sys.stderr.write(" ** qover %.3f sover: %.3f\n"% (q_overlap, s_overlap))
            sys.stderr.write(" ** qq_len: %d qs_len: %d ss_len: %d sq_len %d\n"%(qq_len, qs_len, ss_len, sq_len))

    return (q_overlap, s_overlap)

####
# parse_protein(result_line)
# takes a protein in string format, turns it into a dictionary properly
# looks like:   sp|P30711|GSTT1_HUMAN   up|Q2NL00|GSTT1_BOVIN   86.67   240     32      0       1       240     1       240     1.4e-123        444.0   16VI7DR6IT3IR15KQ3AI6TI11TA7YH8RC12TA3SN10FL10QETM2AT6VMTA2LV2DG4ND6PS24EK6TA11DV14FSPQ5IL3LMML1WK5RQ   |XR:4-76:4-76:s=327;b=134.6;I=0.895;Q=367.8;C=C.Thioredoxin~1|RX:5-82:5-82:s=356;b=146.5;I=0.902;Q=403.3;C=C.Thioredoxin~1|RX:83-93:83-93:s=52;b=21.4;I=0.818;Q=30.9;C=NODOM~0|XR:77-93:77-93:s=86;b=35.4;I=0.882;Q=72.6;C=NODOM~0|RX:94-110:94-110:s=88;b=36.2;I=0.882;Q=75.0;C=vC.GST_C~2v|XR:94-110:94-110:s=88;b=36.2;I=0.882;Q=75.0;C=vC.GST_C~2v|RX:111-201:111-201:s=409;b=168.3;I=0.868;Q=468.3;C=C.GST_C~2|XR:111-201:111-201:s=409;b=168.3;I=0.868;Q=468.3;C=C.GST_C~2|RX:202-240:202-240:s=154;b=63.4;I=0.795;Q=155.9;C=NODOM~0|XR:202-240:202-240:s=154;b=63.4;I=0.795;Q=155.9;C=NODOM~0
#
# returns [data[x] for x in fields] but also data['q/s_dom_list'] and data['q/sinfo_list']
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

    Qdom_list = []
    Sdom_list = []

    Qinfo_list = []
    Sinfo_list = []

    counter = 0

    if ('dom_annot' in data and len(data['dom_annot']) > 0):
        for dom_str in data['dom_annot'].split('|')[1:]:
            if (req_name and not re.search(req_name, dom_str)):
                continue

            counter += 1
            dom = parse_domain(dom_str)
            dom.idnum = counter
            if (dom.rxr == 'RX'):
                Qdom_list.append(dom)
            else:
                Sdom_list.append(dom)

    data['qdom_list'] = Qdom_list
    data['sdom_list'] = Sdom_list

    if ('dom_info' in data and len(data['dom_info']) > 0):
        for info_str in data['dom_info'].split('|')[1:]:
            if (req_name and not re.search(req_name, info_str)):
                continue
            if (not re.search(r'^[DX][XD]',info_str)):
                continue

            dinfo = parse_exon_info(info_str)

            if (dinfo.q_target):
                Qinfo_list.append(dinfo)
            else:
                Sinfo_list.append(dinfo)


        # put links to info_list into dom_list so info_list names can
        # be changed  -- give S/Qinfo's the S/Qdom ids of the overlapping domain

        find_info_overlaps(Qinfo_list, Qdom_list)
        find_info_overlaps(Sinfo_list, Sdom_list)
        
    data['qinfo_list'] = Qinfo_list
    data['sinfo_list'] = Sinfo_list

    return data

################
# "domain" : RX:1-38:3-40:s=37;b=17.0;I=0.289;Q=15.9;C=exon_1~1
# "name"   : like exon_2
# expanded for domain: RX:1-38:3-40:s=37;b=17.0;I=0.289;Q=15.9;C=exon_1{chr1:12345678-123456987}~1
#
def replace_name(domain_text, new_name, new_color_s):
    out = "=".join(domain_text.split("=")[:-1])  # out has everything to last '='

    old_name = domain_text.split(";C=")[-1]
    old_info=""
    
    if (re.search(r'\}~',old_name)):
        (old_info)=re.search(r'(\{[^\}]+\})~',old_name).group(1)

    if (not re.match(r'\d+',new_color_s)):
        new_color_s="0"
    out += "="+new_name+old_info+"~"+new_color_s                 # put it together
    return out

################
# check for overlaps using mid-point
#
def mid_overlaps(qdom_list, sdom_list):

    if (len(qdom_list) != len(sdom_list)):
        return False

    for ix, q_dom in enumerate(qdom_list):
        s_dom = sdom_list[ix]
        q_mid = q_dom.q_start + (q_dom.q_end - q_dom.q_start + 1)/2.0
        if not (q_mid >= s_dom.q_start and q_mid <= s_dom.q_end):
            return False

        q_qfract, q_sfract = overlap_fract(q_dom, s_dom)  # overlap from query perspective
        s_sfract, s_qfract = overlap_fract(s_dom, q_dom)  # overlap from subject perspective

        q_dom.overlap_list.append({"dom": s_dom, "q_over": q_qfract, "s_over": q_sfract})
        s_dom.overlap_list.append({"dom": q_dom, "q_over": s_qfract, "s_over": s_sfract})

    return True

################
# find_overlaps -- populates dom.overlap_list for qdoms, sdoms
#
def find_overlaps(qdom_list, sdom_list, over_thresh):
    # find qdom, sdom overlaps in O(N) time
    #

    if (len(sdom_list) == 0 or len(qdom_list)==0):
        return

    if (len(sdom_list) == len(qdom_list)):  # same number of domains
        if (mid_overlaps(qdom_list, sdom_list)):
            return;
        else:
            for d in qdom_list:
                d.overlap_list = []
            for d in sdom_list:
                d.overlap_list = []

        
    qdom_queue = [x for x in qdom_list]	# build a duplicate list
    sdom_queue = [x for x in sdom_list]

    qdom = qdom_queue.pop(0)	# get the first element of each
    sdom = sdom_queue.pop(0)

    while (True):
        pop_s = pop_q = False

        q_qfract, q_sfract = overlap_fract(qdom, sdom)  # overlap from query perspective
        if (q_qfract > over_thresh  or q_sfract > over_thresh):
            qdom.append_overlap({"dom": sdom, "q_over": q_qfract, "s_over": q_sfract})
            qdom.over_cnt += 1

        s_sfract, s_qfract = overlap_fract(sdom, qdom)  # overlap from query perspective
        if (s_qfract > over_thresh  or s_sfract > over_thresh):
            sdom.append_overlap({"dom": qdom, "q_over": s_qfract, "s_over": s_sfract})
            sdom.over_cnt += 1

        # check to see if we've used up the domain
        if (qdom.s_end >= sdom.s_end):
            pop_s = True
        # else there are more s_dom's that are part of this q_dom

        if (sdom.q_end >= qdom.q_end):
            pop_q = True
        # else there are more q_dom's that are part of this s_dom

#        print 'QS: %s %s\t%s %s' %(pop_q, pop_s, qdom, sdom)

        if (len(qdom_queue) > 0):
            if (pop_q):	# done with this qdom, get next
                qdom = qdom_queue.pop(0)
        elif (pop_q):  # don't break until we try to get the next domain
            break;

        if (len(sdom_queue) > 0):
            if (pop_s):	# done with this sdom, get next
                sdom = sdom_queue.pop(0)
        elif (pop_s):  # don't break until we try to get the next domain
            break;
    ####
    # all done with overlaps

    # # print "overlaps done"
    # for qd in qdom_list:
    #     print qd, qd.over_cnt
    #     for sd in qd.overlap_list:
    #         print " s: q_over %.3f s_over: %.3f %s" % (sd['q_over'], sd['s_over'], str(sd['dom']))
    # print "===="

    # for sd in sdom_list:
    #     print sd, sd.over_cnt
    #     for qd in sd.overlap_list:
    #         print " q: q_over %.3f s_over: %.3f %s" % (qd['q_over'], qd['s_over'], str(qd['dom']))
    # print "===="

################
# info_overlaps -- populates dom.overlap_list for qdoms, sdoms
#
def find_info_overlaps(info_list, dom_list):

    if (len(info_list) == 0 or len(dom_list)==0):
        return

    info_queue = [x for x in info_list]	# build a duplicate list
    dom_queue = [x for x in dom_list]

    info = info_queue.pop(0)	# get the first element of each
    dom = dom_queue.pop(0)

    while (True):
        pop_d = pop_i = False

        if (dom.rxr == 'RX'):  # use dom.q_start/q_end
            if (dom.q_end < info.p_start):
                pop_d = True
            elif (dom.q_end >= info.p_start and dom.q_start <= info.p_end):  # overlap
                dom.info_dom = info
                pop_d = True
                pop_i = True
            elif (info.p_end < dom.q_start):
                pop_i = True
            
        else:			# use dom.s_start/s_end
            if (dom.s_end < info.p_start):
                pop_d = True
            elif (dom.s_end >= info.p_start and dom.s_start <= info.p_end):  # overlap
                dom.info_dom = info
                pop_d = True
                pop_i = True
            elif (info.p_end < dom.s_start):
                pop_i = True

        if (len(info_queue) > 0):
            if (pop_i):	# done with this info, get next
                info = info_queue.pop(0)
        elif (pop_i):  # don't break until we try to get the next domain
            break;

        if (len(dom_queue) > 0):
            if (pop_d):	# done with this dom, get next
                dom = dom_queue.pop(0)
        elif (pop_d):
            break;

################
# build_multi_dict -- builds of dictionaries of multiple overlaps in
#                     qdom.overlap_list or sdom.overlap_list
# returns multi_dict[idnum]
#
def build_multi_dict(dom_list):
    # this code looks for xdom's that are associated with multiple ydoms
    #
    multi_dict = {}  # dict of {qids:/sdom:/qdoms:[]}
    for dom in dom_list:   # for each subject domain
        if (dom.over_cnt <= 1):
            continue

        multi_id_list = []
        multi_dom_list = []
        multi_q_cnt = 0
        for xd_over_yd in dom.overlap_list:  	      # a set of q_doms that overlap the subject
            multi_q_cnt += 1
            multi_id_list.append(xd_over_yd["dom"].idnum)  # these are q_dom idnum's
            multi_dom_list.append(xd_over_yd["dom"])          # these are q_doms

        if (multi_q_cnt > 1):    # only save when two (or more) overlaps
            multi_dict[dom.idnum] = {"yids": multi_id_list, "ydoms":multi_dom_list, 'xdom':dom}

    # # print out current multi_q_list
    # print "--- multi_q dict ---"
    # for db in multi_dict.keys():
    #     print "sdom: %s"%(db)
    #     for ix, qd in enumerate(multi_dict[db]['ydoms']):
    #         print " %d %d: %s"%(ix, multi_dict[db]['yids'][ix], qd)

    # print "--- multi_dict done"

    return multi_dict

################
# find_best_id() -- returns id of domain with longest 'q_over'
#
def find_best_id(overlap_list, over_type):

    max_fract = 0.0
    max_idnum = 0
    for over_d in overlap_list:
        if (over_d[over_type] > max_fract):
            max_idnum = over_d['dom'].idnum
            max_fract = over_d[over_type]

    return max_idnum

################################################################
# final labeling routine -- leave qdom's alone, modify sdoms based on qdoms.
################
# sdom's in more than one qdom are in multi_q_dict[]
# qdom's in more than one sdom are in multi_s_dict[]
# everyone else just gets the qdom name
# returns sdom_displayed_dict{idnum} -- the set of sdoms that have been modified
#
# 13-Nov-2018 -- ensure that there is an info_dom before replacing info_dom.text
#
def label_doms(qdom_list, sdom_list, multi_q_dict, multi_s_dict):

    sdom_displayed_dict = {}
    for qdom in qdom_list:
        # qdom's stay the same
        qdom.out_str = qdom.text

        # check for s_doms with multiple q_doms
        if (qdom.idnum in multi_s_dict):
            # find the best, name it exon_X, find the rest, name it qdom.name
            multi_s_entry = multi_s_dict[qdom.idnum]
            best_id = find_best_id(qdom.overlap_list,'q_over')  # find sdom with most overlap
            for s_over in qdom.overlap_list:   # find the sdom's that overlap this qdom
                sdom = s_over['dom']
                if (sdom.idnum == best_id):
                    sdom.out_str = replace_name(sdom.text, qdom.name, str(qdom.color))
                    if (sdom.info_dom):
                        sdom.info_dom.out_str = replace_name(sdom.info_dom.text,qdom.name, str(qdom.color))
                else:
                    sdom.out_str = replace_name(sdom.text, "exon_X","0")
                    if (sdom.info_dom):
                        sdom.info_dom.out_str = replace_name(sdom.info_dom.text,"exon_X","0")
                sdom_displayed_dict[sdom.idnum] = sdom;
            continue	# prevents re-labeling later

        # check for q_doms with multiple doms
        for sd_over in qdom.overlap_list:
            sdom = sd_over['dom']
            # it might make sense to do this in a second for loop after
            # all the multiple stuff is done
            if (sdom.idnum not in multi_q_dict):
                # this is the simplest case -- sdom.text gets qdom.name
                if (sdom.idnum not in sdom_displayed_dict):
                    sdom.out_str = replace_name(sdom.text, qdom.name, str(qdom.color))
                    if (sdom.info_dom):
                        sdom.info_dom.out_str = replace_name(sdom.info_dom.text,qdom.name, str(qdom.color))
            else:
                # this sdom belongs to multiple q_doms, add each of those q_doms to the name
                exon_str='exon_'
                # "ydoms" here are the qdoms overlapped by sdom
                exon_str += '/'.join([ x.name.split("_")[1] for x in multi_q_dict[sdom.idnum]['ydoms']])
                sdom.out_str = replace_name(sdom.text, exon_str,"0")
                if (sdom.info_dom):
                    sdom.info_dom.out_str = replace_name(sdom.info_dom.text,exon_str,"0")

            sdom_displayed_dict[sdom.idnum]=sdom

    # done with labeling sdoms based on qdoms, but some may be unlabeled
    # check for missing s_doms
    while (len(list(sdom_displayed_dict.keys())) < len(sdom_list)):
        for sdom in sdom_list:
            if (sdom.idnum not in sdom_displayed_dict):
                sdom.out_str = replace_name(sdom.text, "exon_X","0")
                if (sdom.info_dom):
                    sdom.info_dom.out_str = replace_name(sdom.info_dom.text,"exon_X","0")

                sdom_displayed_dict[sdom.idnum] = sdom

    return sdom_displayed_dict

################
#
# aa_to_exon()  --- given a coordinate and the corresponding exon map, return the exon coordinate
# (can only be done for aligned exons)
#
# this version of the function must use an info_list, not an
# align_list, because it uses p_start/p_end rather than q_start/s_start, etc.
# a version using qp_start/sp_start would also need a target argument
#
def aa_to_exon(aa_coords, exon_info_list):

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
# 
def set_data_fields(args, line_data) :

    field_str = 'qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore BTOP dom_annot'
    field_qs_str = 'qseqid qlen sseqid slen pident length mismatch gapopen qstart qend sstart send evalue bitscore BTOP dom_annot'

    if (len(line_data) > 1) :
        if ((not args.have_qslen) and  re.search(r'\d+',line_data[1])):
            args.have_qslen=True

        if ((not args.dom_info) and re.search(r'^\|[DX][XD]\:',line_data[-1])):
            args.dom_info = True

    end_field = -1
    fields = field_str.split(' ')

    if (args.have_qslen):
        fields = field_qs_str.split(' ')

    if (args.dom_info):
        fields.append('dom_info')
        end_field = -2

    return (fields, end_field)

################################################################
#
# main program 
# print "#"," ".join(sys.argv)

def main():

    parser=argparse.ArgumentParser(description='rename_exons.py result_file.m8CB : re-label subject exons to match query')
    parser.add_argument('--have_qslen', help='bl_tab fields include query/subject lengths',dest='have_qslen',action='store_true',default=False)
    parser.add_argument('--dom_info', help='raw domain coordinates included',action='store_true',default=False)
    parser.add_argument('--fill_gcoords', help='fill in genomic coordinates',action='store_true',default=False)
    parser.add_argument('files', metavar='FILE', nargs='*', help='files to read, if empty, stdin is used')

    args=parser.parse_args()

    end_field = -1
    data_fields_reset=False

    (fields, end_field) = set_data_fields(args, [])

    if (args.have_qslen and args.dom_info):
        data_fields_reset=True

    saved_qdom_list = []
    qdom_list = []

    for line in fileinput.input(args.files):
    # pass through comments
        if (line[0] == '#'):
            print(line, end='')	# ',' because have not stripped
            continue

        ################
        # break up tab fields, check for extra fields
        line = line.strip('\n')
        line_data = line.split('\t')
        if (not data_fields_reset):     # look for --have_qslen number, --dom_info data, even if not set
            (fields, end_field) = set_data_fields(args, line_data)
            data_fields_reset = True

        ################
        # get exon annotations
        data = parse_protein(line_data,fields,"exon")	# get score/alignment/domain data

        if (len(data['sdom_list'])==0 and len(data['qdom_list'])==0):
            print(line)	# no domains to be edited, print stripped line and contine
            continue

        # qdom_list=[] outside of loop for cases where the qseqid==sseqid match is not first
        if len(data['qdom_list'])== 0:
            if data['qseqid'] == data['sseqid']:
                saved_qdom_list = [ copy.deepcopy(x) for x in data['sdom_list']]
                max_sdom_id=len(data['sdom_list'])+1
                for qdom in saved_qdom_list:
                    qdom.rxr = 'RX'
                    qdom.idnum = max_sdom_id
                    max_sdom_id += 1

            qdom_list = [copy.deepcopy(x) for x in saved_qdom_list]
        else:
            qdom_list = data['qdom_list']

        # print out non-exon info
    
        if (len(qdom_list) == 0):
            print(line)
            continue

        btab_str = '\t'.join(str(data[x]) for x in fields[:end_field])
        # print  # comment out for single line

        ################
        # find overlaps and multi-overlaps
        #
        find_overlaps(qdom_list,data['sdom_list'], 0.2)

        multi_q_dict = build_multi_dict(data['sdom_list'])  # keys are sdoms hitting multiple qdoms
        multi_s_dict = build_multi_dict(qdom_list)  # keys are qdoms hitting mulitple sdoms

        ################
        # label qdoms, relabel sdoms
        #
        sdom_displayed_dict = label_doms(qdom_list, data['sdom_list'], multi_q_dict, multi_s_dict)

        ################
        # print exon annotations
        #
        q_exon_list = data['qdom_list']

        s_exon_list = [sdom_displayed_dict[x] for x in list(sdom_displayed_dict.keys())]

        ################
        # if args.fill_gcoords, then do the transformations on the current exon lists
        
        if (args.fill_gcoords):
            sa_from_qa = []
            for q_ex in q_exon_list:
                sa_from_qa.append(q_ex.q_start)
                sa_from_qa.append(q_ex.q_end)

            # have list of coordinates, map them to exon
            sex_from_qa2sa = aa_to_exon(sa_from_qa,data['sinfo_list'])

            for iqx, q_ex in enumerate(q_exon_list):
                sg_start = sex_from_qa2sa[2*iqx]
                sg_end = sex_from_qa2sa[2*iqx+1]
                sg_replace="::%s:%d-%d}"%(sg_start['chrom'],sg_start['dpos'],sg_end['dpos'])
                q_ex.text=re.sub(r'\}',sg_replace,q_ex.text)
                q_ex.out_str=re.sub(r'\}',sg_replace,q_ex.out_str)
                
            qa_from_sa = []
            for s_ex in s_exon_list:
                qa_from_sa.append(s_ex.q_start)
                qa_from_sa.append(s_ex.q_end)

            # have list of coordinates, map them to exon
            qex_from_sa2qa = aa_to_exon(qa_from_sa,data['qinfo_list'])

            for isx, s_ex in enumerate(s_exon_list):
                qg_start = sex_from_qa2sa[2*iqx]
                qg_end = sex_from_qa2sa[2*iqx+1]
                qg_replace="{%s:%d-%d::"%(sg_start['chrom'],sg_start['dpos'],sg_end['dpos'])
                s_ex.text=re.sub(r'\{',qg_replace,s_ex.text)
                s_ex.out_str = re.sub(r'\{',qg_replace,s_ex.out_str)

        sorted_exon_list = sorted(q_exon_list+s_exon_list,key = lambda r: r.idnum)

        dom_bar_str = ''
        for exon in sorted_exon_list:
            # print exon.print_bar_str()     # for multi-line output
            dom_bar_str += exon.print_bar_str() 

        info_bar_str = ''
        for info in data['qinfo_list'] + data['sinfo_list']:
            info_bar_str += info.text 

        print('\t'.join((btab_str, dom_bar_str, info_bar_str)))

################
# run the program ...

if __name__ == '__main__':
    main()

