#!/usr/bin/env python
# 
# given a -m8CB file with exon annotations for the query and subject,
# adjust the subject exon names to match the query exon names

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
# main "domain" class that describes a domain/exon alignment annotation
#
class Domain:
    def __init__(self, name, qstart, qend, sstart, send, ident, qscore, RXRState, fulltext):
        self.name = name
        self.q_start = qstart
        self.q_end = qend
        self.s_start = sstart
        self.s_end = send
        self.percid = ident
        self.qscore = qscore
        self.rxr = RXRState
        self.idnum = 0
        self.overlap_list = []
        self.bestfit = False
        self.text = fulltext
        self.out_str = ''
        self.over_cnt = 0
    def append_overlap(self, overlap_dict):
        self.overlap_list.append(overlap_dict)
#    def makebest(self):
#        self.bestfit = True
    def __str__(self):
#        return "[%d]name: %s : %i-%i : %i-%i I=%.1f Q=%.1f %s" % (self.idnum, self.name, self.q_start, self.q_end, self.s_start, self.s_end, self.percid, self.qscore, self.rxr)
        return "[%d:%s] %i-%i:%i-%i::%s" % (self.idnum, self.rxr, self.q_start, self.q_end, self.s_start, self.s_end, self.name)
    def print_bar_str(self):  # checking for 'NADA' 
        if (not self.out_str):
            self.out_str = self.text
        return str("|%s"%(self.out_str))

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

    (RXRState, qstart_s, qend_s, sstart_s, send_s) = re.search(r'^(\w+):(\d+)-(\d+):(\d+)-(\d+)',text).groups()

    m = re.search(r';Q=(\-?\d+\.\d*);',text)
    if (m):
        qscore_s = m.group(1)
    else:
        sys.stderr.write("Error: no Q= score: %s\n" %(text))
        qscore_s = -1.0
    v = re.search(r';I=(\d+\.\d*);',text)
    if (v):
        iscore_s = v.group(1)

    name = re.search(r';C=([\w\.]+)',text).group(1)
    Dom = Domain(name, int(qstart_s), int(qend_s), int(sstart_s), int(send_s), float(iscore_s),float(qscore_s), RXRState, text)

    return Dom

def overlap_fract(qdom, sdom):
    # checks if a query and subject domain overlap
    # if they do, return the amount of overlap with respect to each domain
    # how much of query is covered by subject, how much of subject is covered by query

    q_overlap = 0.0
    s_overlap = 0.0

    qq_len = qdom.q_end-qdom.q_start+1
    qs_len = qdom.s_end-qdom.s_start+1
    sq_len = sdom.q_end-sdom.q_start+1
    ss_len = sdom.s_end-sdom.s_start+1

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
    elif (sdom.q_start <= qdom.q_start and sdom.q_end >= qdom.q_end):
        case = 2
        q_overlap = 1.0
        s_overlap = float(qs_len)/ss_len
        # case (3) query left of subject
    elif (qdom.q_start <= sdom.q_start and sdom.q_end <= qdom.q_end):
        case = 3
        q_overlap = float(qdom.q_end-sdom.q_start+1)/qq_len
        s_overlap = float(sdom.s_end-qdom.s_start+1)/ss_len
        # case (4) query right of subject
    elif (sdom.q_start <= qdom.q_start and qdom.q_end <= sdom.q_end):
        case = 4
        q_overlap = float(sdom.q_end-qdom.q_start+1)/qs_len
        s_overlap = float(qdom.s_end-sdom.s_start+1)/sq_len

    if (q_overlap > 1.0 or s_overlap > 1.0):
        if (0):
            sys.stderr.write("***%i: qdom: %s sdom: %s\n"% (case,str(qdom),str(sdom)))
            sys.stderr.write(" ** qover %.3f sover: %.3f\n"% (q_overlap, s_overlap))
            sys.stderr.write(" ** qq_len: %d qs_len: %d ss_len: %d sq_len %d\n"%(qq_len, qs_len, ss_len, sq_len))
    return (q_overlap, s_overlap)

####
# parse_protein(result_line)
# takes a protein in string format, turns it into a dictionary properly
# looks like:   sp|P30711|GSTT1_HUMAN   up|Q2NL00|GSTT1_BOVIN   86.67   240     32      0       1       240     1       240     1.4e-123        444.0   16VI7DR6IT3IR15KQ3AI6TI11TA7YH8RC12TA3SN10FL10QETM2AT6VMTA2LV2DG4ND6PS24EK6TA11DV14FSPQ5IL3LMML1WK5RQ   |XR:4-76:4-76:s=327;b=134.6;I=0.895;Q=367.8;C=C.Thioredoxin~1|RX:5-82:5-82:s=356;b=146.5;I=0.902;Q=403.3;C=C.Thioredoxin~1|RX:83-93:83-93:s=52;b=21.4;I=0.818;Q=30.9;C=NODOM~0|XR:77-93:77-93:s=86;b=35.4;I=0.882;Q=72.6;C=NODOM~0|RX:94-110:94-110:s=88;b=36.2;I=0.882;Q=75.0;C=vC.GST_C~2v|XR:94-110:94-110:s=88;b=36.2;I=0.882;Q=75.0;C=vC.GST_C~2v|RX:111-201:111-201:s=409;b=168.3;I=0.868;Q=468.3;C=C.GST_C~2|XR:111-201:111-201:s=409;b=168.3;I=0.868;Q=468.3;C=C.GST_C~2|RX:202-240:202-240:s=154;b=63.4;I=0.795;Q=155.9;C=NODOM~0|XR:202-240:202-240:s=154;b=63.4;I=0.795;Q=155.9;C=NODOM~0
#
def parse_protein(line,fields):
    # last part (domain annotions) split('|') and parsed by parse_domain()

    data = {}
    data = dict(zip(fields, line.split('\t')))
    data['qseq_acc'] = data['qseqid'].split('|')[1]
    data['sseq_acc'] = data['sseqid'].split('|')[1]

    Qdom_list = []
    Sdom_list = []
    counter = 0

    if ('dom_annot' in data and len(data['dom_annot']) > 0):
        for dom_str in data['dom_annot'].split('|')[1:]:
            counter += 1
            dom = parse_domain(dom_str)
            dom.idnum = counter
            if (dom.rxr == 'RX'):
                Qdom_list.append(dom)
            else:
                Sdom_list.append(dom)

    data['qdom_list'] = Qdom_list
    data['sdom_list'] = Sdom_list
    return data

# "domain" looks like RX:1-38:3-40:s=37;b=17.0;I=0.289;Q=15.9;C=exon_1~1
# "name" looks like exon_2
def replace_name(domain, name):
    out = "=".join(domain.split("=")[:-1])
    ext = name.split("_")[-1]
    if (not re.match(r'\d+',ext)):
        ext="0"
    out += "="+name+"~"+ext
    return out

################
# find_overlaps -- populates dom.overlap_list for qdoms, sdoms
#
def find_overlaps(qdom_list, sdom_list):
    # find qdom, sdom overlaps in O(N) time
    #

    if (len(sdom_list) == 0 or len(qdom_list)==0):
        return

    qdom_queue = [x for x in qdom_list]	# build a duplicate list
    sdom_queue = [x for x in sdom_list]

    qdom = qdom_queue.pop(0)	# get the first element of each
    sdom = sdom_queue.pop(0)

    while (True):
        pop_s = pop_q = False

        q_qfract, q_sfract = overlap_fract(qdom, sdom)  # overlap from query perspective
        if (q_qfract > 0.0  or q_sfract > 0.0):
            qdom.append_overlap({"dom": sdom, "q_over": q_qfract, "s_over": q_sfract})
            qdom.over_cnt += 1

        s_sfract, s_qfract = overlap_fract(sdom, qdom)  # overlap from query perspective
        if (s_qfract > 0.0  or s_sfract > 0.0):
            sdom.append_overlap({"dom": qdom, "q_over": s_qfract, "s_over": s_sfract})
            sdom.over_cnt += 1

        # check to see if we've used up the domain
        if (qdom.s_end >= sdom.s_end):
            pop_s = True
        # else there are more s_dom's that are part of this q_dom

        if (sdom.s_end >= qdom.s_end):
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
                   sdom.out_str = replace_name(sdom.text, "exon_X")
                else:
                    sdom.out_str = replace_name(sdom.text, qdom.name)
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
                    sdom.out_str = replace_name(sdom.text, qdom.name)
            else:
                # this sdom belongs to multiple q_doms, add each of those q_doms to the name
                exon_str='exon_'
                # "ydoms" here are the qdoms overlapped by sdom
                exon_str += '/'.join([ x.name.split("_")[1] for x in multi_q_dict[sdom.idnum]['ydoms']])
                sdom.out_str = replace_name(sdom.text, exon_str)
            sdom_displayed_dict[sdom.idnum]=sdom

    # done with labeling sdoms based on qdoms, but some may be unlabeled
    # check for missing s_doms
    while (len(sdom_displayed_dict.keys()) < len(sdom_list)):
        for sdom in sdom_list:
            if (sdom.idnum not in sdom_displayed_dict):
                sdom.out_str = replace_name(sdom.text, "exon_X")
                sdom_displayed_dict[sdom.idnum] = sdom

    return sdom_displayed_dict

################################################################
#
# main program 
# print "#"," ".join(sys.argv)

parser=argparse.ArgumentParser(description='scan_exons.py result_file.m8CB : re-label subject exons to match query')
parser.add_argument('--have_qslen', help='bl_tab fields include query/subject lengths',dest='have_qslen',action='store_true',default=False)
parser.add_argument('files', metavar='FILE', nargs='*', help='files to read, if empty, stdin is used')
args=parser.parse_args()

field_str = 'qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore BTOP dom_annot'
field_qs_str = 'qseqid qlen sseqid slen pident length mismatch gapopen qstart qend sstart send evalue bitscore BTOP dom_annot'
fields = field_str.split(' ')
if (args.have_qslen):
    fields = field_qs_str.split(' ')

saved_qdom_list = []

for line in fileinput.input(args.files):
    # pass through commnents
    if (line[0] == '#'):
        print line,	# ',' because have not stripped
        continue

    ################
    # break up tab fields, get exon annotations
    data = parse_protein(line.strip('\n'),fields)	# get score/alignment/domain data

    if len(data['qdom_list'])== 0:
        if data['qseqid'] == data['sseqid']:
            saved_qdom_list = [ copy.deepcopy(x) for x in data['sdom_list']]
            max_sdom_id=len(data['sdom_list'])+1
            for qdom in saved_qdom_list:
                qdom.rxr = 'VX'
                qdom.idnum = max_sdom_id
                max_sdom_id += 1

        qdom_list = [copy.deepcopy(x) for x in saved_qdom_list]
    else:
        qdom_list = data['qdom_list']

    # print out non-exon info
    btab_str = '\t'.join(str(data[x]) for x in fields[:-1])
    # print  # comment out for single line

    ################
    # find overlaps and multi-overlaps
    #
    find_overlaps(qdom_list,data['sdom_list'])
    multi_q_dict = build_multi_dict(data['sdom_list'])  # keys are sdoms hitting multiple qdoms
    multi_s_dict = build_multi_dict(qdom_list)  # keys are qdoms hitting mulitple sdoms

    ################
    # label qdoms, relabel sdoms
    #
    sdom_displayed_dict = label_doms(qdom_list, data['sdom_list'], multi_q_dict, multi_s_dict)

    ################
    # print exon annotations
    #
    exon_list = data['qdom_list']+[sdom_displayed_dict[x] for x in sdom_displayed_dict.keys()]
    sorted_exon_list = sorted(exon_list,key = lambda r: r.idnum)
    bar_str = ''
    for exon in sorted_exon_list:
        # print exon.print_bar_str()     # for multi-line output
        bar_str += exon.print_bar_str() 
    print btab_str+'\t'+ bar_str
