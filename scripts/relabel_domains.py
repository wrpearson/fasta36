#!/usr/bin/env python3

# Given a blast_tabular file with search results from one or more
# protein queries, modify the domain color numbers (e.g. ~1, ~2) so
# that the query and subject domains use the same color numbers, even
# when the different annotation scripts may have assigned different
# numbers.
# 
################################################################
# copyright (c) 2018 by William R. Pearson and The Rector & Visitors
# of the University of Virginia */
# ###############################################################
# Licensed under the Apache License, Version 2.0 (the "License"); you
# may not use this file except in compliance with the License.  You
# may obtain a copy of the License at
# http://www.apache.org/licenses/LICENSE-2.0 Unless required by
# applicable law or agreed to in writing, software distributed under
# this License is distributed on an "AS IS" BASIS, WITHOUT WRRANTIES
# OR CONDITIONS OF ANY KIND, either express or implied.  See the
# License for the specific language governing permissions and
# limitations under the License.
# ###############################################################


import fileinput
import sys
import re
import argparse

from rename_exons import *

def replace_dom_number(line):

    out_str = ''
    if (not re.search(r'~',line)):
        return line

    (info, num, vdom) = re.search(r'^([^~]+)~(\d+)(v?)$',line).groups()
    if (vdom is None):
        vdom=''

    if (num in homolog_dict):
        return "%s~h%s%s" % (info, str(homolog_dict[num]['num']), vdom)
        
    else:
        name = line.split(" ")[-1].split("{")[0]
        if (name == "NODOM"):
            return line
        else:
            if (name in nonhomolog_dict):
                return '~'.join(line.split('~')[:-1]) + "~" + str(nonhomolog_dict[name])
    return out_str


################
# __main__ function
#

e_thresh = 1e-6
q_thresh = 30.0

homolog_dict = {}
nonhomolog_dict = {}

def main():

    # print "#"," ".join(sys.argv)

    hom_color = 1
    n_hom_color = 11

    parser=argparse.ArgumentParser(description='relabel_domains.py result_file.m8CB')

    parser.add_argument('--have_qslen', help='bl_tab fields include query/subject lengths',dest='have_qslen',action='store_true',default=False)
    parser.add_argument('--dom_info', help='raw domain coordinates included',action='store_true',default=False)
    parser.add_argument('files', metavar='FILE', nargs='*', help='files to read, if empty, stdin is used')

    args=parser.parse_args()

    end_field = -1
    data_fields_reset=False

    (fields, end_field) = set_data_fields(args, [])

    if (args.have_qslen and args.dom_info):
        data_fields_reset=True


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
        data = parse_protein(line_data,fields,'')	# get score/alignment/domain data

        if (len(data['sdom_list'])==0 and len(data['qdom_list'])==0):
            print(line)	# no domains to be edited, print stripped line and continue
            continue

        ################
        # have domains, check if significant and new, or old and known
        # goals are: (1) consistent coloring between query and subject for same domain
        #            (2) homologous domains get special labels
        # need dict of good domain names

        ################
        # check to update doms with good E()-value
        if float(data['evalue']) <= e_thresh:
            for q_dom in data['qdom_list']:
                if (float(q_dom.q_score) >= q_thresh and q_dom.name not in homolog_dict ):
                    homolog_dict['q_dom.name'] = q_dom_color
                    dom_color += 1

            for s_dom in data['sdom_list']:
                if (float(s_dom.q_score) >= q_thresh and s_dom.name not in homolog_dict):
                    homolog_dict['s_dom.name'] = s_dom.color
                    hom_color += 1
        else:
            for s_dom in data['sdom_list']:
                if (s_dom.name not in homolog_dict):
                    nonhomolog_dict['s_dom.name'] = s_dom.color
                    n_hom_color += 1


        ################
        # done storing good domains, write things out

        btab_str = '\t'.join(str(data[x]) for x in fields[:end_field])

        for s_dom in data['sdom_list']:
            if (s_dom.name in homolog_dict):
                s_dom.color=homolog_dict[s_dom.name]
            elif (s_dom.name in nonhomolog_dict):
                s_dom.color=nonhomolog_dict[s_dom.name]
        

        dom_bar_str = ''
        for dom in sorted(data['qdom_list']+data['sdom_list'],key=lambda r: r.idnum):
            dom_bar_str += dom.make_bar_str()

        print(btab_str+dom_bar_str)


if __name__ == '__main__':
    main()
