#!/usr/bin/env python

################################################################
# copyright (c) 2014,2015 by William R. Pearson and The Rector &
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
# clustal2fasta.pl 
################################################################
# clustal2fasta.pl takes a standard clustal format alignment file
# and produces the corresponding FASTA file.
#
# if --end_mask or --int_mask are set, then end or internal '-'s are converted to the query (first) sequence
# if --trim is set, then alignments beyond the beginning/end of the query sequence are trimmed
#
################################################################

import argparse
import fileinput
import re

################
#
# python re-write of clustal2fasta.pl
#
# in the future, modify for various query seeding strategies
################    

arg_parse = argparse.ArgumentParser(description='Convert clustal MSA to FASTA library')
arg_parse.add_argument('--query|--query_file', dest='query_file', action='store',help='query sequence file')
arg_parse.add_argument('files', metavar='FILE', nargs='*', help='files to read, if empty, stdin is used')
args=arg_parse.parse_args()

msa = {}
seq_ids = []

is_line1 = True
for line in fileinput.input(args.files):
    if is_line1:
        is_line1 = False
        continue
    line = line.strip()
    if not line:
        continue
    if re.search(r'^[\s:\*\+\.]+$',line):
        continue
    
    (seq_id, align) = re.split(r'\s+',line)

    if seq_id in msa:
        msa[seq_id] += align
    else:
        msa[seq_id] = align
        seq_ids.append(seq_id)

for seq_id in seq_ids:
    fmt_seq = re.sub(r'(.{0,60})',r'\1\n',msa[seq_id])
    print ">%s\n%s" % (seq_id, fmt_seq)
