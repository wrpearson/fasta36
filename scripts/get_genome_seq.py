#!/usr/bin/env python3

################
## get_hg38_bed.py parses an HG38 coordinate into a pseudo-bed entry,
## and runs bedtools getfasta to return the fasta sequence
##

import sys
import re
from subprocess import Popen, PIPE, STDOUT
import shlex
import argparse

## a genome_loc should look like: chr#:start-stop
## if stop < start, coordinates are reversed

genome_dict={'hg38':'genome_dna/hg38/reference.fa',
             'mm10':'genome_dna/mm10/reference.fa',
             'rn6':'genome_dna/rn6/rn6.fa',
             'bosTau9':'genome_dna/bosTau9/bosTau9.fa'}

parser=argparse.ArgumentParser(description='get_genome_seq.py : get fasta sequence from genome coordinates ')
parser.add_argument('--genome', help='genome: hg38 | mm10 | rn6 | bosTau9',dest='genome',action='store',default='hg38')
parser.add_argument('coords', help='genome coordinates chr1:12345-54321', nargs='*')

args=parser.parse_args()

bed_cmd = 'bedtools getfasta -fi $RDLIB2/%s -bed stdin' % (genome_dict[args.genome])

bed_lines = ''
for genome_loc in args.coords:

    chrom, g_range = genome_loc.split(':')
    g_start, g_end = g_range.split('-')

    if (g_start > g_end):
        g_start, g_end = g_end, g_start

    g_start, g_end = int(g_start), int(g_end)
    g_start -= 1

    bed_lines += '%s\t%d\t%d\n' % (chrom, g_start, g_end)

bed_p = Popen(bed_cmd, stdout=PIPE, stdin=PIPE, stderr=STDOUT, shell=True, encoding='utf-8')
out, err = bed_p.communicate(input=bed_lines)

for line in out.split('\n'):
    if (line and line[0]=='>'):
        (chrom, start, stop) = re.search(r'>([^:]+):(\d+)\-(\d+)',line).groups()
        print(line + " @C:%s" % (start))
    elif (line):
        print(line)


