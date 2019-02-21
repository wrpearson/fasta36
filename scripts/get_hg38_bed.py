#!/usr/bin/python

################
## get_hg38_bed.py parses an HG38 coordinate into a pseudo-bed entry,
## and runs bedtools getfasta to return the fasta sequence
##

import sys
import re
from subprocess import Popen, PIPE, STDOUT
import shlex

## a genome_loc should look like: chr#:start-stop
## if stop < start, coordinates are reversed


hg38_ref = 'genome_dna/hg38/reference.fa'
bed_cmd = 'bedtools getfasta -fi $RDLIB2/%s -bed stdin' % (hg38_ref)

bed_lines = ''
for genome_loc in sys.argv[1:]:

    chrom, g_range = genome_loc.split(':')
    g_start, g_end = g_range.split('-')

    if (g_start > g_end):
        g_start, g_end = g_end, g_start

    g_start, g_end = int(g_start), int(g_end)
    g_start -= 1

    bed_lines += '%s\t%d\t%d\n' % (chrom, g_start, g_end)

bed_p = Popen(bed_cmd, stdout=PIPE, stdin=PIPE, stderr=STDOUT, shell=True)
out, err = bed_p.communicate(input=bed_lines)

for line in out.split('\n'):
    if (line and line[0]=='>'):
        (chrom, start, stop) = re.search(r'>([^:]+):(\d+)\-(\d+)',line).groups()
        print line + " @C:%s" % (start)
    else:
        print line


