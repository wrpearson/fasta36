#!/usr/bin/python

################################################################
# copyright (c) 2016 by William R. Pearson and The Rector &
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

import sys
import os
import argparse
import subprocess
import re

# python re-write of psisearch2_msa.pl
#
# logic:
# (1) do initial search
# (2) use results of initial search to produce MSA/PSSM for next search
# (3) do PSSM search
# (4) use results of PSSM search to produce MSA/PSSM for iterative step 3
#
################
#
# command:
# psisearch2_msa.py --query query_file --db database --num_iter N --evalue 0.002 --no_msa --int_mask none/query/random --end_mask none/query/random --tmp_dir results/ --domain --align --suffix --pgm ssearch/psiblast
#
################

################
# locations of required programs:
# (1) m89_btop_msa2.pl
# (2) ssearch
# (3) NCBI blast+ programs: psiblast/makeblastdb
# (4) NCBI datatool (required only for ssearch36 PSSMs)

pgm_bin = "/seqprg/bin"
pgm_data = "/seqprg/data"
ssearch_bin = pgm_bin+"/ssearch36"
psiblast_bin = pgm_bin+"/psiblast"
makeblastdb_bin = pgm_bin+"/makeblastdb"
datatool_bin = "%s/datatool -m %s/NCBI_all.asn" % (pgm_bin,pgm_data)
align2msa_lib = "m89_btop_msa2.pl"

annot_cmds = {'rpd3': '"!../scripts/ann_pfam28.pl --pfacc --db RPD3 --vdoms --split_over"',
              'rpd3nv':'"!../scripts/ann_pfam28.pl --pfacc --db RPD3 --split_over"',
              'pfam':'"!../scripts/ann_pfam30.pl --pfacc --vdoms --split_over"'}

num_iter = 5
evalue = 0.002
dom_flag = 0
align_flag = 0
int_mask = 'none'
end_mask = 'none'
srch_pgm = 'ssearch'
tmp_dir = ''
error_log = 0
rm_flag = 0
annot_type = ''
quiet = 0

################
# log_system()
# run system on string, logging first if error_log
#
def log_system (cmd, error_log):

  if (error_log) :
      sys.stderr.write(cmd+"\n")

  subprocess.call(cmd, shell=True)
#  print cmd

################
# sub get_ssearch_cmd()
# builds an ssearch command line with query, db, and pssm
#
def get_ssearch_cmd(query_file, db_file, pssm_file) :

  search_cmd = '%s -S -m 8CB -d 0 -E "1.0 0" -s BP62' % (ssearch_bin)

  if (annot_type) :
      search_cmd += " -V %s" % (annot_cmds[annot_type])

  if (pssm_file) :
      search_cmd += ' -P "%s 2"' % (pssm_file)

  search_cmd += " %s %s" % (query_file, db_file)

  return search_cmd


################
# sub get_psiblast_cmd()
# builds an ssearch command line with query, db, and pssm
#
def get_psiblast_cmd(query_file, db_file, pssm_file) :

  search_cmd = "%s -num_threads 4 -max_target_seqs 5000 -outfmt '7 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore score btop' -inclusion_ethresh %f -num_iterations 1 -db %s" % (psiblast_bin, evalue, db_file)

  if (pssm_file) :
      search_cmd += " -in_pssm %s" % (pssm_file)
  else :
      search_cmd += " -query %s" % (query_file)

  return search_cmd


################
# sub build_msa_pssm()
#
# given query, search output file (this_file_out), prev_boundary_file
# uses m89_btop_msa2.pl to generate PSSM in .asntxt or .asnbin format, also bound_file_out if align_flag
# (later - optionally deletes intermediate files)
#
# always produce a bound_file_out file to test for convergence
#
def build_msa_pssm(query_file, this_file_out,prev_bound_in, error_log) :

  (this_msa, this_hit_db, this_pssm_asntxt, this_pssm_asnbin, this_psibl_out, this_bound_out) =  (this_file_out+".msa",this_file_out+".hit_db",this_file_out+".asntxt",this_file_out+".asnbin",this_file_out+".psibl_out",this_file_out+".bnd_out")

  blastdb_err = this_file_out+".mkbldb_err"
  aln2msa_cmd = "%s --query %s --evalue %f --masked_lib_out=%s" % (align2msa_lib, query_file, evalue, this_hit_db)

  if (int_mask) :
      aln2msa_cmd += " --int_mask_type %s" % (int_mask)

  if (end_mask) :
      aln2msa_cmd += " --end_mask_type %s" % (end_mask)

  if (dom_flag) :
    aln2msa_cmd += " --domain"

  if (align_flag and prev_bound_in) :
      aln2msa_cmd += " --bound_file_in %s" %(prev_bound_in)

  # always produce this file to check for convergence
  aln2msa_cmd += " --bound_file_out %s" % (this_bound_out)

  log_system("%s %s > %s"%(aln2msa_cmd, this_file_out, this_msa), error_log)

  makeblastdb_cmd = "%s -in %s -dbtype prot -parse_seqids > %s" % (makeblastdb_bin, this_hit_db, blastdb_err)

  log_system(makeblastdb_cmd, error_log)

  buildpssm_cmd = "%s -max_target_seqs 5000 -outfmt 7 -inclusion_ethresh 100.0 -in_msa %s -db %s -out_pssm %s -num_iterations 1 -save_pssm_after_last_round" % (psiblast_bin, this_msa, this_hit_db, this_pssm_asntxt)

  log_system("%s > %s 2> %s.err" % (buildpssm_cmd, this_psibl_out, this_psibl_out), error_log)

  log_system("rm %s.p* %s" % (this_hit_db,blastdb_err), error_log)

  # remove uninformative error logs
  if (not error_log) :
    log_system("rm "+this_psibl_out+".err "+this_file_out+".err",error_log)

  if (srch_pgm != 'psiblast') :
      asn2asn_cmd = "%s -v %s -e %s" % (datatool_bin, this_pssm_asntxt, this_pssm_asnbin)
      log_system(asn2asn_cmd, error_log)
      return (this_pssm_asnbin, this_bound_out)
  else :
      return (this_pssm_asntxt, this_bound_out)

################
# sub has_converged()
# reads two boundary files and compares accessions
#
def has_converged(file1, file2) :

  f1_names = []
  f2_names = []

  with open(file1) as fd:
      for line in fd:
          line = line.rstrip('\n')
          fields = line.split('\t')
          f1_names.append(fields[0])

  with open(file2) as fd:
      for line in fd:
          line = line.rstrip('\n')
          fields = line.split('\t')
          f2_names.append(fields[0])

  # check for same length
  if (len(f1_names) != len(f2_names)) :
      return 0

  f1_names.sort()
  f2_names.sort()

  for i,v in enumerate(f1_names) :
    if (f2_names[i] != v) :
        return 0

  return 1

# main()

srch_subs = {'ssearch' : get_ssearch_cmd,
	     'psiblast': get_psiblast_cmd}

pgm_command =  "# "+" ".join(sys.argv);
if (error_log) :
    sys.stderr.write('pgm_command\n')

arg_parse = argparse.ArgumentParser(description='Iterative search with SSEARCH/PSIBLAST')
arg_parse.add_argument('--query', dest='query_file', action='store',help='query sequence file')
arg_parse.add_argument('--sequence', dest='query_file', action='store',help='query sequence file')
arg_parse.add_argument('--db', dest='db_file', action='store',help='sequence database name')
arg_parse.add_argument('--database', dest='db_file', action='store',help='sequence database name')
arg_parse.add_argument('--dir', dest='tmp_dir', action='store',help='directory for result and tmp_file output')
arg_parse.add_argument('--evalue', dest='evalue', default=0.002, type=float, action='store',help='E()-value threshold for inclusion in PSSM')
arg_parse.add_argument('--annot_db', dest='annot_type', action='store',help='source of domain annotations')
arg_parse.add_argument('--suffix', dest='suffix', action='store',help='suffix for result output')
arg_parse.add_argument('--out_name', dest='file_out', action='store',help='result file name')
arg_parse.add_argument('--iter', dest='num_iter', default=5, type=int, action='store',help='number of iterations')
arg_parse.add_argument('--in_pssm', dest='prev_pssm', action='store',help='initial PSSM')
arg_parse.add_argument('--in_bounds', dest='prev_bound_in', type=str, action='store',help='initial boundaries')
arg_parse.add_argument('--domain', dest='dom_flag', action='store_true',help='use domain annotations')
arg_parse.add_argument('--align', dest='align_flag', action='store_true',help='use alignment boundaries')
arg_parse.add_argument('--pgm', dest='srch_pgm', action='store',default='ssearch',help='search program: ssearch/psiblast')
arg_parse.add_argument('--query_seed', dest='query_mask', action='store_true',help='use query seeding')
arg_parse.add_argument('--int_seed', dest='int_mask', action='store',default='none',help='sequence masking: none/query/random')
arg_parse.add_argument('--end_seed', dest='end_mask', action='store',default='none',help='sequence masking: none/query/random')
arg_parse.add_argument('--int_mask', dest='int_mask', action='store',default='none',help='sequence masking: none/query/random')
arg_parse.add_argument('--end_mask', dest='end_mask', action='store',default='none',help='sequence masking: none/query/random')
arg_parse.add_argument('--save_list', dest='tmp_file_list', action='store',help='temporary extensions saved')
arg_parse.add_argument('--save_all', dest='save_all', action='store_true',help='save all temporary files')
arg_parse.add_argument('--delete_all', dest='delete_tmp', action='store_true',help='delete all temporary files')
arg_parse.add_argument('--delete_bnd', dest='delete_bnd', action='store_true',help='delete boundary temporary file')
arg_parse.add_argument('--quiet', dest='quiet', action='store_true',help='fewer messages')
arg_parse.add_argument('-Q', dest='quiet', action='store_true',help='fewer messages')

args = arg_parse.parse_args()
if (args.quiet) :
    quiet = args.quiet

if (args.srch_pgm) :
    srch_pgm = args.srch_pgm

if (not quiet) :
    print pgm_command

del_file_ext = ["msa","psibl_out","hit_db","asntxt","asnbin"]

if (re.match('psiblast',srch_pgm)) :
    del_file_ext.pop()


if (args.query_mask) :
    if (args.int_mask == 'none') :
        args.int_mask = 'query'
    if (args.end_mask == 'none') :
        args.end_mask = 'query'

delete_bnd = args.delete_bnd
if (args.delete_tmp) :
    delete_bnd = 1
elif (args.save_all) :
    delete_file_ext = ()
    args.tmp_file_list = ''
    delete_bnd = 0

save_file_ext = {}
if (args.tmp_file_list) :
    new_del_file_ext = []
    for ext in re.split(",\s*",args.tmp_file_list) :
        save_file_ext[ext] = 1
    for ext in del_file_ext :
        if (not (ext in save_file_ext)):
            new_del_file_ext.append(ext)
    del_file_ext = new_del_file_ext[:]

this_iter = "it1"

query_pref = query_file = args.query_file
m = re.search(r'([\w\.]+)$',str(args.query_file))
query_pref = m.groups()[0]

if (not args.file_out) :
    file_out = query_pref
else :
    file_out = args.file_out

this_file_out = this_file_pref = file_out+"."+this_iter
if (args.suffix) :
    this_file_out = this_file_pref+"."+args.suffix
if (args.tmp_dir) :
    this_file_out = args.tmp_dir+"/"+this_file_out

####
# parse output to build PSSM
# generate output filenames

prev_file_out = this_file_out

# do the first search
search_str = srch_subs[srch_pgm](args.query_file, args.db_file, args.prev_pssm)
log_system(search_str+" > "+this_file_out+" 2> "+this_file_out+".err", error_log)

prev_file_out = this_file_out

(this_pssm, this_bound_out) = build_msa_pssm(args.query_file, this_file_out, args.prev_bound_in, error_log)
# now have necessary files for next iteration

it=2
while (it <= args.num_iter) :

  prev_pssm = this_pssm
  prev_bound_in = this_bound_out
  ####
  # build filename for this iteration
  this_file_out = this_file_pref = "%s.it%d" % (file_out,it)
  if (args.suffix) :
      this_file_out = this_file_pref+"."+args.suffix
  if (args.tmp_dir) :
      this_file_out = args.tmp_dir+"/"+this_file_out

  search_str = srch_subs[srch_pgm](args.query_file, args.db_file, prev_pssm)
  log_system("%s > %s 2> %s" % (search_str,this_file_out,this_file_out+".err"), error_log)

  if (len(del_file_ext)):
      del_file_list = [ prev_file_out+'.'+ext for ext in del_file_ext]
      log_system('rm '+' '.join(del_file_list),error_log)

  prev_file_out = this_file_out

  (this_pssm, this_bound_out) = build_msa_pssm(query_file, this_file_out, prev_bound_in, error_log)

  if (has_converged(prev_bound_in, this_bound_out)) :
      if (not quiet) :
          sys.stderr.write(" %s %s %s %s converged (%d iterations)\n" % (sys.argv[0], srch_pgm, query_file, args.db_file, it))

      if (len(del_file_ext)):
          del_file_list = [ prev_file_out+'.'+ext for ext in del_file_ext]
          log_system('rm '+' '.join(del_file_list),error_log)

      if (delete_bnd) :
          log_system("rm "+prev_bound_in,error_log)

      exit(0)

  if (delete_bnd) :
      log_system("rm "+prev_bound_in,error_log)

  it += 1

if (len(del_file_ext)):
    del_file_list = [ prev_file_out+'.'+ext for ext in del_file_ext]
    log_system('rm '+' '.join(del_file_list),error_log)

if (delete_bnd):
    log_system("rm "+this_bound_out,error_log)

if (not quiet) :
    sys.stderr.write(" %s %s %s %s finished (%d iterations)\n" % (sys.argv[0], srch_pgm, query_file, args.db_file, it-1))

