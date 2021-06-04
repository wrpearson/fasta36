#!/usr/bin/env perl

################################################################
# copyright (c) 2016,2020 by William R. Pearson and The Rector &
# Visitors of the University of Virginia
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
# 5-July-2020
# modification to allow older blast version, (2.6.0+), as current version (2.10.1+)
# has problem with asnbin/asntxt files
#

use warnings;
use strict;
use Getopt::Long;
use Pod::Usage;

################
# implementation of simple shell script to do iterative searches
#
# logic:
# (1) do initial search or take results from previous search with --prev_m89res
# (2) use results of initial search to produce MSA/PSSM for next search
# (3) do PSSM search
# (4) repeat at step (2)
#
################
#
# command:
# psisearch2_msa.pl --query query.file --in_msa msa.file --db database.file --num_iter N --pssm_evalue 0.002 --int_mask none/query/random --end_mask none/query/random --tmp_dir results/ --domain --align --out_suffix none --pgm ssearch/psiblast --prev_m89res prev_results.itx.m8CB.file --sel_res selected_accs.file --prev_bounds boundary.file
#
################

use vars qw( $query_file $db_file $num_iter $pssm_evalue $srch_evalue $int_mask $end_mask $query_mask $tmp_dir $dom_flag $align_flag $suffix $srch_pgm $file_out $help $shelp $error_log $rm_flag $annot_type $quiet);
use vars qw( $prev_m89res $m_format $prev_sel_res $this_iter $prev_msa $next_msa $prev_hitdb $next_hitdb $prev_pssm $next_pssm $prev_bound $next_bound_out $tmp_file_list $save_all $delete_bnd $delete_tmp $use_stdout);

################
# locations of required programs:
# (1) m89_btop_msa2.pl
# (2) ssearch
# (3) NCBI blast+ programs: psiblast/makeblastdb
# (4) NCBI datatool (required only for ssearch36 PSSMs)

my $pgm_bin = "/seqprg/bin";

my $pgm_data = "/seqprg/data";
my $ssearch_bin = "$pgm_bin/ssearch36";
my $align2msa_lib = "$pgm_bin/m89_btop_msa2.pl";
my $clustal2fasta = "$pgm_bin/clustal2fasta.pl";

## this has been added (5-July-2020) because the script does not work with ncbi-blast-2.10.1+ (or 2.10.0+)
my $ncbi_bin = "~wrp/src/ncbi-blast-2.9.0+/bin";
if (defined($ENV{NCBI_BLAST_BIN})) {
  $ncbi_bin = $ENV{NCBI_BLAST_BIN};
}

my $psiblast_bin = "$ncbi_bin/psiblast";
my $makeblastdb_bin = "$ncbi_bin/makeblastdb";

## ncbi datatool is not part of the BLAST distribution
my $datatool_bin = "$pgm_bin/datatool -m $pgm_data/NCBI_all.asn";

my %srch_subs = ('ssearch' => \&get_ssearch_cmd,
		 'psiblast' => \&get_psiblast_cmd,
		);

my %annot_cmds = ('rpd3' => qq("\!ann_pfam28.pl --pfacc --db RPD3 --vdoms --split_over"),
		  'rpd3nv' => qq("\!ann_pfam28.pl --pfacc --db RPD3 --split_over"),
		  'rpd3nvn' => qq("\!./annot/ann_pfam28.pl --pfacc --db RPD3 --split_over --neg"),
		  'pfam' => qq("\!./annot/ann_pfam30.pl --db pfam33_qfo --vdoms --split_over --neg")
    );

($num_iter, $pssm_evalue, $srch_evalue, $dom_flag, $align_flag, $int_mask, $end_mask, $query_mask, $srch_pgm, $tmp_dir, $error_log, $annot_type, $quiet) =
  ( 5, 0.002, 5.0, 0, 0, 'none', 'none', 0, 'ssearch','',0, 0, "", 0);
($save_all, $tmp_file_list, $delete_bnd, $delete_tmp) = (0, "", 0, 0);
($prev_m89res, $m_format, $prev_sel_res, $prev_bound, $this_iter, $use_stdout) = ("","m8CB", "","", 1, 0);

my $pgm_command =  "# ".join(" ",($0,@ARGV));
print STDERR "# ",join(" ",($0,@ARGV)),"\n" if ($error_log);

GetOptions(
	   'query|sequence=s' => \$query_file,
	   'db|database=s' => \$db_file,
           'suffix|out_suffix=s' => \$suffix,
    	   'dir=s' => \$tmp_dir,
	   'pssm_evalue=f' => \$pssm_evalue,
	   'search_evalue=f' => \$srch_evalue,
	   'annot_db=s' => \$annot_type,
           'out_name=s' => \$file_out,
	   'use_stdout' => \$use_stdout,
    	   'this_iter=s' => \$this_iter,
	   'iter=i' => \$num_iter,
    	   'prev_m89res=s' => \$prev_m89res,
           'sel_res_in=s' => \$prev_sel_res,
           'sel_accs=s' => \$prev_sel_res,
           'sel_file=s' => \$prev_sel_res,
           'sel_file_in=s' => \$prev_sel_res,
	   'in_msa=s' => \$prev_msa,
	   # 'out_msa=s' => \$next_msa,
	   # 'in_hitdb=s' => \$prev_hitdb,
	   # 'out_hitdb=s' => \$next_hitdb,
	   'in_pssm=s' => \$prev_pssm,
	   # 'out_pssm=s' => \$next_pssm,
	   'prev_bounds=s' => \$prev_bound,
	   'out_bounds=s' => \$next_bound_out,
	   'num_iter|max_iter=i' => \$num_iter,
	   'dom|domain' => \$dom_flag,
	   'align' => \$align_flag,
    	   'query_seed|query_mask' => \$query_mask,
	   'int_mask|int-mask|int_seed|int-seed=s' => \$int_mask,
	   'end_mask|end-mask|end_seed|end-seed=s' => \$end_mask,
	   'm_format=s' => \$m_format,
	   'pgm=s' => \$srch_pgm,
    	   'quiet' => \$quiet,
    	   'q' => \$quiet,
    	   'silent' => \$quiet,
	   'h|?' => \$shelp,
	   "help" => \$help,
    	   "errors" => \$error_log,
	   "save_list:s" => \$tmp_file_list,  # files to save (not delete)
	   "save_tmp|save_all" => \$save_all,
	   "del_bnd|delete_bnd" => \$delete_bnd,
	   "del_tmp|del_all|delete_tmp|delete_all" => \$delete_tmp,
	  );

pod2usage(1) if $shelp;
pod2usage(exitstatus => 0, verbose => 2) if $help;

pod2usage(1) unless $query_file && -r $query_file;  # need a query
pod2usage(1) unless $db_file ;        # need a database

my %sel_hits = ();

my @del_file_ext = qw(msa psibl_out hit_db asntxt asnbin);

if ($srch_pgm =~ m/psiblast/) {
  pop(@del_file_ext);
}

if ($query_mask) {
  $int_mask='query' unless $int_mask ne 'none';
  $end_mask='query' unless $end_mask ne 'none';
}

if (!$delete_tmp && $save_all) {
  @del_file_ext = ();
  $tmp_file_list = "";
  $delete_bnd = 0
}

################
# which tmp files should be saved/deleted?
#
my %save_file_ext = ();
if ($tmp_file_list) {
  my @new_del_file_ext = ();
  for my $ext (split(/,\s*/,$tmp_file_list)) {
  # possible values are: msa, asnbin, asntxt, psibl_out, hit_db (see @del_file_ext)
    $save_file_ext{$ext} = 1;
  }

  for my $ext (@del_file_ext) {
    push @new_del_file_ext, $ext unless $save_file_ext{$ext};
  }
  @del_file_ext = @new_del_file_ext;
}

print "$pgm_command\n" unless ($quiet);

####
# generate output filenames

my ($query_pref) = ($query_file =~ m/([\w\.]+)$/);

$file_out = $query_pref unless $file_out;

my $this_file_pref = "$file_out.it$this_iter";
$this_file_pref = "$this_file_pref.$suffix" if ($suffix);
my $this_file_out = $this_file_pref;
$this_file_out = "$tmp_dir/$this_file_out" if ($tmp_dir);

my $prev_file_out = "";

####
# do the first search or use $prev_m89res
my $first_iter = 0;
my $iter_val = $this_iter;
my $search = "";

my @del_err_files = ();

unless ($prev_m89res || $prev_msa) {
  $search = $srch_subs{$srch_pgm}($query_file, $db_file, $prev_pssm);
  unless ($use_stdout) {
    log_system("$search > $this_file_out 2> $this_file_out.err");
  }
  else {
    log_system("$search 2> $this_file_out.err");
  }
  push @del_err_files, "$this_file_out.err";
  $first_iter++;
}
elsif ($prev_m89res) {
  $this_file_out = $prev_m89res;
}
elsif ($prev_msa) {
# build a PSSM, do a search, up the iteration count
  $prev_pssm  = pssm_from_msa($query_file, $prev_msa);
  $search = $srch_subs{$srch_pgm}($query_file, $db_file, $prev_pssm);  
  unless ($use_stdout) {
    log_system("$search > $this_file_out 2> $this_file_out.err");
  }
  else {
    log_system("$search 2> $this_file_out.err");
  }
  push @del_err_files, "$this_file_out.err";
  $first_iter++;
}

my ($this_pssm, $this_bound_out) = ("","");

# now have necessary files for next iteration

for (my $it=$first_iter; $it < $num_iter; $it++) {

  ####
  # build the PSSM for the current search

  my ($this_pssm, $this_bound_out) = build_msa_pssm($query_file, $this_file_out, $prev_bound, $prev_sel_res, $m_format);
  $prev_file_out = $this_file_out;
  $prev_sel_res = "";

  $iter_val = $this_iter + $it;

  ####
  # build filename for this iteration
  $this_file_pref = $this_file_out = "$file_out.it$iter_val";
  $this_file_out = "$this_file_pref.$suffix" if ($suffix);
  $this_file_out = "$tmp_dir/$this_file_out" if ($tmp_dir);

  $search = $srch_subs{$srch_pgm}($query_file, $db_file, $this_pssm);
  unless ($use_stdout) {
    log_system("$search > $this_file_out 2> $this_file_out.err");
  }
  else {
    log_system("$search 2> $this_file_out.err");
  }
  push @del_err_files, "$this_file_out.err";

  ####
  # here, we are done with previous .msa, .asntxt, .asnbin, etc files.  Delete them if desired

  if (@del_file_ext) {
    my @del_file_list = ();
    for my $ext (@del_file_ext) {
      push @del_file_list, "$prev_file_out.$ext";
    }
    log_system("rm ".join(" ",@del_file_list));
  }
  log_system("rm $prev_bound") if ($delete_bnd);

  if (has_converged($prev_bound, $this_bound_out)) {
    print STDERR "$0 $srch_pgm $query_file $db_file converged ($iter_val iterations)\n" unless ($quiet);
    last;
  }
  $prev_bound = $this_bound_out
}

if (@del_err_files) {
  log_system("rm ".join(" ",@del_err_files));
}

log_system("rm $prev_bound") if ($delete_bnd);
unless ($quiet) {
  print STDERR "$0 $srch_pgm $query_file $db_file finished ($num_iter iterations)\n";
}

################
# log_system()
# run system on string, logging first if $error_log
#
sub log_system {

  my ($cmd) = @_;

  print STDERR "# $cmd\n" if $error_log;
  system($cmd);
}

################
# sub get_ssearch_cmd()
# builds an ssearch command line with query, db, and pssm
#
sub get_ssearch_cmd {
  my ($query_file, $db_file, $pssm_file) = @_;

  my $mf_arg = $m_format;
  $mf_arg =~ s/^m//;
  $mf_arg =~ s/\+/ /;

  my $search_cmd = qq($ssearch_bin -S -E "$srch_evalue 0" -s BP62 -m $mf_arg);

  if ($annot_type) {
    $search_cmd .= qq( -V $annot_cmds{$annot_type});
  }
  if ($pssm_file) {
    $search_cmd .= qq( -P "$pssm_file 2");
  }

  $search_cmd .= qq( $query_file $db_file);

  return $search_cmd;
}

################
# sub get_psiblast_cmd()
# builds an ssearch command line with query, db, and pssm
#
sub get_psiblast_cmd {
  my ($query_file, $db_file, $pssm_file) = @_;

  my $search_cmd = qq($psiblast_bin -num_threads 8 -max_target_seqs 5000 -outfmt '7 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore score btop' -inclusion_ethresh $pssm_evalue -evalue $srch_evalue -num_iterations 1 -db $db_file);
  if ($pssm_file) {
    $search_cmd .= qq( -in_pssm $pssm_file);
#    $search_cmd .= qq( -comp_based_stats 0);
  }
  else {
    $search_cmd .= qq( -query $query_file);
  }

  return $search_cmd;
}

################
# sub build_msa_pssm()
#
# given query, search output file ($this_file_out), prev_boundary_file
# uses m89_btop_msa2.pl to generate PSSM in .asntxt or .asnbin format, also bound_file_out if $align_flag
# (later - optionally deletes intermediate files)
#
# always produce a $bound_file_out file to test for convergence
#
sub build_msa_pssm {
  my ($query_file, $this_file_out,$prev_bound_in, $prev_sel_in, $m_format) = @_;

  my ($this_msa, $this_hit_db, $this_pssm_asntxt, $this_pssm_asnbin, $this_psibl_out, $this_bound_out) =
    ("$this_file_out.msa",
     "$this_file_out.hit_db",
     "$this_file_out.asntxt",
     "$this_file_out.asnbin",
     "$this_file_out.psibl_out",
     "$this_file_out.bnd_out",
    );

  my $blastdb_err = "$this_file_out.mkbldb_err";
  my $aln2msa_cmd = qq($align2msa_lib --query $query_file --masked_lib_out=$this_hit_db);

  if ($m_format) {
    $aln2msa_cmd .= qq( --m_format $m_format);
  }

  if ($prev_sel_in) {
    $aln2msa_cmd .= qq( --sel_file_in $prev_sel_in);
  }
  else {
    $aln2msa_cmd .= qq( --evalue $pssm_evalue);
  }

  if ($int_mask) {
    $aln2msa_cmd .= qq( --int_mask_type $int_mask);
  }

  if ($end_mask) {
    $aln2msa_cmd .= qq( --end_mask_type $end_mask);
  }

  if ($dom_flag) {
    $aln2msa_cmd .= qq( --domain);
  }

  if ($align_flag && $prev_bound_in) {
    $aln2msa_cmd .= qq( --bound_file_in $prev_bound_in);
  }

  # always produce this file to check for convergence
  $aln2msa_cmd .= qq( --bound_file_out $this_bound_out);

  log_system("$aln2msa_cmd $this_file_out > $this_msa");

  my $makeblastdb_cmd = "$makeblastdb_bin -in $this_hit_db -dbtype prot -parse_seqids > $blastdb_err";

  log_system($makeblastdb_cmd);

  my $buildpssm_cmd = "$psiblast_bin -max_target_seqs 5000 -outfmt 7 -inclusion_ethresh 100.0 -in_msa $this_msa -db $this_hit_db -out_pssm $this_pssm_asntxt -num_iterations 1 -save_pssm_after_last_round";

  log_system("$buildpssm_cmd  > $this_psibl_out 2> $this_psibl_out.err");

  log_system("rm $this_hit_db.p* $blastdb_err") unless $error_log;

  # remove uninformative error logs
  log_system("rm $this_psibl_out.err") unless $error_log;

  unless ($srch_pgm eq 'psiblast') {
    my $asn2asn_cmd = "$datatool_bin -v $this_pssm_asntxt -e $this_pssm_asnbin";
    log_system($asn2asn_cmd);
    return ($this_pssm_asnbin, $this_bound_out);
  }
  else {
    return ($this_pssm_asntxt, $this_bound_out);
  }
}

################
# pssm_from_msa()
#
# given query, --in_msa Clustal MSA
# use psiblast to generate PSSM in .asntxt or .asnbin format
# (later - optionally deletes intermediate files)
#
# always produce a $bound_file_out file to test for convergence
#
sub pssm_from_msa {
  my ($query_file, $msa_file) = @_;

  my $this_file_out = $query_file;

  my ($this_hit_db, $this_pssm_asntxt, $this_pssm_asnbin, $this_psibl_out, $this_bound_out) =
    ("$this_file_out.hit_db",
     "$this_file_out.asntxt",
     "$this_file_out.asnbin",
     "$this_file_out.psibl_out",
     "$this_file_out.bnd_out",
    );

  my $blastdb_err = "$this_file_out.mkbldb_err";
  ## should not need this, but may need to convert in_msa file to fasta file for equivalence to build_msa_pssm()
  my $clus2fa_cmd = qq($clustal2fasta $msa_file > $this_hit_db);

  log_system($clus2fa_cmd);

  my $makeblastdb_cmd = "$makeblastdb_bin -in $this_hit_db -dbtype prot -parse_seqids > $blastdb_err";
  log_system($makeblastdb_cmd);

  my $buildpssm_cmd = "$psiblast_bin -max_target_seqs 5000 -outfmt 7 -inclusion_ethresh 100.0 -in_msa $msa_file -db $this_hit_db -out_pssm $this_pssm_asntxt -num_iterations 1 -save_pssm_after_last_round";

  log_system("$buildpssm_cmd  > $this_psibl_out 2> $this_psibl_out.err");

  log_system("rm $this_hit_db.p* $blastdb_err") unless $error_log;

  # remove uninformative error logs
  log_system("rm $this_psibl_out.err") unless $error_log;

  unless ($srch_pgm eq 'psiblast') {
    my $asn2asn_cmd = "$datatool_bin -v $this_pssm_asntxt -e $this_pssm_asnbin";
    log_system($asn2asn_cmd);
    return ($this_pssm_asnbin);
  }
  else {
    return ($this_pssm_asntxt);
  }
}

################
# sub has_converged()
# reads two boundary files and compares accessions
#
sub has_converged {
  my ($file1, $file2) = @_;

  return 0 unless ($file1 && $file2);

  my @f1_names = ();
  my @f2_names = ();

  open (my $fd1, '<', $file1) || die "cannot read $file1";
  while (my $line = <$fd1>) {
    chomp($line);
    my @fields = split(/\t/,$line);
    push @f1_names, $fields[0];
  }
  close $fd1;

  open (my $fd2, '<', $file2) || die "cannot read $file2";
  while (my $line = <$fd2>) {
    chomp($line);
    my @fields = split(/\t/,$line);
    push @f2_names, $fields[0];
  }
  close $fd2;

  # check for same length
  return 0  if (scalar(@f1_names) != scalar(@f2_names));

  @f1_names = sort @f1_names;
  @f2_names = sort @f2_names;

  for (my $i=0; $i < scalar(@f1_names); $i++) {
    return 0 if ($f1_names[$i] ne $f2_names[$i]);
  }
  return 1;
}

__END__

=pod

=head1 NAME

psisearch2_msa.pl

=head1 SYNOPSIS

 psisearch2_msa.pl --query q_file --db db_file --pgm ssearch|psiblast --num_iter 5

=head1 OPTIONS

 -h	short help
 --help include description

 --query query file  (also --sequence)
 --db    database file (--database)
 --pgm   [ssearch] program used for searching, ssearch or psiblast
 --num_iter/iter [5] maximum number of iterations (--max_iter)


 --this_iter [0] iteration number for tracking
 --pssm_evalue  [0.002] threshold for inclusion in PSSM
 --search_evalue [5.0] threshold for inclusion in search display
 --annot_db [null] (rpd3, rpd3nv, rpd3nvn, pfam)

 --dir   working directory and location of output
 --out_name/--suffix  result file is "out_name.it#.suffix"
 --in_msa/--out_msa   [not implemented] MSA used to build PSSM, requires --in_hitdb
 --in_hitdb/--out_hitdb [not implemented] used to build PSSM
 --in_pssm/--out_pssm [--out_pssm not implemented]
 --in_bounds/--out_bounds used to control alignment boundaries for PSSM
 --int_mask/--end_mask  none|query|random -  values embeeded in library sequences based on gaps in MSA
 --delete_all remove all tmp files (also --delete_tmp, --del_all, --del_tmp)
 --delete_bnd remove boundary file (included with --delete_all, but not deleted by default)
 --save_all save temporary files (.asnbin, .asntxt, .msa, .psiblout, .hit_db)
 --save_list: comma delimited list of file extensions (above) to save

=head1 DESCRIPTION

C<psisearch2_msa.pl> can perform one, or multiple, successive iterations of C<psiblast> or
C<ssearch36> using different query seeding strategies to reduce PSSM contamination
from alignment over-extension.  C<psisearch2_msa.pl> uses the
C<m89_btop_msa2.pl> program to read BTOP formatted output from
C<psiblast> or C<ssearch36> and produce both a multiple sequence
alignment (MSA) and a fasta formatted custom database of the sequences in the
MSA.  C<psiblast> then produces a PSSM From the MSA and custom database.

Different strategies to reduce PSSM contamination from alignment
overextension can be specified using the C<--int_mask>, C<--end_mask>,
C<--align>, and C<--domain> options.  If C<--int_mask> and
C<--end_mask> are not specified (or set to "none"), then the PSSM is
generated by aligning the MSA with the sequence residues that were
aligned in the previous search. C<--int_mask> and C<--int_mask> are
set to "query", then any gaps in the MSA are filled with the
corresponding aligned residue from the query sequence.  In our
experience, this dramatically reduces alignment over-extension and
false-positives.

=head1 AUTHORS

William R. Pearson (wrp@virginia.edu) and Weizhong Li (wli@ebi.ac.uk)

=cut
