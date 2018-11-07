#!/usr/bin/env perl

################################################################
# copyright (c) 2017,2018 by William R. Pearson and The Rector &
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
# annot_blast_btop2.pl --query query.file --ann_script ann_pfam_www.pl --include_doms blast_tab_btop_file
################################################################
# annot_blast_btop2.pl associates domain annotation information and
# subalignment scores with a blast tabular (-outfmt 6 or -outfmt 7)
# file that contains the raw score and the BTOP alignment encoding
# This file can be generated from "blastp/n" or "blast_formatter"
# using the command:
#   blast_formatter -archive blast_output.asn -outfmt '7 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore score btop'  > blast_output.tab_annot
#
# If the BTOP field or query_file is not available, the script
# produces domain content without sub-alignment scores.
################################################################
## 4-Nov-2018
# add --include_doms, which adds a new field with the coordinates of
# the domains in the protein (independent of alignment)
#
################################################################
## 21-July-2018
# include sequence length (actually alignment end) to produce NODOM's (no NODOM's without length).
#
################################################################
## 13-Jan-2017
# modified to provide query/subject coordinates and identities if no
# query sequence -- does not decrement for reverse-complement fastx/blastx DNA
################################################################
## 16-Nov-2015
# modify to allow multi-query blast searches
################################################################
## 19-Dec-2015
# add -q_annot_script to annotate query sequence
#

use warnings;
use strict;
use IPC::Open2;
use Pod::Usage;
use Getopt::Long;
use File::Temp qw/ tempfile /;

# use Data::Dumper;

# read lines of the form:
# gi|121694|sp|P20432.1|GSTT1_DROME	gi|121694|sp|P20432|GSTT1_DROME	100.00	209	0	0	1	209	1	209	6e-156	433	1113	209
# gi|121694|sp|P20432.1|GSTT1_DROME	gi|1170090|sp|P04907|GSTF3_MAIZE	26.77	198	123	7	4	185	6	197	2e-08	51.2	121	FL1YG ... 1NKRA1YW1
# gi|121694|sp|P20432.1|GSTT1_DROME	gi|81174731|sp|P0ACA5|SSPA_ECO57	39.66	58	32	2	43	100	49	103	8e-06	43.9	102	EDFLLI ... V-I-NEQS3FM
# gi|121694|sp|P20432.1|GSTT1_DROME	gi|121695|sp|P12653|GSTF1_MAIZE	27.62	181	107	7	32	203	34	199	9e-05	40.8	94	LI1LF ... N-1AS1CLLM1

# and report the domain content ala -m 8CC

my ($matrix, $ann_script, $q_ann_script, $show_raw, $shelp, $help) = ("BLOSUM62", "", "", 0, 0, 0);
my ($have_qslen, $dom_info, $sub2query) = (0,0,0);		# blast tabular file has sseqid sseqlen qseqid qseqlen
my ($query_lib_name) = ("");  # if $query_lib_name, do not use $query_file_name
my ($out_field_str) = ("");
my $query_lib_r = 0;

my @blosum62 = ();
my @blosum62_diag = ();
my %aa_map = ();
my ($g_open, $g_ext) = (-11, -1);
init_blosum62();

GetOptions(
    "matrix:s" => \$matrix,
    "ann_script|script:s" => \$ann_script,
    "q_ann_script|q_script:s" => \$q_ann_script,
    "have_qslen|have_sqlen!" => \$have_qslen,
    "domain_info|dom_info!" => \$dom_info,
    "sub2query!" => \$sub2query,
    "query:s" => \$query_lib_name,
    "query_file:s" => \$query_lib_name,
    "query_lib:s" => \$query_lib_name,
    "out_fields:s" => \$out_field_str,
    "raw_score" => \$show_raw,
    "h|?" => \$shelp,
    "help" => \$help,
    );

pod2usage(1) if $shelp;
pod2usage(exitstatus => 0, verbose => 2) if $help;
unless (-f STDIN || -p STDIN || @ARGV) {
 pod2usage(1);
}

if ($query_lib_name) {
  $query_lib_r = parse_query_lib($query_lib_name);
}

my @tab_fields = qw(q_seqid s_seqid percid alen mismatch gopen q_start q_end s_start s_end evalue bits score BTOP);

if ($have_qslen) {
  @tab_fields = qw(q_seqid q_len s_seqid s_len percid alen mismatch gopen q_start q_end s_start s_end evalue bits score BTOP);
}

# the fields that are displayed are listed here.  By default, all fields except score and BTOP are displayed.
my @out_tab_fields = @tab_fields[0 .. $#tab_fields-1];

if ($show_raw) {
  push @out_tab_fields, "raw_score";
}

if ($out_field_str) {
  @out_tab_fields = split(/\s+/,$out_field_str);
}

my @header_lines = ();

# need outer loop to enable multiple queries
while (1) {

  my $next_line = "";
  my $have_data = 0;

  my @hit_list = ();
  my @q_hit_list = ();

  while (my $line = <>) {
    if ($line =~ /^#/) {
      if ($have_data) {
	$next_line = $line;
	$have_data = 0;
	last;
      } else {
	push @header_lines, $line;
      }
      next;
    }

    $have_data = 1;
    my %hit_data = ();
    chomp $line;
    next unless $line;
    @hit_data{@tab_fields} = split(/\t/,$line);

    push @hit_list, \%hit_data;
  }

  # get the query annotations
  if ($q_ann_script) {
      $q_ann_script =~ s/\+/ /g;
  }

  if ($q_ann_script && -x (split(/\s+/,$q_ann_script))[0]) {
    # get the domains for the q_seqid using --q_ann_script
    #
    my ($Reader, $Writer);
    my $pid = open2($Reader, $Writer, $q_ann_script);
    my $hit = $hit_list[0];

    my $q_seq_len = scalar(@{$query_lib_r->{$hit->{q_seqid}}});
    print $Writer $hit->{q_seqid},"\t",$q_seq_len,"\n";
    close($Writer);

    push @q_hit_list,{ s_seq_id=> $hit->{q_seqid}, s_end=> $q_seq_len};

    read_annots($Reader, \@q_hit_list, 0);

    waitpid($pid, 0);
  }

  # get the subject annotations
  if ($ann_script) {
      $ann_script =~ s/\+/ /g;
  }

  if ($ann_script && -x (split(/\s+/,$ann_script))[0]) {
    # get the domains for each s_seqid using --ann_script
    #
    # this does not work currently because only one accession is sent.
    # For mulitple hits, I need to make a tmp_file.

    my ($Reader, $Writer);
    my $pid = open2($Reader, $Writer, $ann_script);

    for my $hit (@hit_list) {
#      print STDERR  $hit->{s_seqid},"\t", $hit->{s_end},"\n";
	#      print $Writer $hit->{s_seqid},"\t", $hit->{s_end},"\n";
      my $s_len = 100000;
      if ($have_qslen) {
	$s_len = $hit->{s_len};
      }
      print $Writer $hit->{s_seqid},"\t", $s_len,"\n";
    }
    close($Writer);

    read_annots($Reader, \@hit_list, 1);

    waitpid($pid, 0);
  }

  for my $line (@header_lines) {
    print $line;
  }
  @header_lines = ($next_line);

  # now get query sequence if available

  if ($sub2query && scalar(@q_hit_list)==0) {
    # copy the information from $hit_list
    for my $tmp_hit ( @hit_list ) {
      if ($tmp_hit->{q_seqid} eq $tmp_hit->{s_seqid}) {
	my %tmp_q_hit = (s_seq_id=> $tmp_hit->{q_seqid}, s_end=> $tmp_hit->{s_len});

	$tmp_q_hit{'domains'} = [];
	for my $dom ( @{$tmp_hit->{domains}} ) {
	  my %new_dom = map { $_ => $dom->{$_} } keys(%$dom);
	  $new_dom{target} = 0;
	  push @{$tmp_q_hit{'domains'}}, \%new_dom;
	}

	$tmp_q_hit{'sites'} = [];
	for my $site ( @{$tmp_hit->{sites}} ) {
	  my %new_site = map { $_ => $site->{$_} } keys(%$site);
	  $new_site{target} = 0;
	  push @{$tmp_q_hit{'sites'}}, \%new_site;
	}
	push @q_hit_list,\%tmp_q_hit;
	last;
      }
    }
  }

  my $q_hit = $q_hit_list[0];

  for my $hit (@hit_list) {
    my @list_covered = ();

    # If I have an encoded aligment {BTOP} and a query sequence $query_lib_r && $query_lib_r->{$hit->{q_seqid}}
    # then I can calculate sub-alignment scores
    if (defined($hit->{BTOP}) && $query_lib_r && $query_lib_r->{$hit->{q_seqid}}) {

      $hit->{raw_score} = 0;  # initialize in case no domains and raw_score requested
      # calculate sub-alignment scores in subject/library coordinates
      if (defined($hit->{domains}) && scalar(@{$hit->{domains}})) {
	($hit->{raw_score}, $hit->{aligned_domains_r}) = 
	  sub_alignment_score($query_lib_r->{$hit->{q_seqid}},
			      $hit, \@blosum62, \@blosum62_diag, $hit->{domains}, 1);
      }

      if (defined($hit->{sites}) && scalar(@{$hit->{sites}})) {
	$hit->{aligned_sites_r} = site_align($query_lib_r->{$hit->{q_seqid}},
					     $hit, \@blosum62, $hit->{sites}, 1);
      }

      # calculate sub-alignment scores in query coordinates
      if (defined($q_hit->{domains}) && scalar(@{$q_hit->{domains}})) {
	($hit->{raw_score}, $hit->{q_aligned_domains_r}) = 
	  sub_alignment_score($query_lib_r->{$hit->{q_seqid}},
			      $hit, \@blosum62, \@blosum62_diag, $q_hit->{domains}, 0);
      }

      if (defined($q_hit->{sites}) && scalar(@{$q_hit->{sites}})) {
	$hit->{q_aligned_sites_r} =  site_align($query_lib_r->{$hit->{q_seqid}},
						$hit, \@blosum62, $q_hit->{sites}, 0);
      }
    }
    elsif (defined($hit->{BTOP})) {
      if (defined($hit->{domains}) && scalar(@{$hit->{domains}})) {
	$hit->{aligned_domains_r} = 
	  sub_alignment_pos($hit, $hit->{domains}, 1);
      }
    }
    else {		   # no alignment info, can provide domain overlap, and subject coordinates
      $hit->{raw_score} = 0;
      for my $dom_r (@{$hit->{domains}}) {
	next if $dom_r->{d_end} < $hit->{s_start}; # before start
	last if $dom_r->{d_pos} > $hit->{s_end}; # after end

	if ($dom_r->{d_pos} <= $hit->{s_end} && $dom_r->{d_end} >= $hit->{s_start}) {
	  push @list_covered, $dom_r->{descr};
	}
      }
    }


    ################
    ## final output display

    print join("\t",@{$hit}{@out_tab_fields}); # show fields from original blast tabular file

    my $merged_annots_r = merge_annots($hit);	# merge the four possible annotation lists into one.

    if (scalar(@$merged_annots_r)) { # show subalignment scores if available
      print "\t";
      print format_annot_info($hit, $merged_annots_r);
      if ($dom_info) {
	  print "\t",format_dom_info($q_hit->{domains}, $hit->{domains});
      }
    }
    elsif (@list_covered) {	# otherwise show domain content
      print "\t",join(";",@list_covered);
      if ($dom_info) {
	  print "\t",format_dom_info($q_hit->{domains}, $hit->{domains});
      }
    }

    print "\n";
  }

  # for my $line (@footer_lines) {
  #   print $line;
  # }
  # @footer_lines = ();

  last if eof(ARGV);
}

for my $line (@header_lines) {
  print $line;
}

################
# read_annots (\@hit_list)
# input: $hit_entry->{s_seq_id, etc}, $target
# output: modified $hit_entry->{domains}
#         modified $hit_entry->{sites}

sub read_annots {
  my ($Reader, $hit_list_r, $target) = @_;

  my $current_domain = "";
  my $hit_ix = 0;
  my @hit_domains = ();
  my @hit_sites = ();

  while (my $line = <$Reader>) {
    next if $line=~ m/^=/;
    chomp $line;

#    print STDERR "$line\n";

    # check for header
    if ($line =~ m/^>/) {
      if ($current_domain) {  # previous domains/sites have already been found and parsed
	if ($hit_list_r->[$hit_ix]{s_seqid} eq $current_domain) {
	    $hit_list_r->[$hit_ix]{domains} = [ @hit_domains ];  # previous domains
	    $hit_list_r->[$hit_ix]{sites} = [ @hit_sites ];      # previous sites
	    $hit_ix++;
	  } else {
	    warn "phase error: $current_domain != $hit_list_r->[$hit_ix]{s_seqid}";
	  }
      }
      @hit_domains = ();   # current domains
      @hit_sites = ();     # current sites
      $current_domain = (split(/\s+/,$line))[0];
      $current_domain =~ s/^>//;
    } else {			# check for data
      my %annot_info = (target=>$target);
      my @a_fields = split(/\t/,$line);
      if ($a_fields[1] eq '-') {
	@annot_info{qw(d_pos type d_end descr)} = @a_fields;
	$annot_info{descr} =~ s/ :(\d+)$/~$1/;
	push @hit_domains, \%annot_info;   # current
      }
      else {
	@annot_info{qw(d_pos type d_val descr)} = @a_fields;
	$annot_info{'d_end'} = $annot_info{'d_pos'};
	push @hit_sites, \%annot_info;   # current
      }
    }
  }
  close($Reader);

  $hit_list_r->[$hit_ix]{domains} = \@hit_domains;
  $hit_list_r->[$hit_ix]{sites} = \@hit_sites;

  # clean up NODOMs in {domains}
  for my $hit ( @$hit_list_r ) {
    # clean-up last NODOM if < 10
    my $tmp_domains = $hit->{domains};
    next unless (scalar(@{$tmp_domains}));
    my ($last_dom, $left_coord) = ($tmp_domains->[-1], $hit->{s_end});
    if ($last_dom->{descr} =~ m/^NODOM/ && (($left_coord - $last_dom->{d_pos} + 1) < 10)) {
      pop @$tmp_domains;
    }
  }
}

# input: a blast BTOP string of the form: "1VA160TS7KG10RK27"
# returns a list_ref of tokens: (1, "VA", 60, "TS", 7, "KG, 10, "RK", 27)
#
sub decode_btop {
  my ($btop_str) = @_;

  my @tokens = split(/(\d+)/,$btop_str);

  shift @tokens unless $tokens[0];

  my @out_tokens = ();

  for my $token (@tokens) {
    if ($token =~ m/^\d+$/) {
      push @out_tokens, $token
    }
    else {
      my @mis_tokens = split(/(..)/,$token);
      for my $mis (@mis_tokens) {
	if ($mis) {push @out_tokens, $mis};
      }	
    }
  }

  return \@out_tokens;
}

sub parse_query_lib {
  my ($query_file) = @_;

  my %query_seqs = ();

  open(my $qfd, $query_file);
  {				# local scope for $/
    local $/ = "\n>";

    while (my $entry = <$qfd>) {  # returns an entire fasta entry
      chomp $entry;
      my ($header, $sequence) = ($entry =~ m/^>?           # ^> only in first entry
                                             ( [^\n]* ) \n # header line
                                             ( .*     )    # the sequence
                                            /osx);    # optimize, multiline, commented
      $sequence =~ s/[^A-Za-z\*]//g;    # remove everything but letters
      $sequence = uc($sequence);
      $header =~ s/\s.*$//;
      my @seq = split(//,$sequence);
      unshift @seq,"";	# @seq is now 1-based
      $query_seqs{$header} =  \@seq;
    }
  }
  return \%query_seqs;
}

sub parse_query_file {
  my ($query_file) = @_;

  my $seq_data = "";

  open(my $qfd, $query_file);
  while (my $line = <$qfd>) {
    next if $line =~ m/^>/;
    next if $line =~ m/^;/;
    chomp $line;
    $line =~ s/[^A-Za-z\*]//g;
    $seq_data .= $line
  }

  $seq_data = uc($seq_data);

  my @seq = split(//,$seq_data); 

  return \@seq;
}

sub init_blosum62 {

  my @ncbi_blaa = qw(    A  R  N  D  C  Q  E  G  H  I  L  K  M  F  P  S  T  W  Y  V  B  Z  X  * );

  $blosum62[ 0] = [ qw(  4 -1 -2 -2  0 -1 -1  0 -2 -1 -1 -1 -1 -2 -1  1  0 -3 -2  0 -2 -1  0 -4) ]; # A
  $blosum62[ 1] = [ qw( -1  5  0 -2 -3  1  0 -2  0 -3 -2  2 -1 -3 -2 -1 -1 -3 -2 -3 -1  0 -1 -4) ]; # R
  $blosum62[ 2] = [ qw( -2  0  6  1 -3  0  0  0  1 -3 -3  0 -2 -3 -2  1  0 -4 -2 -3  3  0 -1 -4) ];
  $blosum62[ 3] = [ qw( -2 -2  1  6 -3  0  2 -1 -1 -3 -4 -1 -3 -3 -1  0 -1 -4 -3 -3  4  1 -1 -4) ];
  $blosum62[ 4] = [ qw(  0 -3 -3 -3  9 -3 -4 -3 -3 -1 -1 -3 -1 -2 -3 -1 -1 -2 -2 -1 -3 -3 -2 -4) ];
  $blosum62[ 5] = [ qw( -1  1  0  0 -3  5  2 -2  0 -3 -2  1  0 -3 -1  0 -1 -2 -1 -2  0  3 -1 -4) ];
  $blosum62[ 6] = [ qw( -1  0  0  2 -4  2  5 -2  0 -3 -3  1 -2 -3 -1  0 -1 -3 -2 -2  1  4 -1 -4) ];
  $blosum62[ 7] = [ qw(  0 -2  0 -1 -3 -2 -2  6 -2 -4 -4 -2 -3 -3 -2  0 -2 -2 -3 -3 -1 -2 -1 -4) ];
  $blosum62[ 8] = [ qw( -2  0  1 -1 -3  0  0 -2  8 -3 -3 -1 -2 -1 -2 -1 -2 -2  2 -3  0  0 -1 -4) ];
  $blosum62[ 9] = [ qw( -1 -3 -3 -3 -1 -3 -3 -4 -3  4  2 -3  1  0 -3 -2 -1 -3 -1  3 -3 -3 -1 -4) ];
  $blosum62[10] = [ qw( -1 -2 -3 -4 -1 -2 -3 -4 -3  2  4 -2  2  0 -3 -2 -1 -2 -1  1 -4 -3 -1 -4) ];
  $blosum62[11] = [ qw( -1  2  0 -1 -3  1  1 -2 -1 -3 -2  5 -1 -3 -1  0 -1 -3 -2 -2  0  1 -1 -4) ];
  $blosum62[12] = [ qw( -1 -1 -2 -3 -1  0 -2 -3 -2  1  2 -1  5  0 -2 -1 -1 -1 -1  1 -3 -1 -1 -4) ];
  $blosum62[13] = [ qw( -2 -3 -3 -3 -2 -3 -3 -3 -1  0  0 -3  0  6 -4 -2 -2  1  3 -1 -3 -3 -1 -4) ];
  $blosum62[14] = [ qw( -1 -2 -2 -1 -3 -1 -1 -2 -2 -3 -3 -1 -2 -4  7 -1 -1 -4 -3 -2 -2 -1 -2 -4) ];
  $blosum62[15] = [ qw(  1 -1  1  0 -1  0  0  0 -1 -2 -2  0 -1 -2 -1  4  1 -3 -2 -2  0  0  0 -4) ];
  $blosum62[16] = [ qw(  0 -1  0 -1 -1 -1 -1 -2 -2 -1 -1 -1 -1 -2 -1  1  5 -2 -2  0 -1 -1  0 -4) ];
  $blosum62[17] = [ qw( -3 -3 -4 -4 -2 -2 -3 -2 -2 -3 -2 -3 -1  1 -4 -3 -2 11  2 -3 -4 -3 -2 -4) ];
  $blosum62[18] = [ qw( -2 -2 -2 -3 -2 -1 -2 -3  2 -1 -1 -2 -1  3 -3 -2 -2  2  7 -1 -3 -2 -1 -4) ];
  $blosum62[19] = [ qw(  0 -3 -3 -3 -1 -2 -2 -3 -3  3  1 -2  1 -1 -2 -2  0 -3 -1  4 -3 -2 -1 -4) ];
  $blosum62[20] = [ qw( -2 -1  3  4 -3  0  1 -1  0 -3 -4  0 -3 -3 -2  0 -1 -4 -3 -3  4  1 -1 -4) ];
  $blosum62[21] = [ qw( -1  0  0  1 -3  3  4 -2  0 -3 -3  1 -1 -3 -1  0 -1 -3 -2 -2  1  4 -1 -4) ];
  $blosum62[22] = [ qw(  0 -1 -1 -1 -2 -1 -1 -1 -1 -1 -1 -1 -1 -1 -2  0  0 -2 -1 -1 -1 -1 -1 -4) ];
  $blosum62[23] = [ qw( -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4  1) ];

  die "blosum62 length mismatch $#blosum62 != $#ncbi_blaa" if (scalar(@blosum62) != scalar(@ncbi_blaa));

  for (my $i=0; $i < scalar(@ncbi_blaa); $i++) {
    $aa_map{$ncbi_blaa[$i]} = $i;
    $blosum62_diag[$i] = $blosum62[$i][$i];
  }

  ($g_open, $g_ext) = (-11, -1);
}

# given: (1) a query sequence; (2) an encoded alignment; (3) a scoring matrix
# calculate a score

sub alignment_score {
  my ($query_r, $query_start, $btop_align_r, $matrix_2d) = @_;

  my ($gap0, $gap1) = (0,0);

  my $qix = $query_start-1; # start from zero

  my ($score, $m_score) = 0;
  my ($seq0, $seq1) = ("","");
  for my $btop (@{$btop_align_r}) {
    if ($btop =~ m/^\d+$/) {  # matching query sequence, add it up
      for (my $i=0; $i < $btop; $i++) {
	$score += $matrix_2d->{$query_r->[$qix]}{$query_r->[$qix]};
	$qix++;
      }
    }
    else {
      ($seq0, $seq1) = split(//,$btop);
      if ($btop=~ m/\-/) {
	if ($seq0 eq '-') {
	  if ($gap0) { $score += $g_ext;}
	  else { $score += $g_open+$g_ext;}
	  $gap0 = 1;
	}
	else {
	  if ($gap1) { $score += $g_ext;}
	  else { $score += $g_open+$g_ext;}
	  $gap1 = 1;
	  $qix++;
	}
      }
      else {
	$score += $matrix_2d->{$seq0}{$seq1};
	$gap1=$gap0 = 0;
	$qix++;
      }
    }
  }
  return $score;
}

################################################################
# sub_alignment_score()
# input: $query_r : a query sequence;
#        $hit_r->{BTOP} : an encoded alignment;
#        $matrix_2d, $matrix_diag : a scoring matrix
#        $domain_r : domain boundaries in query (target=0) or subject (target=1)
#        $target : 0=query, 1=target
#
# calculate a score
# updates $domain_r in place with new values:
# domain_r->[]->{ident} (as fraction identical),
#             ->{score} --matrix raw similarity score
#             ->{qa_start,qa_end}  domain boundaries in query
#             ->{sa_start, sa_end} domain boundaries in subject
#
sub sub_alignment_score {
  my ($query_r, $hit_r, $matrix_2d, $matrix_diag, $domain_r, $target) = @_;

  return (0, $domain_r) unless ($domain_r && scalar(@$domain_r));

  my $btop_enc_r = decode_btop($hit_r->{BTOP});

  my ($gap0, $gap1) = (0,0);

  my @active_dom_list = ();
  my @aligned_domains = ();

  my $left_active_end = $domain_r->[-1]->{d_end}+1;	# as far right as possible
  my $left_align_end = $hit_r->{q_end};
  if ($target) {
    $left_align_end = $hit_r->{s_end};
  }

  if ($left_active_end > $left_align_end ) {
    $left_active_end = $left_align_end ;
  }

  my ($q_start, $s_start, $h_start, $h_end) = @{$hit_r}{qw(q_start s_start s_start s_end)};
  my ($qix, $six)  = ($q_start, $s_start); # $qix now starts from 1, like $six;

  my $ds_ix = \$six;	# use to track the subject position
  # reverse coordinate names if $target==0
  unless ($target) {
    $ds_ix = \$qix;	# track query position
    $h_start = $hit_r->{q_start};
    $h_end = $hit_r->{q_end};
  }

  my ($score, $m_score) = 0;
  my ($seq0, $seq1) = ("","");

  # find the first overlapping domain
  my ($dom_ix, $dom_nx) = (0,scalar(@$domain_r));
  my $dom_r = $domain_r->[0];

  # skip over domains that do not overlap alignment
  # capture first domain that alignment overlaps
  for ($dom_ix=0; $dom_ix < $dom_nx; $dom_ix++) {
    if ($domain_r->[$dom_ix]->{d_end} >= $h_start) {  # if {d_end} < $_start, cannot overlap
      $dom_r = $domain_r->[$dom_ix];
      if ($dom_r->{d_pos} <= $h_start) {  # {d_pos} is less, {d_end} is greater, overlap
	push @aligned_domains, $dom_r;
	$left_active_end = push_annot_match(\@active_dom_list, $dom_r, $q_start, $s_start, 0, 0);
      }
      else { last; }
    }
  }

  my ($dom_score, $id_cnt) = (0,0);

  for my $btop (@{$btop_enc_r}) {

    if ($btop =~ m/^\d+$/) {  # matching query sequence, add it up
      for (my $i=0; $i < $btop; $i++) {  # $i is used to count through BTOP, not to index anything.

	my $seq0_map = $aa_map{'X'};
	unless ($query_r->[$qix]) {
	  warn "qix: $qix out of range";
	}
	else {
	  $seq0_map = $aa_map{$query_r->[$qix]} if exists($aa_map{$query_r->[$qix]});
#	  print "$qix:$six : ",$query_r->[$qix],"\n";
	}

	$m_score = $matrix_diag->[$seq0_map];
	$score += $m_score;

	if ($dom_ix < $dom_nx && $$ds_ix == $dom_r->{d_pos}) {
	  push @aligned_domains, $dom_r;
	  $left_active_end = push_annot_match(\@active_dom_list, $dom_r, $qix, $six, $id_cnt, $dom_score);
	  $dom_ix++;
	  $dom_r = $domain_r->[$dom_ix];
	  ($dom_score, $id_cnt) = (0,0);
	}

	if (@active_dom_list) {
	  $dom_score += $m_score;
	  $id_cnt++;
	  if ($$ds_ix == $left_active_end) {
	    $left_active_end = pop_annot_match(\@active_dom_list, $qix, $six, $$ds_ix, $id_cnt, $dom_score);
	    $dom_score = $id_cnt = 0;
	  }
	}

	$qix++;
	$six++;
	$gap0 = $gap1 = 0;
      }
    }
    else {
      ($seq0, $seq1) = split(//,$btop);

#      print "$qix:$six : $btop\n";

      if ($btop=~ m/\-/) {
	if ($seq0 eq '-') {  # gap in seq0
	  if ($gap0) {
	    $m_score = $g_ext;
	  }
	  else {
	    $m_score = $g_open+$g_ext;
	    $gap0 = 1;
	  }

	  $score += $m_score;

	  if ($target) {	# subject domains
	    if ($dom_ix < $dom_nx && $$ds_ix == $dom_r->{d_pos}) {
	      push @aligned_domains, $dom_r;
	      $left_active_end = push_annot_match(\@active_dom_list, $dom_r, $qix, $six, $id_cnt, $dom_score);
	      $dom_ix++;
	      $dom_r = $domain_r->[$dom_ix];
	      ($dom_score, $id_cnt) = (0,0);
	    }
	    if (@active_dom_list) {
	      $dom_score += $m_score;
	      if ($$ds_ix == $left_active_end) {
		$left_active_end = pop_annot_match(\@active_dom_list, $qix, $six, $$ds_ix, $id_cnt, $dom_score);
		$dom_score = $id_cnt = 0;
	      }
	    }
	  }
	  $six++;
	}
	else {  # gap in seq1
	  if ($gap1) {
	    $m_score = $g_ext;
	  }
	  else {
	    $m_score = $g_open+$g_ext;
	    $gap1 = 1;
	  }
	  $score += $m_score;

	  unless ($target) {	# query domains
	    if ($dom_ix < $dom_nx && $$ds_ix == $dom_r->{d_pos}) {
	      push @aligned_domains, $dom_r;
	      $left_active_end = push_annot_match(\@active_dom_list, $dom_r, $qix, $six, $id_cnt, $dom_score);
	      $dom_ix++;
	      $dom_r = $domain_r->[$dom_ix];
	      ($dom_score, $id_cnt) = (0,0);
	    }

	    if (@active_dom_list) {
	      $dom_score += $m_score;
	      if ($$ds_ix == $left_active_end) {
		$left_active_end = pop_annot_match(\@active_dom_list, $qix, $six, $$ds_ix, $id_cnt, $dom_score);
		$dom_score = $id_cnt = 0;
	      }
	    }
	  }
	  $qix++;
	}
      }
      else {	# mismatch
	my ($seq0_map, $seq1_map) = ($aa_map{$seq0},$aa_map{$seq1});
	$seq0_map = $aa_map{'X'} unless defined($seq0_map);
	$seq1_map = $aa_map{'X'} unless defined($seq1_map);

	$m_score = $matrix_2d->[$seq0_map][$seq1_map];
	$score += $m_score;
	if ($dom_ix < $dom_nx && $$ds_ix == $dom_r->{d_pos}) {
	  push @aligned_domains, $dom_r;
	  $left_active_end = push_annot_match(\@active_dom_list, $dom_r, $qix, $six, $id_cnt, $dom_score);
	  $dom_ix++;
	  $dom_r = $domain_r->[$dom_ix];
	  ($dom_score, $id_cnt) = (0,0);
	}

	if (@active_dom_list) {
	  $dom_score += $m_score;
	  if ($$ds_ix == $left_active_end) {
	    $left_active_end = pop_annot_match(\@active_dom_list, $qix, $six, $$ds_ix, $id_cnt, $dom_score);
	    $dom_score = $id_cnt = 0;
	  }
	}
	$qix++;
	$six++;
	$gap0 = $gap1 = 0;
      }
    }
#    print join(":",($qix, $six, $score)),"\n";
  }

  # all done, finish any domain stuff
  if (@active_dom_list) {
    last_annot_match(\@active_dom_list, $hit_r->{q_end}, $hit_r->{s_end}, $id_cnt, $dom_score);
  }

  return ($score, \@aligned_domains);
}

################################################################
# sub_alignment_pos
# input: $hit_r->{BTOP} : an encoded alignment;
#        $domain_r : domain boundaries in query (target=0) or subject (target=1)
#        $target : 0=query, 1=target
#
# updates $domain_r in place with new values:
# domain_r->[]->{ident} (as fraction identical),
#             ->{sa_start, sa_end} domain boundaries in subject
#
sub sub_alignment_pos {
  my ($hit_r, $domain_r, $target) = @_;

  return (0, $domain_r) unless ($domain_r && scalar(@$domain_r));

  my $btop_enc_r = decode_btop($hit_r->{BTOP});

  my ($gap0, $gap1) = (0,0);

  my @active_dom_list = ();
  my @aligned_domains = ();

  my $left_active_end = $domain_r->[-1]->{d_end}+1;	# as far right as possible
  my ($q_start, $s_start, $h_start) = @{$hit_r}{qw(q_start s_start s_start)};
  my ($qix, $six)  = ($q_start, $s_start); # $qix now starts from 1, like $ssix;

  my $ds_ix = \$six;	# use to track the subject position
  # reverse coordinate names if $target==0
  unless ($target) {
    $ds_ix = \$qix;	# track query position
    $h_start = $hit_r->{q_start};
  }

  my ($score, $m_score) = 0;
  my ($seq0, $seq1) = ("","");

  # find the first overlapping domain
  my ($dom_ix, $dom_nx) = (0,scalar(@$domain_r));
  my $dom_r = $domain_r->[0];

  # skip over domains that do not overlap alignment
  # capture first domain that alignment overlaps
  for ($dom_ix=0; $dom_ix < $dom_nx; $dom_ix++) {
    if ($domain_r->[$dom_ix]->{d_end} >= $h_start) {  # if {d_end} < $_start, cannot overlap
      $dom_r = $domain_r->[$dom_ix];
      if ($dom_r->{d_pos} <= $h_start) {  # {d_pos} is less, {d_end} is greater, overlap
	push @aligned_domains, $dom_r;
	$left_active_end = push_annot_match(\@active_dom_list, $dom_r, $q_start, $s_start, 0, 0);
      }
      else { last; }
    }
  }

  my ($dom_score, $id_cnt) = (0,0);

  for my $btop (@{$btop_enc_r}) {

    if ($btop =~ m/^\d+$/) {  # matching query sequence, add it up
      for (my $i=0; $i < $btop; $i++) {  # $i is used to count through BTOP, not to index anything.

	if ($dom_ix < $dom_nx && $$ds_ix == $dom_r->{d_pos}) {
	  push @aligned_domains, $dom_r;
	  $left_active_end = push_annot_match(\@active_dom_list, $dom_r, $qix, $$ds_ix, $id_cnt, $dom_score);
	  $dom_ix++;
	  $dom_r = $domain_r->[$dom_ix];
	  ($dom_score, $id_cnt) = (0,0);
	}
	if (@active_dom_list) {
	  $id_cnt++;
	  if ($$ds_ix == $left_active_end) {
	    $left_active_end = pop_annot_match(\@active_dom_list, $qix, $six, $$ds_ix, $id_cnt, $dom_score);
	    $dom_score = $id_cnt = 0;
	  }
	}

	$qix++;
	$six++;
	$gap0 = $gap1 = 0;
      }
    }
    else {
      ($seq0, $seq1) = split(//,$btop);

#      print "$qix:$six : $btop\n";

      if ($btop=~ m/\-/) {
	if ($seq0 eq '-') {  # gap in seq0

	  if ($target) {	# subject domains
	    if ($dom_ix < $dom_nx && $$ds_ix == $dom_r->{d_pos}) {
	      push @aligned_domains, $dom_r;
	      $left_active_end = push_annot_match(\@active_dom_list, $dom_r, $qix, $$ds_ix, $id_cnt, $dom_score);
	      $dom_ix++;
	      $dom_r = $domain_r->[$dom_ix];
	      ($dom_score, $id_cnt) = (0,0);
	    }
	    if (@active_dom_list) {
	      if ($dom_ix < $dom_nx && $$ds_ix == $left_active_end) {
		$left_active_end = pop_annot_match(\@active_dom_list, $qix, $six, $$ds_ix, $id_cnt, $dom_score);
		$dom_score = $id_cnt = 0;
	      }
	    }
	  }
	  $six++;
	}
	else {  # gap in seq1

	  unless ($target) {	# query domains
	    if ($dom_ix < $dom_nx && $$ds_ix == $dom_r->{d_pos}) {
	      push @aligned_domains, $dom_r;
	      $left_active_end = push_annot_match(\@active_dom_list, $dom_r, $qix, $$ds_ix, $id_cnt, $dom_score);
	      $dom_ix++;
	      $dom_r = $domain_r->[$dom_ix];
	      ($dom_score, $id_cnt) = (0,0);
	    }
	    if (@active_dom_list) {
	      $dom_score += $m_score;
	      if ($dom_ix < $dom_nx && $$ds_ix == $left_active_end) {
		$left_active_end = pop_annot_match(\@active_dom_list, $qix, $six, $$ds_ix, $id_cnt, $dom_score);
		$dom_score = $id_cnt = 0;
	      }
	    }
	  }
	  $qix++;
	}
      }
      else {	# mismatch
	my ($seq0_map, $seq1_map) = ($aa_map{$seq0},$aa_map{$seq1});
	if ($dom_ix < $dom_nx && $$ds_ix == $dom_r->{d_pos}) {
	  push @aligned_domains, $dom_r;
	  $left_active_end = push_annot_match(\@active_dom_list, $dom_r, $qix, $$ds_ix, $id_cnt, $dom_score);
	  $dom_ix++;
	  $dom_r = $domain_r->[$dom_ix];
	  ($dom_score, $id_cnt) = (0,0);
	}
	if (@active_dom_list) {
	  if ($$ds_ix == $left_active_end) {
	    $left_active_end = pop_annot_match(\@active_dom_list, $qix, $six, $$ds_ix, $id_cnt, $dom_score);
	    $dom_score = $id_cnt = 0;
	  }
	}
	$qix++;
	$six++;
	$gap0 = $gap1 = 0;
      }
    }
#    print join(":",($qix, $six, $score)),"\n";
  }

  # all done, finish any domain stuff
  if (@active_dom_list) {
    last_annot_match(\@active_dom_list, $hit_r->{q_end}, $hit_r->{s_end}, $id_cnt, $dom_score);
  }

  return ($score, \@aligned_domains);
}

################
# push_annot_match - adds domain to set of @$active_doms_r,
#		     update ->{score}, ->{ident} for existing @$active_doms_r
#		     initialize ->{score}, ->{ident} to zero for new domain
#		     insert (splice) new domain in list ordered left-to-right by ->{d_end}
#                    returns current left-most {d_end} boundary
#
sub push_annot_match {
  my ($active_doms_r, $dom_r, $q_pos, $s_pos, $c_ident, $c_score) = @_;

  $dom_r->{ident} = 0;
  $dom_r->{score} = 0;
  $dom_r->{qa_start} = $dom_r->{qa_pos} = $q_pos;
  $dom_r->{sa_start} = $dom_r->{sa_pos} = $s_pos;

  # no previous domains, just initialize
  unless (scalar(@$active_doms_r)) {
    push @$active_doms_r, $dom_r;
    return $dom_r->{d_end};
  }

  # some previous domains, update score, identity for domains in list
  # also find insertion point
  my $nx = scalar(@$active_doms_r);
  my $min_ix = $nx;
  for (my $ix=0; $ix < $nx; $ix++) {
    $active_doms_r->[$ix]->{ident} += $c_ident;
    $active_doms_r->[$ix]->{score} += $c_score;
    if ($dom_r->{d_end} < $active_doms_r->[$ix]->{d_end}) {
      $min_ix = $ix;
    }
  }

  # now have location for insert
  splice(@$active_doms_r, $min_ix, 0, $dom_r);
  return $active_doms_r->[0]->{d_end};
}

################
# pop_annot_match - update domains in @$active_doms_r 
#		    update: ->{ident}, ->{score}
#		    add: ->{qa_end},->{sa_end}
#                   remove all domains that end at $s_ix and convert {ident} count to fraction
#                   return left-most right boundary

sub pop_annot_match {
  my ($active_doms_r, $q_pos, $s_pos, $d_pos, $c_ident, $c_score) = @_;

  my $nx = scalar(@$active_doms_r);

  # we know the left most (first) domain matches,
  my $pop_count = 0;
  for my $cur_r (@$active_doms_r) {
    $cur_r->{ident} += $c_ident;
    $cur_r->{score} += $c_score;
    $pop_count++ if ($cur_r->{d_end} == $d_pos);
  }

  while ($pop_count-- > 0) {
    my $cur_r = shift @$active_doms_r;
    # convert identity count to identity fraction
    $cur_r->{percid} = $cur_r->{ident}/($cur_r->{d_end} - $cur_r->{d_pos}+1);
    $cur_r->{qa_end} = $cur_r->{qa_pos} = $q_pos;
    $cur_r->{sa_end} = $cur_r->{sa_pos} = $s_pos;
  }

  if (scalar(@$active_doms_r)) {
    my $leftmost_end = $active_doms_r->[0]->{d_end};
    for (my $lix = 1; $lix < scalar(@$active_doms_r); $lix++) {
      if ($active_doms_r->[$lix]->{d_end} < $leftmost_end) {
	$leftmost_end = $active_doms_r->[$lix]->{d_end};
      }
    }
    return $leftmost_end;
  }
  else {
    return -1;
  }
}

sub last_annot_match {
  my ($active_doms_r, $q_pos, $s_pos, $c_ident, $c_score) = @_;

  my $nx = scalar(@$active_doms_r);

  # we know the left most (first) domain matches,
  my $pop_count = 0;
  for my $cur_r (@$active_doms_r) {
    $cur_r->{ident} += $c_ident;
    $cur_r->{score} += $c_score;
    $cur_r->{percid} = $cur_r->{ident}/($cur_r->{d_end} - $cur_r->{d_pos}+1);
    $cur_r->{qa_end} = $cur_r->{qa_pos} = $q_pos;
    $cur_r->{sa_end} = $cur_r->{sa_pos} = $s_pos;
  }

  $active_doms_r = [];
}

# given: (1) a query sequence; (2) an encoded alignment; (3) a scoring matrix
# report matches/mismatches on annotated sites
# updates $site_r->[]->{q_coord, s_coord}
#                    ->{q_res, s_res}

sub site_align {
  my ($query_r, $hit_r, $matrix_2d, $site_r, $target) = @_;

  return [] unless ($site_r && scalar(@$site_r));

  my @aligned_sites = ();

  my $btop_enc_r = decode_btop($hit_r->{BTOP});

  my ($q_start, $q_end, $s_start, $s_end) = @{$hit_r}{qw(q_start q_end s_start s_end)};
  my ($qix, $six)  = ($q_start, $s_start); # $qix, $six 1-based
  my $ds_ix = \$six;	# use to track the subject position

  unless ($target) {
    ($q_start, $q_end, $s_start, $s_end) = @{$hit_r}{qw(s_start s_end q_start q_end)};
    $ds_ix = \$qix;	# track query position
  }

  my ($seq0, $seq1) = ("","");

  # find the first overlapping domain

  my ($site_ix, $site_nx) = (0,scalar(@$site_r));
  my $s_r = $site_r->[0];

  # skip over sites that do not overlap alignment
  for ($site_ix=0; $site_ix < $site_nx; $site_ix++) {
    if ($site_r->[$site_ix]->{d_pos} >= $s_start) {  # find the first site inside alignment
      $s_r = $site_r->[$site_ix];
      last;
    }
  }

  return [] unless $site_ix < $site_nx;

  for my $btop (@{$btop_enc_r}) {
    last if ($site_ix >= $site_nx);
    if ($btop =~ m/^\d+$/) {  # matching query sequence, check for sites within current region
      my $bt_end = $$ds_ix + $btop - 1;
      if ($bt_end < $s_r->{d_pos}) {  # no site in identical region
	$qix += $btop;
	$six += $btop;
      }
      else {      # yes site in region, jump to it
	my $c_pos;
	while ($site_ix < $site_nx && $s_r->{d_pos} <= $bt_end) {
	  $c_pos = $$ds_ix;	# must be inside loop because $ds_ix points to $qix or $six
	  $qix += $s_r->{d_pos} - $c_pos;	# jump forward to site
	  $six += $s_r->{d_pos} - $c_pos;	# jump forward to site
	  $seq0 = $query_r->[$qix];

	  @{$s_r}{qw(annot_ix qa_pos sa_pos q_res s_res m_symb d_end)} = ($site_ix, $qix, $six, $seq0, $seq0, match_symb($seq0, $seq0, $matrix_2d));
	  push @aligned_sites, $s_r;
	  $site_ix++;
	  $s_r=$site_r->[$site_ix];
	}
	# past the last site annotation, but not done with $btop;
	$c_pos = ($bt_end - $$ds_ix + 1);
	$qix += $c_pos;
	$six += $c_pos;
      }
    }
    else {	# sequence does not match -- must check each position
      ($seq0, $seq1) = split(//,$btop);
      if ($btop =~ m/\-/) {
	if ($seq0 eq '-') {
	  if ($target) {
	    while ($site_ix < $site_nx && $s_r->{d_pos} == $six) {
	      @{$s_r}{qw(annot_ix qa_pos sa_pos q_res s_res m_symb)} = ($site_ix, $qix, $six, $seq0, $seq1, match_symb($seq0, $seq1, $matrix_2d));
	      push @aligned_sites, $s_r;
	      $site_ix++;
	      $s_r=$site_r->[$site_ix];
	    }
	  }
	  $six++;
	}
	else {  # gap in seq1, cannot match domain
	  unless ($target) {
	    while ($site_ix < $site_nx && $s_r->{d_pos} == $qix) {
	      @{$s_r}{qw(annot_ix qa_pos sa_pos q_res s_res m_symb)} = ($site_ix, $qix, $six, $seq0, $seq1, match_symb($seq0, $seq1, $matrix_2d));
	      push @aligned_sites, $s_r;
	      $site_ix++;
	      $s_r=$site_r->[$site_ix];
	    }
	  }
	  $qix++;
	}
      }
      else {	# mismatch; $btop string is twice length of covered region
	while ($s_r->{d_pos} == $$ds_ix && $site_ix < $site_nx ) {
	  @{$s_r}{qw(annot_ix qa_pos sa_pos q_res s_res m_symb)} = ($site_ix, $qix, $six, $seq0, $seq1, match_symb($seq0, $seq1, $matrix_2d));
	  push @aligned_sites, $s_r;
	  $site_ix++;  $s_r=$site_r->[$site_ix];
	}
	$qix++;
	$six++;
      }
    }
  }

  return (\@aligned_sites);
}

sub match_symb {
  my ($seq0, $seq1, $matrix_2d) = @_;

  if (uc($seq0) eq uc($seq1)) {
    return "=";
  }
  else {
    my $seq0_map = $aa_map{$seq0};
    $seq0_map = $aa_map{'X'} unless defined($seq0_map);

    my $seq1_map = $aa_map{$seq1};
    $seq1_map = $aa_map{'X'} unless defined($seq1_map);

    my $m_score = $matrix_2d->[$seq0_map][$seq1_map];

    if ($m_score < 0) {return "<";}
    elsif ($m_score > 0) {return ">";}
    else {return "z";}
  }
}



# merge up to four lists of annotations into a single list, and return
# a reference to the list
# input: $hit references, possibly with {aligned_domains_r}, {aligned_sites_r}
#                                       {q_aligned_domains_r}, {q_aligned_sites_r}
#
sub merge_annots {
  my ($hit_r) = @_;

  my @merged_array = ();

  # merge the sites arrays first, so that conserved annotated sites are juxtaposed

  my ($qs_ix, $ss_ix, $qs_nx, $ss_nx) = (0,0,0,0);

  $ss_nx = scalar(@{$hit_r->{aligned_sites_r}}) if (exists($hit_r->{aligned_sites_r}));
  $qs_nx = scalar(@{$hit_r->{q_aligned_sites_r}}) if (exists($hit_r->{q_aligned_sites_r}));

  if ($ss_nx && $qs_nx) {  # have sites on both sequences
    # find out how many positions match between {q_aligned_sites_r} and {aligned_sites_r}

    my @uniq_sites = ();

    for my $qs_ref (@{$hit_r->{q_aligned_sites_r}}) {
      $qs_ref->{merged} = 0;
      for my $ss_ref ( @{$hit_r->{aligned_sites_r}} ) {
	next if ($ss_ref->{qa_pos} < $qs_ref->{qa_pos});
	last if ($ss_ref->{qa_pos} > $qs_ref->{qa_pos});
	if ($qs_ref->{qa_pos} == $ss_ref->{qa_pos} && $qs_ref->{type} eq $ss_ref->{type}) {
	  $qs_ref->{merged} = $ss_ref->{merged} = 1;
	  $qs_ref->{target} = $ss_ref->{target} = 2;
	  # save match
	  push @uniq_sites, $qs_ref;
	}
      }
    }

    # save merged sites
    push @merged_array, @uniq_sites;

    # save unmerged subject
    @uniq_sites = ();
    for my $ss_ref ( @{$hit_r->{aligned_sites_r}} ) {
      push @uniq_sites, $ss_ref if (!defined($ss_ref->{merged}) || $ss_ref->{merged} == 0);
    }
    push @merged_array, @uniq_sites;

    # save unmerged query
    @uniq_sites = ();
    for my $qs_ref ( @{$hit_r->{aligned_sites_r}} ) {
      push @uniq_sites, $qs_ref if (!defined($qs_ref->{merged}) || $qs_ref->{merged} == 0);
    }
    push @merged_array, @uniq_sites;
  }
  elsif ($ss_nx) {
    push @merged_array, @{$hit_r->{aligned_sites_r}};
  }
  elsif ($qs_nx) {
    push @merged_array, @{$hit_r->{q_aligned_sites_r}};
  }

#  for my $ann_r ( @merged_array) {
#    unless ($ann_r->{qa_pos}) {
#      print STDERR "missing qa_pos:",join(":",@{$ann_r}{qw(q_seqid s_seqid)}),"\n";
#    }
#  }

  @merged_array = sort { $a->{qa_pos} <=> $b->{qa_pos} } @merged_array;


  push @merged_array, @{$hit_r->{aligned_domains_r}} if (exists($hit_r->{aligned_domains_r}));
  push @merged_array, @{$hit_r->{q_aligned_domains_r}} if (exists($hit_r->{q_aligned_domains_r}));

  @merged_array = sort { $a->{qa_pos} <=> $b->{qa_pos} } @merged_array;

  return \@merged_array;
}

####
# print raw domain info:
# |DX:%d-%d;C=dom_info|XD:%d-%d:C=dom_info
#
sub format_dom_info {
  my ($q_dom_r, $dom_r) = @_;

  my $dom_str = "";
  for my $dom ( @$q_dom_r ) {
    $dom_str .= sprintf("|DX:%d-%d;C=%s",@{$dom}{qw(d_pos d_end descr)});
  }
  for my $dom ( @$dom_r ) {
    $dom_str .= sprintf("|XD:%d-%d;C=%s",@{$dom}{qw(d_pos d_end descr)});
  }

  return $dom_str;
}

# merged annot output formatter
sub format_annot_info {
  my ($hit_r, $annot_list_r) = @_;

  my $raw_score = 0;

  if  ($hit_r->{raw_score} ) {
    $raw_score = $hit_r->{raw_score};
  }
  else {
#    warn "no raw_score at: ".$hit_r->{s_seqid}."\n";
    $raw_score = $hit_r->{score};
  }

  my $score_scale = $hit_r->{score}/$raw_score;

  my $annot_str = "";

  # two types of annotations, domains and sites.

  for my $annot_r ( @$annot_list_r ) {

    if ($annot_r->{type} eq '-') { # domain with scores
      my $fsub_score = $annot_r->{score}/$raw_score;

      my ($ns_score, $s_bit) = (int($annot_r->{score} * $score_scale + 0.5),
				int($hit_r->{bits} * $fsub_score + 0.5),
			       );
      my $qval = 0.0;
      if ($hit_r->{evalue} == 0.0) {
	if ($s_bit > 50) {
	  $qval = 3000.0
	}
	else {
	  $qval = -10.0 * (log(400.0 * 400.) + $s_bit)/log(10.0);
	}
      } else {
	$qval = -10.0*log($hit_r->{evalue})*$fsub_score/(log(10.0))
      }

      $qval = 0 if $qval < 0;

      $annot_str .= join(";",(sprintf("|%s:%d-%d:%d-%d:s=%d",
				      $annot_r->{target} ? "XR" : "RX",
				      $annot_r->{qa_start},$annot_r->{qa_end},
				      $annot_r->{sa_start},$annot_r->{sa_end},$ns_score),
			      sprintf("b=%.1f",$s_bit),
			      sprintf("I=%.3f",$annot_r->{percid}),
			      sprintf("Q=%.1f",$qval),"C=".$annot_r->{descr}));
    }
    else {	# site annotation
      my $ann_type = $annot_r->{type};
      my $site_str = "|".$ann_type . "X";
      if ($annot_r->{target} == 1) {
	$site_str = "|X".$ann_type;
      }
      elsif ($annot_r->{target} == 2) {
	$site_str = "|$ann_type$ann_type";
      }

      $annot_str .= "$site_str:" . sprintf("%d%s%s%d%s",
					   $annot_r->{qa_pos}, $annot_r->{q_res}, $annot_r->{m_symb}, $annot_r->{sa_pos}, $annot_r->{s_res});

    }
  }
  return $annot_str;
}

__END__

=pod

=head1 NAME

annot_blast_btop2.pl

=head1 SYNOPSIS

 annot_blast_btop2 --ann_script ann_pfam_www_e.pl [--query_file query.fasta] --out_fields "q_seqid s_seqid percid evalue" blast_tabular_file

=head1 OPTIONS

 -h	short help
 --help include description

 --ann_script -- annotation script returning site/domain locations for subject sequences
              -- same as --script

 --have_qslen -- use a blast tabular format that includes the query and subject sequence lengths:
	      -- q_seqid q_len s_seqid s_len ...

 --q_ann_script -- annotation script for query sequences
                -- same as --q_script

 --query_file -- fasta query sequence
              -- same as --query, --query_lib
                 (can contain multiple sequences for multi-sequence search)

 --out_fields -- blast tabular fields shown before domain information

 --raw_score -- add the raw_score used to normalized domain scores to
                tabular output (raw_scores are only calculated for domains)

=head1 DESCRIPTION

C<annot_blast_btop2.pl> runs the script specified by
C<--ann_script/--q_ann_script> to annotate functional sites domain
content of the sequences specified by the subject/query seqid field of
blast tabular format (-outfmt 6 or 7) or FASTA blast tabular format
(-m 8).  The C<--ann_script/--q_ann_script> script produces domain
boundary coordinates, which are mapped to the alignment.  For searches
against SwissProt sequences, C<--ann_script ann_feats_up_www2.pl> will
acquire features and domains from Uniprot.  C<--ann_script
ann_pfam_www.pl --neg> will get domain information from Pfam, and
score non-domain (NODOM) regions.

The tab file is read and parsed, and then the subject/query seqid is used to
capture domain locations in the subject/query sequence.  If the domains
overlap the aligned region, the domain names are appended to the
intput.

If a C<--query_file> is specified and two additional fields, C<score>
and C<btop> are present, C<annot_blast_btop2.pl> calculates
sub-alignment scores, including fraction identity, bit score, and
Q-value (-log10(E-value)), partitioning the alignment score, identity,
and bit score across the overlapping domains.

The C<--out_fields> specifies the blast tabular fields that can be
returned.  By default, C<q_seqid s_seqid percid alen mismatch gopen
q_start q_end s_start s_end evalue bits> (but not C<score> and
C<BTOP>) are shown.

Currently, this program is fully functional only for blastp (or
blastn) searches.  For translated searches (blastx) domain content,
location and identity is provided, but not bit-scores or Q-values.

=head1 AUTHOR

William R. Pearson, wrp@virginia.edu

=cut
