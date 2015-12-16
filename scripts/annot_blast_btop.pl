#!/usr/bin/perl -w

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
# annot_blast_btop.pl --query query.file --ann_script ann_pfam_www.pl blast_tab_btop_file
################################################################
# annot_blast_btop.pl associates domain annotation information and
# subalignment scores with a blast tabular (-outfmt 6 or -outfmt 7)
# file that contains the raw score and the BTOP alignment encoding
# This file can be generated from "blastp/n" or "blast_formatter"
# using the command:
#   blast_formatter -archive blast_output.asn -outfmt '7 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore score btop'  > blast_output.tab_annot
#
# If the BTOP field or query_file is not available, the script
# produces domain content without sub-alignment scores.
################################################################

## 16-Nov-2015
# modify to allow multi-query blast searches

use strict;
use IPC::Open2;
use Pod::Usage;
use Getopt::Long;
# use Data::Dumper;

# read lines of the form:
# gi|121694|sp|P20432.1|GSTT1_DROME	gi|121694|sp|P20432|GSTT1_DROME	100.00	209	0	0	1	209	1	209	6e-156	433	1113	209
# gi|121694|sp|P20432.1|GSTT1_DROME	gi|1170090|sp|P04907|GSTF3_MAIZE	26.77	198	123	7	4	185	6	197	2e-08	51.2	121	FL1YG ... 1NKRA1YW1
# gi|121694|sp|P20432.1|GSTT1_DROME	gi|81174731|sp|P0ACA5|SSPA_ECO57	39.66	58	32	2	43	100	49	103	8e-06	43.9	102	EDFLLI ... V-I-NEQS3FM
# gi|121694|sp|P20432.1|GSTT1_DROME	gi|121695|sp|P12653|GSTF1_MAIZE	27.62	181	107	7	32	203	34	199	9e-05	40.8	94	LI1LF ... N-1AS1CLLM1

# and report the domain content ala -m 8CC

my ($matrix, $ann_script, $shelp, $help) = ("BLOSUM62", "ann_feats_up_www2.pl --neg --no-feats", 0, 0);
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
    "ann_script:s" => \$ann_script,
    "query:s" => \$query_lib_name,
    "query_file:s" => \$query_lib_name,
    "query_lib:s" => \$query_lib_name,
    "out_fields:s" => \$out_field_str,
    "script" => \$ann_script,
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

# the fields that are displayed are listed here.  By default, all fields except score and BTOP are displayed.
my @out_tab_fields = @tab_fields[0 .. $#tab_fields-1];
push @out_tab_fields, "raw_score";

if ($out_field_str) {
  @out_tab_fields = split(/\s+/,$out_field_str);
}

my @header_lines = ();

# need outer loop to enable multiple queries
while (1) {

  my $next_line = "";
  my $have_data = 0;

  my @hit_list = ();
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
    @hit_data{@tab_fields} = split(/\t/,$line);

    push @hit_list, \%hit_data;
  }

  # get the current query sequence
  if ($ann_script && -x (split(/\s+/,$ann_script))[0]) {
    # get the domains for each s_seqid using --ann_script
    #
    local (*Reader, *Writer);
    my $pid = open2(\*Reader, \*Writer, $ann_script);
    for my $hit (@hit_list) {
      print Writer $hit->{s_seqid},"\n";
    }
    close(Writer);

    my $current_domain = "";
    my $hit_ix = 0;
    my @hit_domains = ();

    while (my $line = <Reader>) {
      chomp $line;
      if ($line =~ m/^>/) {
	if ($current_domain) {
	  if ($hit_list[$hit_ix]{s_seqid} eq $current_domain) {
	    $hit_list[$hit_ix]{domains} = [ @hit_domains ];
	    $hit_ix++;
	  } else {
	    warn "phase error: $current_domain != $hit_list[$hit_ix]{s_seqid}";
	  }
	}
	@hit_domains = ();
	$current_domain = $line;
	$current_domain =~ s/^>//;
      } else {
	next if $line=~ m/^=/;
	my %dom_info = ();
	@dom_info{qw(sd_start dash sd_end descr)} = split(/\t/,$line);
	next unless $dom_info{dash} eq '-';
	$dom_info{descr} =~ s/ :(\d+)$/~$1/;
	delete($dom_info{dash});
	push @hit_domains, \%dom_info;
      }
    }
    close(Reader);
    waitpid($pid, 0);

    if (@hit_domains) {
      $hit_list[$hit_ix]{domains} = [ @hit_domains ];
    }
  }

  for my $line (@header_lines) {
    print $line;
  }
  @header_lines = ($next_line);

  # now get query sequence if available

  for my $hit (@hit_list) {
    my @list_covered = ();

    # If I have an encoded aligment and a query sequence, I can calculate sub-alignment scores
    if (defined($hit->{BTOP}) && 
	$query_lib_r && $query_lib_r->{$hit->{q_seqid}}
       && $hit->{domains}) {
      ($hit->{raw_score}, $hit->{aligned_domains_r}) = 
	sub_alignment_score($query_lib_r->{$hit->{q_seqid}},
			    $hit, \@blosum62, \@blosum62_diag, $hit->{domains});
    } else {		   # no alignment info, just check for overlap
      $hit->{raw_score} = 0;
      for my $dom_r (@{$hit->{domains}}) {
	next if $dom_r->{sd_end} < $hit->{s_start}; # before start
	last if $dom_r->{sd_start} > $hit->{s_end}; # after end

	if ($dom_r->{sd_start} <= $hit->{s_end} && $dom_r->{sd_end} >= $hit->{s_start}) {
	  push @list_covered, $dom_r->{descr};
	}
      }
    }

    ################
    ## final output display

    print join("\t",@{$hit}{@out_tab_fields}); # show fields from original blast tabular file

    if (defined($hit->{aligned_domains_r})) { # show subalignment scores if available
      print "\t";
      for my $dom_r ( @{$hit->{aligned_domains_r}} ) {
	print format_dom_info($hit, $hit->{raw_score}, $dom_r);
      }
    } elsif (@list_covered) {	# otherwise show domain content
      print "\t",join(";",@list_covered);
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

  my @ncbi_blaa = qw(A  R  N  D  C  Q  E  G  H  I  L  K  M  F  P  S  T  W  Y  V  B  Z  X * );

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

# given: (1) a query sequence; (2) an encoded alignment; (3) a scoring matrix
# calculate a score
# updates $domain_r in place with new values:
# domain_r->[]->{ident} (as fraction identical),
#             ->{score} --matrix raw similarity score
#             ->{qa_start,qa_end}  domain boundaries in query
#             ->{sa_start, sa_end} domain boundaries in subject

sub sub_alignment_score {
  my ($query_r, $hit_r, $matrix_2d, $matrix_diag, $domain_r) = @_;

  return (0, $domain_r) unless ($domain_r && scalar(@$domain_r));

  my $btop_enc_r = decode_btop($hit_r->{BTOP});

  my ($gap0, $gap1) = (0,0);

  my @active_dom_list = ();
  my @aligned_domains = ();

  my $left_active_end = $domain_r->[-1]->{sd_end}+1;	# as far right as possible

  my ($q_start, $q_end, $s_start, $s_end) = @{$hit_r}{qw(q_start q_end s_start s_end)};

  my ($qix, $six)  = ($q_start-1, $s_start); # $qix starts from zero, but $six stays 1-based

  my ($score, $m_score) = 0;
  my ($seq0, $seq1) = ("","");

  # find the first overlapping domain

  my ($sdom_ix, $sdom_nx) = (0,scalar(@$domain_r));
  my $dom_r = $domain_r->[0];

  # skip over domains that do not overlap alignment
  # capture first domain that alignment overlaps
  for ($sdom_ix=0; $sdom_ix < $sdom_nx; $sdom_ix++) {
    if ($domain_r->[$sdom_ix]->{sd_end} >= $s_start) {  # if {sd_end} < $_start, cannot overlap
      $dom_r = $domain_r->[$sdom_ix];
      if ($dom_r->{sd_start} <= $s_start) {  # {sd_start} is less, {sd_end} is greater, overlap
	push @aligned_domains, $dom_r;
	$left_active_end = push_annot_match(\@active_dom_list, $dom_r, $q_start, $s_start, 0, 0);
      }
      else { last; }
    }
  }

  my ($s_dom_score, $s_id_cnt) = (0,0);

  for my $btop (@{$btop_enc_r}) {

    if ($btop =~ m/^\d+$/) {  # matching query sequence, add it up
      for (my $i=0; $i < $btop; $i++) {

	my $seq0_map = $aa_map{$query_r->[$qix]};
	$seq0_map = $aa_map{'X'} unless defined($seq0_map);

	$m_score = $matrix_diag->[$seq0_map];
	$score += $m_score;

	if ($sdom_ix < $sdom_nx && $six == $dom_r->{sd_start}) {
	  push @aligned_domains, $dom_r;
	  $left_active_end = push_annot_match(\@active_dom_list, $dom_r, $qix+1, $six, $s_id_cnt, $s_dom_score);
	  $sdom_ix++;
	  $dom_r = $domain_r->[$sdom_ix];
	  ($s_dom_score, $s_id_cnt) = (0,0);
	}
	if (@active_dom_list) {
	  $s_dom_score += $m_score;
	  $s_id_cnt++;
	  if ($six == $left_active_end) {
	    $left_active_end = pop_annot_match(\@active_dom_list, $qix+1, $six, $s_id_cnt, $s_dom_score);
	    $s_dom_score = $s_id_cnt = 0;
	  }
	}

	$qix++;
	$six++;
	$gap0 = $gap1 = 0;
      }
    }
    else {
      ($seq0, $seq1) = split(//,$btop);

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

	  if ($sdom_ix < $sdom_nx && $six == $dom_r->{sd_start}) {
	    push @aligned_domains, $dom_r;
	    $left_active_end = push_annot_match(\@active_dom_list, $dom_r, $qix+1, $six, $s_id_cnt, $s_dom_score);
	    $sdom_ix++;
	    $dom_r = $domain_r->[$sdom_ix];
	    ($s_dom_score, $s_id_cnt) = (0,0);
	  }
	  if (@active_dom_list) {
	    $s_dom_score += $m_score;
	    if ($sdom_ix < $sdom_nx && $six == $left_active_end) {
	      $left_active_end = pop_annot_match(\@active_dom_list, $qix+1, $six, $s_id_cnt, $s_dom_score);
	      $s_dom_score = $s_id_cnt = 0;
	    }
	  }
	  $six++;
	}
	else {  # gap in seq1, cannot match domain
	  if ($gap1) {
	    $m_score = $g_ext;
	  }
	  else {
	    $m_score = $g_open+$g_ext;
	    $gap1 = 1;
	  }
	  $score += $m_score;
	  $qix++;
	}
      }
      else {	# mismatch
	my ($seq0_map, $seq1_map) = ($aa_map{$seq0},$aa_map{$seq1});
	$seq0_map = $aa_map{'X'} unless defined($seq0_map);
	$seq1_map = $aa_map{'X'} unless defined($seq1_map);

	$m_score = $matrix_2d->[$seq0_map][$seq1_map];
	$score += $m_score;
	if ($sdom_ix < $sdom_nx && $six == $dom_r->{sd_start}) {
	  push @aligned_domains, $dom_r;
	  $left_active_end = push_annot_match(\@active_dom_list, $dom_r, $qix+1, $six, $s_id_cnt, $s_dom_score);
	  $sdom_ix++;
	  $dom_r = $domain_r->[$sdom_ix];
	  ($s_dom_score, $s_id_cnt) = (0,0);
	}
	if (@active_dom_list) {
	  $s_dom_score += $m_score;
	  if ($six == $left_active_end) {
	    $left_active_end = pop_annot_match(\@active_dom_list, $qix+1, $six, $s_id_cnt, $s_dom_score);
	    $s_dom_score = $s_id_cnt = 0;
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
    last_annot_match(\@active_dom_list, $q_end, $s_end, $s_id_cnt, $s_dom_score);
  }

  return ($score, \@aligned_domains);
}

################
# push_annot_match - adds domain to set of @$active_doms_r,
#		     update ->{score}, ->{ident} for existing @$active_doms_r
#		     initialize ->{score}, ->{ident} to zero for new domain
#		     insert (splice) new domain in list ordered left-to-right by ->{sd_end}
#                    returns current left-most {sd_end} boundary
#
sub push_annot_match {
  my ($active_doms_r, $dom_r, $q_pos, $s_pos, $c_ident, $c_score) = @_;

  $dom_r->{ident} = 0;
  $dom_r->{score} = 0;
  $dom_r->{qa_start} = $q_pos;
  $dom_r->{sa_start} = $s_pos;

  # no previous domains, just initialize
  unless (scalar(@$active_doms_r)) {
    push @$active_doms_r, $dom_r;
    return $dom_r->{sd_end};
  }

  # some previous domains, update score, identity for domains in list
  # also find insertion point
  my $nx = scalar(@$active_doms_r);
  my $min_ix = $nx;
  for (my $ix=0; $ix < $nx; $ix++) {
    $active_doms_r->[$ix]->{ident} += $c_ident;
    $active_doms_r->[$ix]->{score} += $c_score;
    if ($dom_r->{sd_end} < $active_doms_r->[$ix]->{sd_end}) {
      $min_ix = $ix;
    }
  }

  # now have location for insert
  splice(@$active_doms_r, $min_ix, 0, $dom_r);
  return $active_doms_r->[0]->{sd_end};
}

################
# pop_annot_match - update domains in @$active_doms_r 
#		    update: ->{ident}, ->{score}
#		    add: ->{qa_end},->{sa_end}
#                   remove all domains that end at $s_ix and convert {ident} count to fraction
#                   return left-most right boundary

sub pop_annot_match {
  my ($active_doms_r, $q_pos, $s_pos, $c_ident, $c_score) = @_;

  my $nx = scalar(@$active_doms_r);

  # we know the left most (first) domain matches,
  my $pop_count = 0;
  for my $cur_r (@$active_doms_r) {
    $cur_r->{ident} += $c_ident;
    $cur_r->{score} += $c_score;
    $pop_count++ if ($cur_r->{sd_end} == $s_pos);
  }

  while ($pop_count-- > 0) {
    my $cur_r = shift @$active_doms_r;
    # convert identity count to identity fraction
    $cur_r->{ident} = $cur_r->{ident}/($cur_r->{sd_end} - $cur_r->{sd_start}+1);
    $cur_r->{qa_end} = $q_pos;
    $cur_r->{sa_end} = $s_pos;
  }
  if (scalar(@$active_doms_r)) {
    return $active_doms_r->[0]->{sd_end};
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
    $cur_r->{ident} = $cur_r->{ident}/($cur_r->{sd_end} - $cur_r->{sd_start}+1);
    $cur_r->{qa_end} = $q_pos;
    $cur_r->{sa_end} = $s_pos;

  }

  $active_doms_r = [];
}

# domain output formatter

sub format_dom_info {
  my ($hit_r, $raw_score, $dom_r) = @_;

  unless ($raw_score) {
    warn "no raw_score at: ".$hit_r->{s_seqid}."\n";
    $raw_score = $hit_r->{score};
  }

  my ($score_scale, $fsub_score) = ($hit_r->{score}/$raw_score, $dom_r->{score}/$raw_score);

  my $qval = 0.0;
  if ($hit_r->{evalue} == 0.0) {
    $qval = 3000.0
  }
  else {
    $qval = -10.0*log($hit_r->{evalue})*$fsub_score/(log(10.0))
  }

  my ($ns_score, $s_bit) = (int($dom_r->{score} * $score_scale+0.5),
			      int($hit_r->{bits} * $fsub_score +0.5),
			     );
  $qval = 0 if $qval < 0;

  #	print join(":",($dom_r->{asd_start},$dom_r->{asd_end},$ns_score, $s_bit, sprintf("%.1f",$qval))),"\n";
  return join(";",(sprintf("|XR:%d-%d:%d-%d:s=%d",
			   $dom_r->{qa_start},$dom_r->{qa_end},
			   $dom_r->{sa_start},$dom_r->{sa_end},$ns_score),
		   sprintf("b=%.1f",$s_bit),
		   sprintf("I=%.3f",$dom_r->{ident}),
		   sprintf("Q=%.1f",$qval),$dom_r->{descr}));
}

__END__

=pod

=head1 NAME

annot_blast_btop.pl

=head1 SYNOPSIS

 annot_blast_btop --ann_script ann_pfam_www_e.pl [--query_file query.fasta] --out_fields "q_seqid s_seqid percid evalue" blast_tabular_file

=head1 OPTIONS

 -h	short help
 --help include description

 --ann_script -- annotation script returning domain locations
 --query_file -- fasta query sequence
 --out_fields -- blast tabular fields shown before domain information

=head1 DESCRIPTION

C<annot_blast_btop.pl> runs the script specified by C<--ann_script> to
annotate the domain content of the sequences specified by the subject
seqid field of blast tabular format (-outfmt 6 or 7) or FASTA blast
tabular format (-m 8).  The C<--ann_script> file is run to produce
domain boundary coordinates; if no C<--ann_script> is provided, domains
are downloaded from UniProt (ann_feats_up_www2.pl).

The tab file is read and parsed, and then the subject seqid is used to
capture domain locations in the subject sequence.  If the domains
overlap the aligned region, the domain names are appended to the
intput.

If a C<--query_file> is specified and two additional fields, C<score>
and C<btop> are present, C<annot_blast_btop.pl> calculates
sub-alignment scores, including fraction identity, bit score, and
Q-value (-log10(E-value)), partitioning the alignment score, identity,
and bit score across the overlapping domains.

The C<--out_fields> specifies the blast tabular fields that can be
returned.  By default, C<q_seqid s_seqid percid alen mismatch gopen
q_start q_end s_start s_end evalue bits> (but not C<score> and
C<BTOP>) are shown.

=head1 AUTHOR

William R. Pearson, wrp@virginia.edu

=cut
