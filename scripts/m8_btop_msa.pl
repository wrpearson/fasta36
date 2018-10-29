#!/usr/bin/env perl

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
# m8_btop_msa.pl --query query.file blast_tab_btop_file
################################################################
# m8_btop_msa.pl takes a query sequence and a blast tabular format
# file with a BTOP field, and constructs a query-driven multiple
# sequence alignment of the subject sequences that can be used as
# input to psiblast with the "--in_msa msa.file" option.
#
# (because BLAST BTOP encoding provides the mismatched residues, the
# library sequences are not required to produce the MSA -- they are
# available in the BTOP string)
# 
# The BTOP alignment encoding file generated from "blastp/n" or
# "blast_formatter" using the command: blast_formatter -archive
# blast_output.asn -outfmt '7 qseqid sseqid pident length mismatch
# gapopen qstart qend sstart send evalue bitscore score btop' >
# blast_output.tab_annot
#
# the raw score shown above is used by the annot_blast_btop2.pl
# program, if present, but is not required by m8_btop_msa.pl
#
################################################################

use warnings;
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

my ($shelp, $help, $evalue) = (0, 0, 0.001);
my ($query_file, $bound_file) = ("","");
my ($out_field_str) = ("");
my $query_lib_r = 0;

GetOptions(
    "query=s" => \$query_file,
    "query_file=s" => \$query_file,
    "evalue=f" => \$evalue,
    "bound_file=s" => \$bound_file,
    "bound=s" => \$bound_file,
    "seqbdr=s" => \$bound_file,
    "h|?" => \$shelp,
    "help" => \$help,
    );

pod2usage(1) if $shelp;
pod2usage(exitstatus => 0, verbose => 2) if $help;
unless (-f STDIN || -p STDIN || @ARGV) {
 pod2usage(1);
}

my @hit_list = ();
my @multi_align = ();
my @multi_names = ();

################
# get query sequence, and insert into MSA
#
my ($query_acc, $query_seq_r, $query_len);
if ($query_file) {
  ($query_acc, $query_seq_r) = parse_query_lib($query_file);
  $query_len = scalar(@$query_seq_r)-1;  # -1 for ' ' 1: offset
}

if (! $query_file || !$query_len) {
  die "query sequence required";
}

push @multi_names, $query_acc;
push @multi_align, btop2alignment($query_seq_r, $query_len, {BTOP=>$query_len, q_start=>1, q_end=>$query_len});
my $max_sseqid_len = length($query_acc);

################
# get sequence boundaries if available
#
my $seq_bound_hr = 0;

if ($bound_file) {
  $seq_bound_hr = parse_bound_file($bound_file)
}

my @tab_fields = qw(q_seqid s_seqid percid alen mismatch gopen q_start q_end s_start s_end evalue bits score BTOP);

while (my $line = <>) {
  if ($line =~ m/^# Fields:/ && $line !~ m/bit score, score, BTOP/)  {
    # raw score missing, edit @tab_fields
    pop @tab_fields;
    pop @tab_fields;
    push @tab_fields, "BTOP";
    next;
  }
  next if ($line =~ m/^#/);
  chomp $line;
  next unless $line;

  my %hit_data = ();

  @hit_data{@tab_fields} = split(/\t/,$line);

  next if ($hit_data{evalue} > $evalue);

#  push @hit_list, \%hit_data;
  if (length($hit_data{s_seqid}) > $max_sseqid_len) {
    $max_sseqid_len = length($hit_data{s_seqid});
  }

  if ($bound_file) {
    if (defined($seq_bound_hr->{$hit_data{subj_acc}})) {
      push @multi_names, $hit_data{s_seqid};
      push @multi_align, bound_btop2alignment($query_seq_r, $query_len, \%hit_data, @{$seq_bound_hr->{$hit_data{subj_acc}}}{qw(start end)});
    }
  }
  else {  # no sequence boundaries
    push @multi_names, $hit_data{s_seqid};
    push @multi_align, btop2alignment($query_seq_r, $query_len, \%hit_data);
  }
}

# final MSA output
$max_sseqid_len += 4;

print "SSEARCHm8 multiple sequence alignment\n\n\n";

my $i_pos = 0;
for (my $j = 0; $j < $query_len/60; $j++) {
  my $i_end = $i_pos + 59;
  if ($i_end > $query_len) {$i_end = $query_len-1;}
  for (my $n = 0; $n < scalar(@multi_names); $n++) {
    printf("%-".$max_sseqid_len."s %s\n",$multi_names[$n],join("",@{$multi_align[$n]}[$i_pos .. $i_end]));
  }
  $i_pos += 60;
  print "\n\n";
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


sub btop2alignment {
  my ($query_seq_r, $query_len, $hit_data_hr) = @_;

  # $query_seq_r is 1: based
  my @alignment = ();

  # make a local copy
  # my @query_seq = @{$query_seq_r};

  # the left unaligned region gets " ";
  for (my $i=1; $i < $hit_data_hr->{q_start}; $i++) {
    push @alignment, "-";
  }

  my $btop_align_r = decode_btop($hit_data_hr->{BTOP});

  my ($seq0, $seq1) = ("","");
  my $qix = $hit_data_hr->{q_start};
  for my $btop (@{$btop_align_r}) {
    if ($btop =~ m/^\d+$/) {  # matching query sequence, add it up
      for (my $i=0; $i < $btop; $i++) {
	push @alignment, $query_seq_r->[$qix++];
      }
    }
    else {  # could be: TS/-S/T-
      ($seq0, $seq1) = split(//,$btop);
      if ($seq0 ne '-') {
	push @alignment, $seq1;
	$qix++;
      }
    }
  }
  # all done with alignment, double check that $qix = $hit_data_hr->{q_end}
  unless ($qix == $hit_data_hr->{q_end}+1) {
    warn "$qix != ".$hit_data_hr->{q_end}+1;
  }

  for (my $i = $hit_data_hr->{q_end}+1; $i <= $query_len; $i++) {
    push @alignment, "-";
  }

  return \@alignment;
}

sub bound_btop2alignment {
  my ($query_seq_r, $query_len, $hit_data_hr, $sb_start, $sb_end) = @_;

  # $query_seq_r is 1: based
  my @alignment = ();

  # the left unaligned region gets " ";
  for (my $i=1; $i < $hit_data_hr->{q_start}; $i++) {
    push @alignment, "-";
  }

  my $btop_align_r = decode_btop($hit_data_hr->{BTOP});

  my ($seq0, $seq1) = ("","");
  my ($qix, $six) = @{$hit_data_hr}{qw(q_start s_start)};

  for my $btop (@{$btop_align_r}) {
    if ($btop =~ m/^\d+$/) {  # matching query sequence, add it up
      for (my $i=0; $i < $btop; $i++) {
	if ($six >= $sb_start && $six <= $sb_end) {
	  push @alignment, $query_seq_r->[$qix];
	}
	else {
	  push @alignment, '-';
	}
	$qix++; $six++;
      }
    }
    else {  # could be: TS/-S/T-
      ($seq0, $seq1) = split(//,$btop);
      if ($seq1 eq '-') {  # gap in subject
	push @alignment, '-';
	$qix++;
      }
      elsif ($seq0 ne '-') {  # mismatch
	if ($six >= $sb_start && $six <= $sb_end) {
	  push @alignment, $seq1;
	}
	else {
	  push @alignment, '-';
	}
	$qix++;
	$six++;
      }
      else {  # gap in query, consume $six
	$six++;
      }
    }
  }
  # all done with alignment, double check that $qix = $hit_data_hr->{q_end}
  unless ($qix == $hit_data_hr->{q_end}+1) {
    warn $qix." != ".$hit_data_hr->{q_end}+1;
  }

  for (my $i = $hit_data_hr->{q_end}+1; $i <= $query_len; $i++) {
    push @alignment, "-";
  }

  return \@alignment;
}

sub parse_query_lib {
  my ($query_file) = @_;

  my %query_seqs = ();

  open(my $qfd, $query_file);


  my ($header, $sequence) = ("","");
  while (my $entry = <$qfd>) {  # returns an entire fasta entry
    chomp $entry;
    if ($entry =~ m/^>/) {
      $header = $entry;
    }
    else {
      $sequence .= $entry
    }
  }

  $sequence =~ s/[^A-Za-z\*]//g;    # remove everything but letters
  $sequence = uc($sequence);

  $header =~ s/^>//;
  $header =~ s/\s.*$//;
  my @seq = split(//,$sequence);
  unshift @seq,"";	# @seq is now 1-based

  return ($header, \@seq);
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
(-m 8).  The C<--ann_script/--q_ann_script> file is run to produce
domain boundary coordinates.  For searches against SwissProt
sequences, C<--ann_script ann_feats_up_www2.pl> will acquire features
and domains from Uniprot.  C<--ann_script ann_pfam_www.pl --neg> will
get domain information from Pfam, and score non-domain (NODOM)
regions.

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

=head1 AUTHOR

William R. Pearson, wrp@virginia.edu

=cut
