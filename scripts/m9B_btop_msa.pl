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
# m9B_btop_msa.pl --query query.file blast_tab_btop_file
################################################################
# m9B_btop_msa.pl takes a query sequence and a fasta -m 9b
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
################################################################

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

my ($shelp, $help, $evalue, $qvalue, $domain_bound) = (0, 0, 0.001, 30.0,0);
my ($query_file, $bound_file_in, $bound_file_only, $bound_file_out, $masked_lib_out) = ("","","","","");
my $query_lib_r = 0;

GetOptions(
    "query=s" => \$query_file,
    "query_file=s" => \$query_file,
    "evalue=f" => \$evalue,
    "expect=f" => \$evalue,
    "qvalue=f" => \$qvalue,
    "bound_file_in=s" => \$bound_file_in,
    "bound_file_only=s" => \$bound_file_only,
    "bound_file_out=s" => \$bound_file_out,
    "masked_library_out=s" => \$masked_lib_out,
    "masked_lib_out=s" => \$masked_lib_out,
    "domain_bound" => \$domain_bound,
    "domain" => \$domain_bound,
    "bound_in=s" => \$bound_file_in,
    "bound_only=s" => \$bound_file_only,
    "bound_out=s" => \$bound_file_out,
    "h|?" => \$shelp,
    "help" => \$help,
    );

pod2usage(1) if $shelp;
pod2usage(exitstatus => 0, verbose => 2) if $help;
unless (-f STDIN || -p STDIN || @ARGV) {
 pod2usage(1);
}

my @m9_field_names = qw(percid perc_sim raw_score a_len q_start q_end qc_start qc_end s_start s_end sc_start sc_end gap_q gap_l fs);

my @hit_list = ();
my %multi_align = ();
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
$multi_align{$query_acc} = btop2alignment($query_seq_r, $query_len, {BTOP=>$query_len, q_start=>1, q_end=>$query_len}, 0);
my $max_sseqid_len = length($query_acc);

################
# get sequence boundaries if available
#
my $seq_bound_hr = 0;
my @seq_bound_accs = ();

if ($bound_file_in) {
  $seq_bound_hr = parse_bound_file($bound_file_in);
}
elsif ($bound_file_only) {
  $seq_bound_hr = parse_bound_file($bound_file_only);
}
elsif ($domain_bound) {
  my %seq_bound = ();
  $seq_bound_hr = \%seq_bound;
}
elsif ($bound_file_out) {
  my %seq_bound = ();
  $seq_bound_hr = \%seq_bound;
}

################
# skip down to "The best scores are:"
#
my ($q_num, $query_descr, $q_len, $lib_cnt, $lib_len, $best_yes) = skip_to_results();
warn "Cannot find the best scores are:"  unless $query_descr;

my ($tmp, $gi, $q_db, $q_acc, $q_id);

if ($query_descr =~ /^gi\|\d+\|/) {
    ($tmp, $gi, $q_db,$q_acc, $q_id) = split(/\|/,$query_descr);
}
else {
    ($q_db,$q_acc, $q_id) = split(/\|/,$query_descr);
}

$q_acc =~ s/\.\d+$//;

while (my $line = <>) {
  chomp $line;
  next unless ($line);

  last if $line =~ m/>>>/;
  next if $line =~ m/^\+\-/; # skip over HSPs
  chomp ($line);

  my %hit_data =();

  my ($left, $right, $align_f, $annot_f) = split(/\t/,$line);

  $align_f= 'NULL' unless $align_f;
  $annot_f= 'NULL' unless $annot_f;

  my @fields = split(/\s+/,$left);
  my ($ldb, $l_id, $l_acc) = ("","","");
  if ($fields[0] =~ m/:/) {
      ($ldb, $l_id) = split(/:/,$fields[0]);
      ($l_acc) = $fields[1];
  }
  else {
      ($ldb, $l_acc,$l_id) = split(/\|/,$fields[0]);
  }

  @hit_data{@m9_field_names} = split(/\s+/,$right);
  @hit_data{qw(bits evalue)} = @fields[-2,-1];

  #
  # currently preselbdr files have $ldb|$l_acc, not full s_seqid, so construct it
  #
  my ($s_seqid, $subj_acc) = (join('|',($ldb, $l_acc, $l_id)), "$ldb|$l_acc");
  @hit_data{qw(s_seqid subj_acc)} = ($s_seqid, $subj_acc);
  @hit_data{qw(query_id query_acc)} = ($query_descr, $q_acc);
  $hit_data{BTOP} = $align_f;

  next if ($hit_data{evalue} > $evalue);

  if (length($s_seqid) > $max_sseqid_len) {
    $max_sseqid_len = length($s_seqid);
  }

  my $have_dom = 0;
  if ($domain_bound) {
    my $hit_doms_ar = parse_hit_domains($annot_f);
    # scan from left to right to make domain boundaries based on $qvalue
    my ($left_bound, $right_bound) = @hit_data{qw(s_end s_start)};
    foreach my $dom_r ( @$hit_doms_ar ) {
      next unless $dom_r->{target} eq 'subj';
      next if $dom_r->{virtual};
      next unless $dom_r->{qval} > $qvalue;

      if ($dom_r->{s_start} < $left_bound) {
	$left_bound = $dom_r->{s_start};
	$have_dom = 1;
      }

      if ($dom_r->{s_end} > $right_bound) {
	$right_bound = $dom_r->{s_end};
	$have_dom = 1;
      }
    }

    if ($have_dom) {
      if (exists($seq_bound_hr->{$subj_acc})) {
	@{$seq_bound_hr->{$subj_acc}}{qw(start end)} = ($left_bound, $right_bound);
      }
      else {
	$seq_bound_hr->{$subj_acc} = {start=>$left_bound, end=>$right_bound};
	push @seq_bound_accs, $subj_acc;
      }
    }
  }

  # must have separate @hit_list that can be sorted, for searches with multiple alignment results

  if ($bound_file_only || $have_dom) {
    if (exists($seq_bound_hr->{$subj_acc})) {
      my ($status, $alignment) = bound_btop2alignment($query_seq_r, $query_len, \%hit_data, @{$seq_bound_hr->{$subj_acc}}{qw(start end)});
      if ($status) {	# aligment is within boundary
	push @multi_names, $s_seqid;
	$multi_align{$s_seqid} = $alignment;
      }
      # do not delete entry, because it needs to be preserved 
    }
  }
  elsif ($bound_file_in) {
    if (exists($seq_bound_hr->{$subj_acc})) {
      my ($status, $alignment) = bound_btop2alignment($query_seq_r, $query_len, \%hit_data, @{$seq_bound_hr->{$subj_acc}}{qw(start end)});
      if ($status) {
	push @multi_names, $s_seqid;
	$multi_align{$s_seqid} = $alignment;
#	push @multi_align, $alignment;
      }
    }
    else {
      push @multi_names, $s_seqid;
#      push @multi_align, btop2alignment($query_seq_r, $query_len, \%hit_data, );
      $multi_align{$s_seqid} = btop2alignment($query_seq_r, $query_len, \%hit_data);
      @{$seq_bound_hr->{$subj_acc}}{qw(start end)} = @hit_data{qw(s_start s_end)};
      push @seq_bound_accs, $subj_acc;
    }
  }
  else {  # no sequence boundaries
    push @multi_names, $s_seqid;
    $multi_align{$s_seqid} = btop2alignment($query_seq_r, $query_len, \%hit_data);
#    push @multi_align, btop2alignment($query_seq_r, $query_len, \%hit_data);
    if (!$have_dom && ($bound_file_out)) {
      @{$seq_bound_hr->{$subj_acc}}{qw(start end)} = @hit_data{qw(s_start s_end)};
      push @seq_bound_accs, $subj_acc;
    }
  }
}

# final MSA output
$max_sseqid_len += 2;

print "SSEARCHm9B multiple sequence alignment\n\n\n";

my $i_pos = 0;
for (my $j = 0; $j < $query_len/60; $j++) {
  my $i_end = $i_pos + 59;
  if ($i_end > $query_len) {$i_end = $query_len-1;}
  for my $acc (@multi_names) {
    next unless $acc;
    printf("%-".$max_sseqid_len."s %s\n",$acc,join("",@{$multi_align{$acc}}[$i_pos .. $i_end]));
  }
  $i_pos += 60;
  print "\n\n";
}

################
# if bound_file_out provide it
if ($bound_file_out) {
  open(my $bound_fd, ">", $bound_file_out) || die "cannot open $bound_file_out";
  for my $s_acc ( @seq_bound_accs ) {
    print $bound_fd join("\t", ($s_acc, @{$seq_bound_hr->{$s_acc}}{qw(start end)})),"\n";
  }
  close($bound_fd);
}

if ($masked_lib_out) {
  open(my $masked_fd, ">", $masked_lib_out) || die "cannot open $masked_lib_out";

  for my $s_acc ( @multi_names ) {
    print $masked_fd ">$s_acc\n";
    my $seq_lines = join('',@{$multi_align{$s_acc}});

    # currently, simply remove all '-' insertions --
    # other options would be to make external '-'s 'X's
    $seq_lines =~ s/\-//g;

    $seq_lines =~ s/(.{60})/$1\n/g; 
    print $masked_fd "$seq_lines\n";
  }
  close($masked_fd);
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

sub parse_hit_domains {
  my ($annot_str) = @_;

## annot_str looks like: "|RX:6-65:6-65:s=311;b=125.4;I=1.000;Q=339.6;C=C.HTH~1
#                         |XR:6-65:6-65:s=311;b=125.4;I=1.000;Q=339.6;C=C.HTH~1
#                         |RX:66-297:66-297:s=1200;b=483.7;I=1.000;Q=1409.6;C=NODOM~0
#                         |XR:66-297:66-297:s=1200;b=483.7;I=1.000;Q=1409.6;C=NODOM~0

  return 0 unless ($annot_str);

  my @hit_annots = ();

  my @annots = split(/\|/,$annot_str);
  shift @annots;	# remove first blank

  for my $annot ( @annots ) {
    my %dom_info = ();

    # parse an entry:
    # |RX:6-65:6-65:s=311;b=125.4;I=1.000;Q=339.6;C=C.HTH~1
    my @d_fields = split(";",$annot);

    ($dom_info{dom}) = ($d_fields[4] =~ m/C=(.+?)~?\d*v?$/);   # also remove virtual domain symbols
    next if ($dom_info{dom} =~ m/NODOM/);

    ################
    # parse @d_fields
    if ($d_fields[4] =~ m/v$/) {
      $dom_info{virtual} = 1;
    }
    else {
      $dom_info{virtual} = 0;
    }

    ($dom_info{bits}) = ($d_fields[1] =~ m/b=(\-?\d+\.?\d*)/);
    unless (defined($dom_info{bits})) {
	warn "missing score info - annot: $annot\n  annot_str: $annot_str";
	$dom_info{bits} = '\N';
    }
    ($dom_info{percid}) = ($d_fields[2] =~ m/I=(\-?[\d\.]+)/);
    unless (defined($dom_info{percid})) {
	warn "missing percid info - annot: $annot\n  annot_str: $annot_str";
	$dom_info{percid} = '\N';
    }

    ($dom_info{qval}) = ($d_fields[3] =~ m/Q=([\d\.]+)/);

    ################
    # parse @c_fields
    my @c_fields = split(":",$d_fields[0]);

    if ($c_fields[0] =~ m/RX/) {$dom_info{target} = 'query';}
    else {$dom_info{target} = 'subj';}

    @dom_info{qw(q_start q_end)} = ($c_fields[1] =~ m/(\d+)\-(\d+)/);
    @dom_info{qw(s_start s_end)} = ($c_fields[2] =~ m/(\d+)\-(\d+)/);
    ($dom_info{score}) = ($c_fields[3] =~ m/s=(\-?\d+)/);
    unless (defined($dom_info{score})) {
	warn "missing score info - annot: $annot\n  annot_str: $annot_str";
	$dom_info{score} = '\N';
    }

    push @hit_annots, \%dom_info;
  }

  return \@hit_annots;
}


sub btop2alignment {
  my ($query_seq_r, $query_len, $hit_data_hr, $seq_bound_hr) = @_;

  # $query_seq_r is 1: based
  my @alignment = ();

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

################
# generates MSA alignment entry between $sb_start and $sb_end
# if there are no aligned residues between these locations, return $status=0

sub bound_btop2alignment {
  my ($query_seq_r, $query_len, $hit_data_hr, $sb_start, $sb_end) = @_;

  # $query_seq_r is 1: based
  my @alignment = ();

  my $have_aligned_res = 0;

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
	  $have_aligned_res=1;
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
	  $have_aligned_res=1;
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

  return ($have_aligned_res, \@alignment);
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

sub parse_bound_file {
  my ($bound_file) = @_;

  my %seq_bound = ();

  open(my $qfd, $bound_file) || return 0;

  while (my $line = <$qfd>) {
    next if ($line =~ m/^#/);
    chomp $line;
    my @data = split(/\t/,$line);
    if (!defined($seq_bound{$data[0]})) {
      $seq_bound{$data[0]} = {start=>$data[1], end=>$data[2]};
      push @seq_bound_accs, $data[0];
    }
    else {
      warn "multiple boundaries for $data[0]";
    }
  }

  return \%seq_bound;
}

sub skip_to_results {

  my ($q_num, $query_desc, $q_start, $q_stop, $q_len, $l_num, $l_len, $best_yes);

  while (my $line = <>) {
    if ($line =~ m/^\s*(\d+)>>>(\S+)\s.+ \- (\d+) aa$/) {
      ($q_num,$query_desc, $q_len) = ($1,$2,$3);
#      ($q_len) = ($line =~ m/(\d+) aa$/);
      $line = <>;	# skip Library:
      $line = <>;	#   153571012 residues in 291716 sequences
      ($l_len, $l_num) = ($line =~ m/^\s+(\d+)\s+residues in\s+(\d+)/);
      goto have_query;
    }
    elsif ($line =~ m/>>>\/\/\//) {goto done;}
  }
 done:
  return (0,"");

 have_query:
  while (my $line = <>) {
    $best_yes = 0;

    if ($line =~ m/^The best scores are:/) {
      $best_yes = 1;
      last;
    }
    last if ($line =~ m/^!! No sequences/);
  }
  return ($q_num, $query_desc,$q_start, $q_stop, $q_len, $l_num, $l_len, $best_yes);
}

__END__

=pod

=head1 NAME

 m9B_btop_msa.pl

=head1 SYNOPSIS

 m9B_btop_msa.pl --query_file query.fasta [--bound_file seqbdr.tab] fasta_m9_output.file

=head1 OPTIONS

 -h	short help
 --help include description

 --query_file -- query sequence file
              -- same as --query
                 (only one sequence per file)

 --bound_file_in -- tab delimited accession<tab>start<tab>end that
                    specifies MSA boundaries WITHIN alignment.
                    Additional hits use alignment (or domain)
                    boundaries.

 --bound_file_only -- tab delimited accession<tab>start<tab>end that
                      specifies MSA boundaries WITHIN alignment.
                      Only sequences in --bound_file_only will be in the MSA.

 --bound_file_out -- "--bound_file" for next iteration of psisearch2

 --domain_bound  parse domain annotations (-V) from m9B file
 --domain

 --masked_lib_out -- FASTA format library of MSA sequences

=head1 DESCRIPTION

C<m9B_btop_msa.pl> takes a fasta36/ssearch36 -m 9B ouput file, which
includes a BTOP encoded alignment string, and produces the multiple
sequence alignment (MSA) implied by the query sequence, alignment
boundaries, and pairwise alignments.  The alignment does not allow
gaps in the query sequence, only in the subject sequences.

The C<--query_file> must be specified, and the query sequence is
provided as the first sequence in the MSA.

If a C<--bound_file> is provided, then the ends of the alignments are
reduced to the coordinates specified by the C<bound_file>.  In
addition, only sequences included in the C<bound_file> are included in
the MSA.

Output: A clustal-like interleaved multiple sequence alignment that
can be used as input (using the C<-in_msa> option) to C<psiblast>.

=head1 AUTHOR

William R. Pearson, wrp@virginia.edu

=cut
