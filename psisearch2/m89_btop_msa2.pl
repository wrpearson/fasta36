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
# m89_btop_msa2.pl --query query.file blast_tab_btop_file
################################################################
# m89_btop_msa2.pl takes a query sequence and either a BLAST -outfmt
# 7/fasta -m 8CB output file (default) or fasta -m 8B (--m_format m9)
# output file with a BTOP field, and constructs a query-driven
# multiple sequence alignment of the subject sequences that can be
# used as input to psiblast with the "--in_msa msa.file" option.
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

use warnings;
use strict;
use Pod::Usage;
use Getopt::Long;

# read lines of the form:
# gi|121694|sp|P20432.1|GSTT1_DROME	gi|121694|sp|P20432|GSTT1_DROME	100.00	209	0	0	1	209	1	209	6e-156	433	1113	209
# gi|121694|sp|P20432.1|GSTT1_DROME	gi|1170090|sp|P04907|GSTF3_MAIZE	26.77	198	123	7	4	185	6	197	2e-08	51.2	121	FL1YG ... 1NKRA1YW1
# gi|121694|sp|P20432.1|GSTT1_DROME	gi|81174731|sp|P0ACA5|SSPA_ECO57	39.66	58	32	2	43	100	49	103	8e-06	43.9	102	EDFLLI ... V-I-NEQS3FM
# gi|121694|sp|P20432.1|GSTT1_DROME	gi|121695|sp|P12653|GSTF1_MAIZE	27.62	181	107	7	32	203	34	199	9e-05	40.8	94	LI1LF ... N-1AS1CLLM1


my ($shelp, $help, $m_format, $evalue, $qvalue, $domain_bound) = (0, 0, "m8CB", 0.001, 30.0,0);
my ($query_file, $sel_file, $bound_file_in, $bound_file_only, $bound_file_out, $masked_lib_out,$mask_type_end, $mask_type_int) = ("","","","","","","","");
my ($clustal_id,$trunc_acc,$min_align) = (0,0,0.0);
my $query_lib_r = 0;
my ($eval2_fmt, $eval2) = (0,"");

GetOptions(
    "query=s" => \$query_file,
    "query_file=s" => \$query_file,
    "eval2=s" => \$eval2,		# change the evalue used for inclusion
    "evalue|expect=f" => \$evalue,
    "qvalue=f" => \$qvalue,
    "format=s" => \$m_format,
    "clustal!" => \$clustal_id,
    "trunc_acc!" => \$trunc_acc,
    "selected_file_in|sel_file_in|sel_accs=s" => \$sel_file,
    "m_format|mformat=s" => \$m_format,
    "bound_file_in=s" => \$bound_file_in,
    "bound_file_only=s" => \$bound_file_only,
    "bound_file_out=s" => \$bound_file_out,
    "bound_in=s" => \$bound_file_in,
    "bound_only=s" => \$bound_file_only,
    "bound_out=s" => \$bound_file_out,
    "masked_library_out=s" => \$masked_lib_out,
    "masked_lib_out=s" => \$masked_lib_out,
    "mask_lib_out=s" => \$masked_lib_out,
    "mask_out=s" => \$masked_lib_out,
    "end_mask_type=s" => \$mask_type_end,
    "end_mask=s" => \$mask_type_end,
    "domain_bound" => \$domain_bound,
    "domain" => \$domain_bound,
    "int_mask_type=s" => \$mask_type_int,
    "int_mask=s" => \$mask_type_int,
    "min_align=f" => \$min_align,
    "h|?" => \$shelp,
    "help" => \$help,
    );

pod2usage(1) if $shelp;
pod2usage(exitstatus => 0, verbose => 2) if $help;
unless (-f STDIN || -p STDIN || @ARGV) {
 pod2usage(1);
}

################
# initialization
my @random_res = ();
if ($mask_type_end =~ m/^rand/i) {
  @random_res = init_random_res();
}

my @m9_field_names = qw(percid perc_sim raw_score a_len q_start q_end qc_start qc_end s_start s_end sc_start sc_end gap_q gap_l fs);
my @m8_field_names = qw(q_seqid s_seqid percid a_len mismatch gopen q_start q_end s_start s_end evalue bits);

my @hit_list = ();

my %seq_bound = ();  # boundaries for each accession
my %acc_names = ();  # generate uniq s_seq_id names
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

my %sel_accs = ();
my $have_sel_accs = 0;
if ($sel_file) {
  if (open (my $sfd, $sel_file)) {
    while (my $line = <$sfd>) {
      chomp $line;
      next if $line =~ m/^#/;
      if ($line =~ m/\t/) {
	my @data = split(/\t/,$line);
	$sel_accs{$data[0]} = $data[1];
      }
      else {
	$sel_accs{$line} = 1;
      }
    }
    close($sfd);
    $evalue = 1000000.0;
    $have_sel_accs = 1;
  }
  else {
    warn "Cannot open selected sequence file: $sel_file";
  }
}

push @multi_names, $query_acc;
$acc_names{$query_acc} = 1;
$multi_align{$query_acc} = btop2alignment($query_seq_r, $query_len, {BTOP=>$query_len, q_start=>1, q_end=>$query_len}, 0);
my $max_sseqid_len = length($query_acc);

################
# get sequence boundaries if available
#
my $seq_bound_hr = 0;
my @seq_bound_accs = ();

$seq_bound_hr = \%seq_bound;

if ($bound_file_in) {
  $seq_bound_hr = parse_bound_file($bound_file_in);
}
elsif ($bound_file_only) {
  $seq_bound_hr = parse_bound_file($bound_file_only);
}

################
# skip down to "The best scores are:"
#
my ($q_num, $query_descr, $q_len, $lib_cnt, $lib_len, $best_yes, $last_fields_r);

if ($m_format =~ m/^m9/i) {
  ($q_num, $query_descr, $q_len, $lib_cnt, $lib_len, $best_yes) = skip_to_m9results();
  warn "Cannot find the best scores are:"  unless $query_descr;
}
elsif ($m_format =~ m/^m8/) {
  ($query_descr, $best_yes, $last_fields_r) = skip_to_m8results();
  warn "Cannot find the best scores are: in $ARGV" unless $query_descr;

  if (scalar(@{$last_fields_r})) {
    push @m8_field_names, @{$last_fields_r};
  }
  push @m8_field_names, "annot";
}
else {
  die "cannot recognize format: $m_format";
}

my $eval_fptr = \&eval_func;
if ($eval2_fmt && $eval2) {
  if($eval2 eq 'eval2') {
    $eval_fptr = \&eval2_func;
  }
  elsif ($eval2 eq 'ave') {
    $eval_fptr = \&eval_ave;
  }
}

my ($tmp, $gi, $q_db, $q_acc, $q_id);

if ($query_descr =~ /^gi\|\d+\|/) {
  ($tmp, $gi, $q_db, $q_acc, $q_id) = split(/\|/,$query_descr);
}
elsif ($query_descr eq 'unnamed') {
  $q_acc = 'unnamed'
}
else {
  ($q_db,$q_acc, $q_id) = split(/\|/,$query_descr);
}

unless ($q_acc) {
  $q_acc = $query_descr;
}

$acc_names{$q_acc} = $q_acc;	# this is necessary for the new acc-only NCBI SwissProt libraries

$q_acc =~ s/\.\d+$//;

while (my $line = <>) {

  chomp $line;
  next unless ($line);
  my %hit_data =();
  my ($s_seqid, $subj_acc, $s_seqid_u);
  my $annot_f='NULL';

  if ($m_format =~ m/^m9/i) {
    last if $line =~ m/>>>/ || $line =~ m/^<\/pre>/;
    next if $line =~ m/^\+\-/; # skip over HSPs
    my ($left, $right, $align_f) = ("","",'NULL');
    ($left, $right, $align_f, $annot_f) = split(/\t/,$line);

    $align_f= 'NULL' unless $align_f;
    $annot_f= 'NULL' unless $annot_f;

    if ($left =~ m/<font/) {
      $left =~ s/<font color="darkred">//;
      $left =~ s/<\/font>//;
    }

    my @fields = split(/\s+/,$left);
    $subj_acc = $s_seqid = $fields[0];

    # my ($ldb, $l_id, $l_acc) = ("","","");
    # if ($fields[0] =~ m/:/) {
    #   ($ldb, $l_id) = split(/:/,$fields[0]);
    #   ($l_acc) = $fields[1];
    # } else {
    #   ($ldb, $l_acc,$l_id) = split(/\|/,$fields[0]);
    # }

    @hit_data{@m9_field_names} = split(/\s+/,$right);

    if ($eval2_fmt) {
      @hit_data{qw(bits evalue eval2)} = @fields[-3, -2,-1];
    }
    else {
      @hit_data{qw(bits evalue)} = @fields[-2,-1];
    }

    #
    # currently preselbdr files have $ldb|$l_acc, not full s_seqid, so construct it
    #
    # ($s_seqid, $subj_acc) = (join('|',($ldb, $l_acc, $l_id)), "$ldb|$l_acc");
    @hit_data{qw(s_seqid subj_acc)} = ($s_seqid, $subj_acc);
    @hit_data{qw(query_id query_acc)} = ($query_descr, $q_acc);
    $hit_data{BTOP} = $align_f;
    $hit_data{annot} = $annot_f;
  }
  else {
    last if $line =~ m/^#/;
    @hit_data{@m8_field_names} = split(/\t/,$line);
    $subj_acc = $hit_data{'s_seqid'};
    # remove gi number
    if ($subj_acc =~ m/^gi|\d+\|/) {
      $subj_acc =~ s/^gi\|\d+\|//;
    }
  }

  if ($have_sel_accs) {
    next unless ($sel_accs{$hit_data{'s_seqid'}});
  }

  # a better solution would be to rename the q_seqid, or at least to
  # check for identity
  # next if ($hit_data{q_seqid} eq $hit_data{s_seqid});
  next if ($hit_data{a_len} == $query_len && $hit_data{BTOP} =~ m/^$query_len$/);

  $s_seqid_u = $hit_data{'s_seqid'};
  if ($acc_names{$s_seqid_u}) {
    next;  # skip additional HSPs
#    $acc_names{$s_seqid_u}++;
#    $s_seqid_u .= "_". $acc_names{$subj_acc};
  }
  else {
    my $tr_acc = $hit_data{'s_seqid'};
    $acc_names{$hit_data{'s_seqid'}} = 1;
  }

  # must be after duplicate seqid check because blast HSP's have bad E-values after good.
  next if ($eval_fptr->(\%hit_data) > $evalue);

  next if (($hit_data{q_end}-$hit_data{q_start}+1)/$query_len < $min_align);

  $hit_data{s_seqid_u} = $s_seqid_u;

  my $have_dom = 0;
  if ($domain_bound && $hit_data{annot}) {
    my $hit_doms_ar = parse_hit_domains($hit_data{annot});
    # scan from left to right to make domain boundaries based on $qvalue
    # the following seems reversed, but it is putting upper (and lower) limits on boundaries
    my ($left_bound, $right_bound) = @hit_data{qw(s_end s_start)};
    foreach my $dom_r ( @$hit_doms_ar ) {
      next unless $dom_r->{target} eq 'subj';
      # next if $dom_r->{virtual};	# should be controlled by annotation process
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
	push @multi_names, $hit_data{s_seqid_u};
	$multi_align{$hit_data{s_seqid_u}} = $alignment;
      }
      # do not delete entry, because it needs to be preserved
    }
  }
  elsif ($bound_file_in) {
    if (exists($seq_bound_hr->{$subj_acc})) {
      my ($status, $alignment) = bound_btop2alignment($query_seq_r, $query_len, \%hit_data, @{$seq_bound_hr->{$subj_acc}}{qw(start end)});
      if ($status) {
	push @multi_names, $hit_data{s_seqid_u};
	$multi_align{$hit_data{s_seqid_u}} = $alignment;
#	push @multi_align, $alignment;
      }
    }
    else {
      push @multi_names, $hit_data{s_seqid_u};
#      push @multi_align, btop2alignment($query_seq_r, $query_len, \%hit_data, );
      $multi_align{$hit_data{s_seqid_u}} = btop2alignment($query_seq_r, $query_len, \%hit_data);
      @{$seq_bound_hr->{$subj_acc}}{qw(start end)} = @hit_data{qw(s_start s_end)};
      push @seq_bound_accs, $subj_acc;
    }
  }
  else {  # no sequence boundaries
    push @multi_names, $hit_data{s_seqid_u};
    $multi_align{$hit_data{s_seqid_u}} = btop2alignment($query_seq_r, $query_len, \%hit_data);
#    push @multi_align, btop2alignment($query_seq_r, $query_len, \%hit_data);
    if (!$have_dom && ($bound_file_out)) {
      @{$seq_bound_hr->{$subj_acc}}{qw(start end)} = @hit_data{qw(s_start s_end)};
      push @seq_bound_accs, $subj_acc;
    }
  }
}

$max_sseqid_len = 10;
for my $acc ( @multi_names) {
  my $this_len = length($acc);
  if ($trunc_acc && ($acc=~m/\|\w+\|(\w+)$/)) {
    $this_len = length($1);
  }
  if ($this_len > $max_sseqid_len) {
    $max_sseqid_len = $this_len;
  }
}

# final MSA output
$max_sseqid_len += 2;

if (! $clustal_id) {
  printf "BTOP%s multiple sequence alignment\n\n\n",$m_format;
}
else {
  print "CLUSTALW (1.8) multiple sequence alignment\n\n\n";
}

my $i_pos = 0;
for (my $j = 0; $j < $query_len/60; $j++) {
  my $i_end = $i_pos + 59;
  if ($i_end >= $query_len) {$i_end = $query_len-1;}
  for my $acc (@multi_names) {
    next unless $acc;

    my $this_acc = $acc;
    if ($trunc_acc && ($acc=~m/\|\w+\|(\w+)$/)) {
      $this_acc = $1;
    }
    printf("%-".$max_sseqid_len."s %s\n",$this_acc,join("",@{$multi_align{$acc}}[$i_pos .. $i_end]));
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

    # here we have four choices for masking:
    # (1) simply delete the '-'s
    # (2) delete the leading/trailing '-',s replace interal '-'s with 'X'
    # (3) replace leading/trailing '-' with 'X', remove internal
    # (4) replace leading/trailing '-' with query sequence, remove internal
    # (5) replace leading/trailing '-' with random, remove internal
    # (6) replace leading/trailing '-' with random, internal with 'X'

#    my @masked_seq = @{$multi_align{$s_acc}};
    my $seq = join('',@{$multi_align{$s_acc}});

    my @masked_seq = @{$multi_align{$s_acc}};
    my $n_res = scalar(@masked_seq);
    my $n_rand_res = scalar(@random_res);
    if ($mask_type_end =~ m/x/i) {
      for (my $i=0; $i < $n_res; $i++) {
	last if ($masked_seq[$i] ne '-') ;
	$masked_seq[$i] = 'X';
      }
      for (my $i=$n_res-1; $i >= 0; $i--) {
	last if ($masked_seq[$i] ne '-') ;
	$masked_seq[$i] = 'X';
      }
    }
    elsif ($mask_type_end =~ m/^q/i) {
      if ($mask_type_int =~ m/^q/i) {
	for (my $i=0; $i < $n_res; $i++) {
	  if ($masked_seq[$i] eq '-') {
	    $masked_seq[$i] = $multi_align{$query_acc}[$i];
	  }
	}
      }
      else {
	my $li = 0;
	for ( ; $li < $n_res; $li++) {
	  last if ($masked_seq[$li] ne '-') ;
	  $masked_seq[$li] = $multi_align{$query_acc}[$li];
	}
	my $ri = $n_res-1;
	for ( ; $ri >= 0; $ri--) {
	  last if ($masked_seq[$ri] ne '-') ;
	  $masked_seq[$ri] = $multi_align{$query_acc}[$ri];
	}
	if ($mask_type_int =~ m/^rand/i) {
	  for (my $i=$li; $i <= $ri; $i++) {
	    if ($masked_seq[$i] eq '-') {
	      $masked_seq[$i] = $random_res[int(rand($n_rand_res))];
	    }
	  }
	}
      }
    } elsif ($mask_type_end =~ m/^rand/i) {
      if ($mask_type_int =~ m/^rand/i) {
	for (my $i=0; $i < $n_res; $i++) {
	  if ($masked_seq[$i] eq '-') {
	    $masked_seq[$i] = $random_res[int(rand($n_rand_res))];
	  }
	}
      }
      else {
	for (my $i=0; $i < $n_res; $i++) {
	  last if ($masked_seq[$i] ne '-') ;
	  $masked_seq[$i] = $random_res[int(rand($n_rand_res))];
	}
	for (my $i=$n_res-1; $i >= 0; $i--) {
	  last if ($masked_seq[$i] ne '-') ;
	  $masked_seq[$i] = $random_res[int(rand($n_rand_res))];
	}
      }
    }

    my $masked_seq = join("",@masked_seq);
    if ($mask_type_int =~ m/X/) {
      $masked_seq =~ s/\-/X/g;
    }
    else {
      $masked_seq =~ s/\-//g;
    }

    $masked_seq =~ s/(.{60})/$1\n/g;
    print $masked_fd "$masked_seq\n";
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

  close($qfd);

  return ($header, \@seq);
}

sub parse_bound_file {
  my ($bound_file) = @_;

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

sub skip_to_m9results {

  my ($q_num, $query_desc, $q_start, $q_stop, $q_len, $l_num, $l_len, $best_yes);

  while (my $line = <>) {
    if ($line =~ m/^\s*(\d+)>>>(\S+)\s.*\- (\d+) aa$/) {
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
      $eval2_fmt = 1 if ($line =~ m/E2()/);
      last;
    }
    last if ($line =~ m/^!! No sequences/);
  }
  return ($q_num, $query_desc,$q_start, $q_stop, $q_len, $l_num, $l_len, $best_yes);
}

sub skip_to_m8results {

  my ($query_desc, $best_yes);

  $best_yes = 0;
  $eval2_fmt = 0;

  my @last_fields = ();

  while (my $line = <>) {
    if ($line =~ m/^# Query:/) {  # Query:
      ($query_desc) = ($line =~ m/^# Query:\s+(\S+)/);
      $query_desc = 'unnamed' unless ($query_desc);
#      ($q_len) = ($line =~ m/(\d+) aa$/);
      $line = <>;	# Database:
      $line = <>;	# Fields:

      unless ($line =~ m/# Fields:/) {
	warn "!!! warning !!!: # Fields not found: $line";
      }

      if ($line =~ m/,\seval2/) {		# only with FASTA
	push @last_fields, "eval2";
	$eval2_fmt = 1;
      }

      if ($line =~ m/,\sscore,\s+BTOP/i) {	# only with BLAST
	push @last_fields, qw(score BTOP);
      }
      elsif ($line =~ m/\BTOP,\s+score/i) {
	  push @last_fields, qw(BTOP score);
	}
      elsif ($line =~ m/,\s+BTOP/) {
	push @last_fields, qw(BTOP);
      }

      $line = <>; 	# NNN fits found or
      if ($line =~ m/^#\s+\d+\s+hits found/) {
	$best_yes = 1;
      }
      goto have_query;
    }
    elsif ($line =~ m/>>>\/\/\//) {goto done;}
  }
 done:
  return (0,"");

 have_query:
  return ($query_desc, $best_yes, \@last_fields);
}

################
# eval_func -- return evalue
sub eval_func {
    my ($hit) = @_;
    return $hit->{evalue};
}

# eval_func -- return evalue
sub eval2_func {
    my ($hit) = @_;
    return $hit->{eval2};
}
# eval_func -- return evalue
sub eval_ave {
    my ($hit) = @_;
    return sqrt($hit->{evalue}*$hit->{eval2});
}


################
# init_random_res initializes a 1000 element array of amino acid residues with Robinson/Robinson frequencies

sub init_random_res {
  my @rr_res = qw(A R N D C Q E G H I L K M F P S T W Y V);
  my @rr_counts = (35155, 23105, 20212, 24161,  8669,  19208,
		   28354, 33229,  9906, 23161, 40625,  25872,
		   10101, 17367, 23435, 32070, 26311,   5990,
		   14488, 29012);
  my $rr_total = 450431;

  my $rr_seq = "";
  for (my $i=0; $i < 20; $i++) {
    $rr_seq = $rr_seq . $rr_res[$i] x int(1000.0 *$rr_counts[$i]/$rr_total + 0.5);
  }
  return split(//,$rr_seq);
}

__END__

=pod

=head1 NAME

 m89_btop_msa2.pl

=head1 SYNOPSIS

 m89_btop_msa2.pl --query_file query.fasta fasta_m8CB_output.file
 m89_btop_msa2.pl --query_file query.fasta blast_outfmt7_BTOP.output.file
 m89_btop_msa2.pl --query_file query.fasta [--m_format m9] fasta_m9C_output.file

=head1 OPTIONS

 -h	short help
 --help include description

 --query_file -- query sequence file
 --query      -- same as --query_file
 (only one sequence per file)

 --expect|evalue: 0.001 -- maximum e-value to be include in output

 --eval2 : "": use E()-value, "eval2": use E2()/eval2, "ave": use geom. mean

 --qvalue: 30.0  -- minimum qvalue for domain to be considered

 --bound_file_in -- tab delimited accession<tab>start<tab>end that
                    specifies MSA boundaries WITHIN alignment.
                    Additional hits use alignment (or domain)
                    boundaries.

 --bound_file_only -- tab delimited accession<tab>start<tab>end that
                      specifies MSA boundaries WITHIN alignment.
                      Only sequences in --bound_file_only will be in the MSA.

 --bound_file_out -- "--bound_file" for next iteration of psisearch2

 --clustal -- use "CLUSTALW (1.8)" multiple alignment string

 --trunc_acc -- remove db, acc from db|acc|ident, e.g. sp|P0948|GSTM1_HUMAN becomes GSTM1_HUMAN

 --domain_bound  parse domain annotations (-V) from m9B file
 --domain

 --masked_lib_out -- FASTA format library of MSA sequences

 --min_align:0.0  -- minimum fractional alignment (q_end-q_start+1)/q_len

 --int_mask_type = "query", "rand", "X", "none"
 --end_mask_type = "query", "rand", "X", "none"
specify the residues to be inserted into output library

=head1 DESCRIPTION

C<m89_btop_msa2.pl> takes a fasta36/ssearch36 -m 9B ouput file, which
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

If an C<--masked_lib_out> filename.fasta is specified, a version of the MSA
in FASTA format is written to filename.fasta.  This file can be
converted to BLAST format (C<makeblastdb --in filename.fasta>) and
and the converted blast library can be used to rebuild a PSSM with
C<psiblast -num_iterations 2 -db filename.fasta -in_msa filename.msa
-out_pssm filename.asn_txt>.

The sequences in the C<--masked_lib_out> fasta file can be modifed
where gaps are present in the MSA as specified by the C<--int_mask_type> and C<--end_mask_type>
options. If no --int_mask_type/--end_mask_type is specified
("none"), then the subject sequence in the output library matches the
aligned part of the subject sequenc (gaps characters are deleted).  If
the --end/int_mask_type is "query", "rand", or "X", then either the
aligned query residue, a random residue, or an "X" substituted at
each gap, producing a library of subject sequences that may differ
from the original subject sequences.  These different options can be
used to force C<psiblast> to build a PSSM that more accurately
reflects the original C<ssearch36> alignment.

=head1 AUTHOR

William R. Pearson, wrp@virginia.edu

=cut
