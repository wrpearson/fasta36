#!/usr/bin/env perl

################################################################
# copyright (c) 2014, 2015 by William R. Pearson and The Rector &
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

# ann_pdb_vast.pl gets an annotation file from fasta36 -V with a line of the form:

# gi|62822551|sp|P00502|GSTA1_RAT Glutathione S-transfer\n  (at least from pir1.lseg)
#
# it must:
# (1) read in the line
# (2) parse it to get the up_acc
# (3) return the tab delimited features
#

# this version is designed for various formats of the pdbaa/pdbaa_off NCBI files with the lines:
# >gi|4388890|pdb|1GTUA|sp|P09488  or 
# >gi|4388890|pdb|1GTU|A
# if I can find |sp|P09488, I will use that, otherwise I will use
# |pdb|1GTU|A (concatenated) and a different part of the cath_dom
# database
#

use warnings;
use strict;

use LWP::Simple;
use Getopt::Long;
use HTML::TableExtract;
use Pod::Usage;

use vars qw($host $db $port $user $pass);

my ($neg_doms, $lav, $shelp, $help, $class) = (0, 0, 0, 0, 0);
my ($min_nodom) = (10);

my $color_sep_str = " :";
$color_sep_str = '~';

GetOptions(
    "host=s" => \$host,
    "db=s" => \$db,
    "lav" => \$lav,
    "neg" => \$neg_doms,
    "neg_doms" => \$neg_doms,
    "neg-doms" => \$neg_doms,
    "min_nodom=i" => \$min_nodom,
    "h|?" => \$shelp,
    "help" => \$help,
    );

pod2usage(1) if $shelp;
pod2usage(exitstatus => 0, verbose => 2) if $help;
pod2usage(1) unless (@ARGV || -p STDIN || -f STDIN);

################################################################
# strategy for connecting to NCBI to get list of domains

my $db = "structure";
my $report = "FASTA";
my $pdb_acc_chain = "";

my $utils = "http://www.ncbi.nlm.nih.gov/entrez/eutils";
my $vast_url = "http://www.ncbi.nlm.nih.gov/Structure/vastplus/vastplus.cgi?cmd=v&uid=";

my $esearch = "$utils/esearch.fcgi?db=structure&retmax=1&term=";

################################################################

my %domains = (NODOM=>0);
my $domain_cnt = 0;

my ($tmp, $gi, $sdb, $acc, $id, $use_acc);

# get the query
my ($query, $seq_len) = @ARGV;
$seq_len = 1 unless $seq_len;

$query =~ s/^>// if ($query);

my @annots = ();

#if it's a file I can open, read and parse it
# unless ($query && $query =~ m/\|/) {

if (! $query || -r $query) {
  while (my $a_line = <>) {
    $a_line =~ s/^>//;
    chomp $a_line;
    push @annots, show_annots($lav,$a_line);
  }
}
else {
  push @annots, show_annots($lav, "$query $seq_len");
}

for my $seq_annot (@annots) {
  print ">",$seq_annot->{seq_info},"\n";
  for my $annot (@{$seq_annot->{list}}) {
    if (!$lav && defined($domains{$annot->[-1]})) {
      $annot->[-1] .= $color_sep_str.$domains{$annot->[-1]};
    }
    print join("\t",@$annot),"\n";
  }
}

exit(0);

sub show_annots {
  my ($lav, $annot_line) = @_;

  my ($query, $seq_len) = split(/\s+/,$annot_line);

  my %annot_data = (seq_info => $query);

  my ($tmp, $gi, $pdb, $pdb_acc, $pdb_id, $pdb_chain, $sdb, $up_acc, $off_flag);

  $off_flag = 0;
  if ($query =~ m/^gi\|/) {
    if ($query =~ m/\|sp\|/) {
      ($tmp, $gi, $pdb, $pdb_acc, $sdb, $up_acc) = split(/\|/,$query);
      $up_acc =~ s/\.\d+$//;
      $off_flag = 1;
    }
    elsif ($query =~ m/\|pdb\|/) {
      ($tmp, $gi, $pdb, $pdb_id, $pdb_chain) = split(/\|/,$query);
      $pdb_acc = $pdb_id . $pdb_chain;
    }
  }
  elsif ($query =~ m/^sp\|/) {
    ($pdb, $pdb_acc) = split(/\|/,$query);
  }
  elsif  ($query =~ m/^pdb\|(\w{4})\|(\w)/) {
    $pdb_acc = $1 . $2;
  }
  else {
    $pdb_acc = $query;
  }

# only get the first res_beg because it is used to calculate pdbaa_off @c:xxx

  $annot_data{list} = get_vast_annots($lav, $pdb_acc, $seq_len, $off_flag);

  return \%annot_data;
}

sub get_vast_annots {
  my ($lav, $pdb_acc, $seq_length) = @_;

  my $domain_href = get_vast_info($pdb_acc);

  return unless (scalar(@$domain_href));

  # do a consistency check
  for (my $i=1; $i < scalar(@$domain_href); $i++) {
    if ($domain_href->[$i]{seq_start} <= $domain_href->[$i-1]{seq_end}) {
      my $delta = $domain_href->[$i]{seq_start} - $domain_href->[$i-1]{seq_end};
      $domain_href->[$i-1]{seq_end} -= $delta/2;
      $domain_href->[$i]{seq_start} = $domain_href->[$i-1]{seq_end}+1;
    }
  }

  if ($neg_doms) {
    my @nvast_domains;
    my $prev_dom={seq_end=>0};
    for my $cur_dom ( @$domain_href) {
      if ($cur_dom->{seq_start} - $prev_dom->{seq_end} > $min_nodom) {
	my %new_dom = (seq_start=>$prev_dom->{seq_end}+1, seq_end => $cur_dom->{seq_start}-1, info=>'NODOM');
	push @nvast_domains, \%new_dom;
      }
      push @nvast_domains, $cur_dom;
      $prev_dom = $cur_dom;
    }
    my %new_dom = (seq_start=>$prev_dom->{seq_end}+1, seq_end=>$seq_length, info=>'NODOM');
    if ($new_dom{seq_end} > $new_dom{seq_start}) {push @nvast_domains, \%new_dom;}

    $domain_href = \@nvast_domains;
  }

  my @feats = ();

  if ($lav) {
    for my $d_ref (@$domain_href) {
      push @feats, [$d_ref->{seq_start}, $d_ref->{seq_end}, $d_ref->{info} ];
    }
  }
  else {
    for my $d_ref (@$domain_href) {
      push @feats, [$d_ref->{seq_start}, '-', $d_ref->{seq_end}, $d_ref->{info} ];
    }
  }

  return \@feats;

}

################################################################
# get_vast_info ( pdb_acc_chain pdb|1ABC|D
#

sub get_vast_info {
  my $pdb_acc_chain = shift @_;

  $pdb_acc_chain =~ s/^pdb\|//;

  my ($pdb_acc, $pdb_chain) = ($pdb_acc_chain =~ m/(\w{4})(\w)/);

  my $esearch_result = get($esearch . $pdb_acc . "[pdb+accession]");

  #  print "\nESEARCH RESULT: $esearch_result\n";

  my ($Count) =  ($esearch_result =~  m|<Count>(\d+)</Count>|s);

  my $mmdb_id = 0;

  ($mmdb_id) = ($esearch_result =~  m|<Id>(\d+)</Id>|s);

  my $vast_dom_html = get($vast_url.$mmdb_id);

  my $te = HTML::TableExtract->new(depth=>0, count=>0);

  $te->parse($vast_dom_html);

  my ($ts) = $te->tables();

  my @ts_rows = $ts->rows();

  my $row_header = shift(@ts_rows);

  # print "header:\t\t",join("\t",@$row_header),"\n";

  my $current_chain = "";

  my @domain_list = ();

  for my $row ( @ts_rows ) {
    if ($row->[0]) {
      ($current_chain) = ($row->[0] =~ m/^\[(\w+)\]/);
    }
    next unless ($current_chain eq $pdb_chain);
    next if ($row->[1] =~ m/^Entire/);
    my $range = $row->[2];

    # $range can actually be a list of ranges.  Need to remove the spaces
    $range =~ s/\s+//g;

    push @domain_list, $range;
  }

  # here we have the list of domain, indexed, but some domains will have multiple parts

  my @part_list = ();

  my $dom_cnt=1;
  for my $range (@domain_list) {
    my @parts = split(/,/,$range);
    for my $part ( @parts ) {
      my ($start, $end) = ($part =~ m/(\d+)\s*\-\s*(\d+)/);
      push @part_list, {info=>"VASTdom".$dom_cnt, seq_start=>$start, seq_end=>$end};
    }
    $dom_cnt++;
  }

  @part_list = sort { $a->{seq_start} <=> $b->{seq_start} } @part_list;

  return \@part_list;
}

__END__

=pod

=head1 NAME

ann_feats.pl

=head1 SYNOPSIS

 ann_pdb_vast.pl --neg 'sp|P09488|GSTM1_NUMAN' | accession.file

=head1 OPTIONS

 -h	short help
 --help include description

 --lav  produce lav2plt.pl annotation format, only show domains/repeats
 --neg-doms,  -- report domains between annotated domains as NODOM
                 (also --neg, --neg_doms)
 --min_nodom=10  -- minimum length between domains for NODOM

=head1 DESCRIPTION

C<ann_pdb_cath.pl> extracts domain information from the pfam26 msyql
database.  Currently, the program works with database sequence
descriptions in several formats:

 >pdb|3F6F|A	-- database|accession

>gi|262118558|pdb|3F6F|A  -- with GI number

C<ann_pdb_vast.pl> is designed to be used by the B<FASTA> programs with
the C<-V \!ann_pfam26.pl> or C<-V "\!ann_pfam26.pl --neg"> option.

=head1 AUTHOR

William R. Pearson, wrp@virginia.edu

=cut
