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

# ann_exons_up_www.pl gets an annotation file from fasta36 -V with a line of the form:

# gi|23065544|ref|NP_000552.2| 
#
# and returns the exons present in the protein from NCBI gff3 tables (human and mouse only)
#
# it must:
# (1) read in the line
# (2) parse it to get the acc
# (3) get exon information from EBI/Uniprot
# (4) return the tab delimited exon boundaries

# 22-May-2017 -- use get("http://"), not get_https("https://"), because EBI does not have LWP::Protocol:https

use strict;

use Getopt::Long;
use LWP::Simple;
use LWP::UserAgent;
# use LWP::Protocol::https;
use Pod::Usage;
use JSON qw(decode_json);

use vars qw($host $db $port $user $pass);

my ($lav, $shelp, $help) = (0, 0,0);

my $color_sep_str = " :";
$color_sep_str = '~';

GetOptions(
    "lav" => \$lav,
    "h|?" => \$shelp,
    "help" => \$help,
    );

pod2usage(1) if $shelp;
pod2usage(exitstatus => 0, verbose => 2) if $help;
pod2usage(1) unless (@ARGV || -p STDIN || -f STDIN);

my %domains = (NODOM=>0);
my @domain_list = (0);
my $domain_cnt = 0;

my ($tmp, $gi, $sdb, $acc, $id, $use_acc);

my $get_annot_sub = \&get_up_www_exons;

my $ua = LWP::UserAgent->new(ssl_opts=>{verify_hostname => 0});
my $uniprot_url = 'http://www.ebi.ac.uk/proteins/api/coordinates/';
my $uniprot_suff = ".json";

# get the query
my ($query, $seq_len) = @ARGV;
$seq_len = 0 unless defined($seq_len);

$query =~ s/^>// if ($query);

my @annots = ();
my %annot_set = (); # re-use annotations if they are available (not yet implemented)

#if it's a file I can open, read and parse it
unless ($query && ($query =~ m/[\|:]/ ||
		   $query =~ m/^[OPQ][0-9][A-Z0-9]{3}[0-9]|[A-NR-Z][0-9]([A-Z][A-Z0-9]{2}[0-9]){1,2}\s/)) {

  while (my $a_line = <>) {
    $a_line =~ s/^>//;
    chomp $a_line;
    push @annots, show_annots($a_line, $get_annot_sub);
  }
}
else {
  push @annots, show_annots("$query\t$seq_len", $get_annot_sub);
}

for my $seq_annot (@annots) {
  print ">",$seq_annot->{seq_info},"\n";
  for my $annot (@{$seq_annot->{list}}) {
    print join("\t",@$annot),"\n";
  }
}

exit(0);

sub show_annots {
  my ($query_len, $get_annot_sub) = @_;

  my ($annot_line, $seq_len) = split(/\t/,$query_len);

  my %annot_data = (seq_info=>$annot_line);

  $use_acc = 1;

  if ($annot_line =~ m/^gi\|/) {
    ($tmp, $gi, $sdb, $acc, $id) = split(/\|/,$annot_line);
  }
  elsif ($annot_line =~ m/^(sp|tr|up)\|/) {
    ($sdb, $acc, $id) = split(/\|/,$annot_line);
  }
  elsif ($annot_line =~ m/^(SP|TR):(\w+) (\w+)/) {
      ($sdb, $id, $acc) = (lc($1), $2, $3);
  }
  elsif ($annot_line =~ m/^(SP|TR):(\w+)/) {
    ($sdb, $id, $acc) = (lc($1), $2, "");
    warn "*** $0 - accession required: $annot_line";
  }
  elsif ($annot_line =~ m/\|/) {
    ($sdb, $acc) = split(/\|/,$annot_line);
  }
  else {
    ($acc) = ($annot_line =~ m/^(\S+)/);
  }

  $acc =~ s/\.\d+$//;

  my  $exon_json = get($uniprot_url.$acc.$uniprot_suff);

  unless (!$exon_json || $exon_json =~ m/errorMessage/ || $exon_json =~ m/Can not find/) {
    $annot_data{list} = parse_json_up_exons($exon_json);
  }

  return \%annot_data;
}

sub parse_json_up_exons {
  my ($exon_json) = @_;

  my @exons = ();

  my $acc_exons = decode_json($exon_json);

  my $exon_num = 1;
  for my $exon ( @{$acc_exons->{'gnCoordinate'}[0]{'genomicLocation'}{'exon'}} ) {
    my ($p_begin, $p_end) = ($exon->{'proteinLocation'}{'begin'}{'position'},$exon->{'proteinLocation'}{'end'}{'position'});
    if ($p_end >= $p_begin) {
      push @exons, {
		    info=>"exon_".$exon_num.$color_sep_str.$exon_num,
		    seq_start=>$p_begin,
		    seq_end=>$p_end,
		   };
      $exon_num++;
    }
  }

  # check for domain overlap, and resolve check for domain overlap
  # (possibly more than 2 domains), choosing the domain with the best
  # evalue

  my @ex_feats = ();

  for my $d_ref (@exons) {
    if ($lav) {
      push @ex_feats, [$d_ref->{seq_start}, $d_ref->{seq_end}, $d_ref->{info}];
    }
    else {
      push @ex_feats, [$d_ref->{seq_start}, '-', $d_ref->{seq_end},  $d_ref->{info} ];
    }
  }
  return \@ex_feats;
}

sub get_https {
  my ($url) = @_;

  my $result = "";
  my $response = $ua->get($url);

  if ($response->is_success) {
    $result = $response->decoded_content;
  } else {
    $result = '';
  }
  return $result;
}

sub domain_name {

  my ($value) = @_;

  if (!defined($domains{$value})) {
    $domain_cnt++;
    $domains{$value} = $domain_cnt;
  }
  return $value;
}

__END__

=pod

=head1 NAME

ann_exons_up_www.pl

=head1 SYNOPSIS

 ann_exons_up_www.pl 'sp|P09488|GSTM1_HUMAN' | accession.file

=head1 OPTIONS

 -h	short help
 --help include description
 --lav  produce lav2plt.pl annotation format

=head1 DESCRIPTION

C<ann_exons_up_www.pl> extracts exon coordinates for proteins using
the EBI Proteins REST API described here:
C<https://www.ebi.ac.uk/proteins/api/doc/#coordinatesApi>.  Exon
intron boundaries, in protein coordinates, are available for Uniprot
proteins with Ensembl entries.

C<ann_pfam.pl> is designed to be used by the B<FASTA> programs with
the C<-V \!ann_exons_up_www.pl> or C<-V "q\!ann_exons_up_www.plg"> option.

=head1 AUTHOR

William R. Pearson, wrp@virginia.edu

=cut
