#!/usr/bin/perl -w

################################################################
# copyright (c) 2014 by William R. Pearson and The Rector &
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

# ann_pfam_www_e.pl gets an annotation file from fasta36 -V with a line of the form:

# gi|62822551|sp|P00502|GSTA1_RAT Glutathione S-transfer\n  (at least from pir1.lseg)
#
# it must:
# (1) read in the line
# (2) parse it to get the up_acc
# (3) return the tab delimited features
#

# This version uses the Pfam RESTful interface, rather than a local database
# >pf26|164|O57809|1A1D_PYRHO
# and only provides domain information

# use strict;

use Getopt::Long;
use Pod::Usage;
use LWP::Simple;
use XML::Twig;

my ($auto_reg,$rpd2_fams, $neg_doms, $lav, $no_doms, $pf_acc, $shelp, $help, $no_over, $pfamB) = (0, 0, 0, 0,0, 0,0,0,0,0);
my ($min_nodom) = (10);

GetOptions(
    "lav" => \$lav,
    "min_nodom=i" => \$min_nodom,
    "neg" => \$neg_doms,
    "neg_doms" => \$neg_doms,
    "neg-doms" => \$neg_doms,
    "no-over" => \$no_over,
    "no_over" => \$no_over,
    "pfamB" => \$pfamB,
    "pfacc" => \$pf_acc,
    "pfam_acc" => \$pf_acc,
    "acc" => \$pf_acc,
    "h|?" => \$shelp,
    "help" => \$help,
    );

pod2usage(1) if $shelp;
pod2usage(exitstatus => 0, verbose => 2) if $help;
pod2usage(1) unless @ARGV;

my %annot_types = ();
my %domains = (NODOM=>0);
my $domain_cnt = 0;

my $loc="http://pfam.xfam.org/";
my $url;

my @pf_domains;
my $pf_seq_length=0;

my $get_annot_sub = \&get_pfam_annots;

my $twig = new XML::Twig();

my ($tmp, $gi, $sdb, $acc, $id, $use_acc);

# get the query
my ($query, $seq_len) = @ARGV;
$seq_len = 0 unless defined($seq_len);

$query =~ s/^>//;

my $ANN_F;

my @annots = ();

#if it's a file I can open, read and parse it
if ($query !~ m/\|/ && open($ANN_F, $query)) {

  while (my $a_line = <$ANN_F>) {
    $a_line =~ s/^>//;
    chomp $a_line;
    push @annots, show_annots($a_line, $get_annot_sub);
  }
}
else {
  push @annots, show_annots("$query $seq_len", $get_annot_sub);
}

for my $seq_annot (@annots) {
  print ">",$seq_annot->{seq_info},"\n";
  for my $annot (@{$seq_annot->{list}}) {
    if (!$lav && defined($domains{$annot->[-1]})) {
      $annot->[-1] .= " :".$domains{$annot->[-1]};
    }
    print join("\t",@$annot),"\n";
  }
}

exit(0);

sub show_annots {
  my ($query_len, $get_annot_sub) = @_;

  my ($annot_line, $seq_len) = split(/\s+/,$query_len);

  my $pfamA_acc;

  my %annot_data = (seq_info=>$annot_line);

  $use_acc = 1;
  if ($annot_line =~ m/^pf26\|/) {
    ($sdb, $gi, $acc, $id) = split(/\|/,$annot_line);
  }
  elsif ($annot_line =~ m/^gi\|/) {
    ($tmp, $gi, $sdb, $acc, $id) = split(/\|/,$annot_line);
  }
  elsif ($annot_line =~ m/^sp\|/) {
    ($sdb, $acc, $id) = split(/\|/,$annot_line);
  }
  elsif ($annot_line =~ m/^tr\|/) {
    ($sdb, $acc, $id) = split(/\|/,$annot_line);
  }
  elsif ($annot_line =~ m/^SP:/i) {
    ($sdb, $id) = split(/:/,$annot_line);
    $use_acc = 0;
  }

  # remove version number
  unless ($use_acc) {
    $annot_data{list} = get_pfam_www($id);
  }
  else {
    $acc =~ s/\.\d+$//;
    $annot_data{list} = get_pfam_www($acc);
  }

  return \%annot_data;
}

sub get_length {
    my ($t, $elt) = @_;
    $pf_seq_length = $elt->{att}->{length};
}

sub push_match {
    my ($t, $elt) = @_;
#    return unless ($elt->{att}->{type} =~ m/Pfam-A/);
    my $attr_ref = $elt->{att};
    my $loc_ref = $elt->first_child('location')->{att};
    push @pf_domains, { %$attr_ref, %$loc_ref };
}

sub get_pfam_www {
  my ($acc, $seq_length) = @_;

#  if ($acc =~ m/_/) {$url = "protein?id=$acc&output=xml"; }
#  else {$url = "protein/$acc?output=xml"; }

  $url = "protein/$acc?output=xml";

  my $res = get($loc . $url);

  @pf_domains = ();

  my $twig = XML::Twig->new(twig_roots => {matches => 1, sequence => 1},
#			    start_tag_handlers => {
#						   'sequence' => \&get_length,
#						  },
			    twig_handlers => {
					      'match' => \&push_match,
					      'sequence' => \&get_length,
					     },
			    pretty_print => 'indented');
  my $xml = $twig->parse($res);

  $seq_length = $pf_seq_length;


  @pf_domains = sort { $a->{start} <=> $b->{start} } @pf_domains;

  unless ($pfamB) {
    @pf_domains = grep { $_->{type} !~ m/Pfam-B/ } @pf_domains;
  }

  # check for domain overlap, and resolve check for domain overlap
  # (possibly more than 2 domains), choosing the domain with the best
  # evalue

  for my $dom_ref (@pf_domains) {
    if ($pf_acc) {
      $dom_ref->{info} = $dom_ref->{accession};
    }
    else {
      $dom_ref->{info} = $dom_ref->{id};
    }
  }

  if($no_over && scalar(@pf_domains) > 1) {

    my @tmp_domains = @pf_domains;
    my @save_domains = ();

    my $prev_dom = shift @tmp_domains;

    while (my $curr_dom = shift @tmp_domains) {

      my @overlap_domains = ($prev_dom);

      my $diff = $prev_dom->{end} - $curr_dom->{start};
      # check for overlap > domain_length/3

      my ($prev_len, $cur_len) = ($prev_dom->{end}-$prev_dom->{start}+1, $curr_dom->{end}-$curr_dom->{start}+1);
      my $inclusion = ((($curr_dom->{start} >= $prev_dom->{start}) && ($curr_dom->{end} <= $prev_dom->{end})) ||
		       (($curr_dom->{start} <= $prev_dom->{start}) && ($curr_dom->{end} >= $prev_dom->{end})));

      my $longer_len = ($prev_len > $cur_len) ? $prev_len : $cur_len;

      while ($inclusion || ($diff > 0 && $diff > $longer_len/3)) {
	push @overlap_domains, $curr_dom;
	$curr_dom = shift @tmp_domains;
	last unless $curr_dom;
	$diff = $prev_dom->{end} - $curr_dom->{start};
	($prev_len, $cur_len) = ($prev_dom->{end}-$prev_dom->{start}+1, $curr_dom->{end}-$curr_dom->{start}+1);
	$longer_len = ($prev_len > $cur_len) ? $prev_len : $cur_len;
	$inclusion = ((($curr_dom->{start} >= $prev_dom->{start}) && ($curr_dom->{end} <= $prev_dom->{end})) ||
		      (($curr_dom->{start} <= $prev_dom->{start}) && ($curr_dom->{end} >= $prev_dom->{end})));
      }

      # check for overlapping domains; >1 because $prev_dom is always there
      if (scalar(@overlap_domains) > 1 ) {
	# if $rpd2_fams, check for a chosen one

	for my $dom ( @overlap_domains) {
	  $dom->{evalue} = 1.0 unless defined($dom->{evalue});
	}

	@overlap_domains = sort { $a->{evalue} <=> $b->{evalue} } @overlap_domains;
	$prev_dom = $overlap_domains[0];
      }

      # $prev_dom should be the best of the overlaps, and we are no longer overlapping > dom_length/3
      push @save_domains, $prev_dom;
      $prev_dom = $curr_dom;
    }
    if ($prev_dom) {push @save_domains, $prev_dom;}

    @pf_domains = @save_domains;

    # now check for smaller overlaps
    for (my $i=1; $i < scalar(@pf_domains); $i++) {
      if ($pf_domains[$i-1]->{end} >= $pf_domains[$i]->{start}) {
	my $overlap = $pf_domains[$i-1]->{end} - $pf_domains[$i]->{start};
	$pf_domains[$i-1]->{end} -= int($overlap/2);
	$pf_domains[$i]->{start} = $pf_domains[$i-1]->{end}+1;
      }
    }
  }

  if ($neg_doms) {
    my @npf_domains;
    my $prev_dom={end=>0};
    for my $curr_dom ( @pf_domains) {
      if ($curr_dom->{start} - $prev_dom->{end} > $min_nodom) {
	my %new_dom = (start=>$prev_dom->{end}+1, end => $curr_dom->{start}-1, info=>'NODOM');
	push @npf_domains, \%new_dom;
      }
      push @npf_domains, $curr_dom;
      $prev_dom = $curr_dom;
    }
    if ($seq_length - $prev_dom->{end} > $min_nodom) {
      my %new_dom = (start=>$prev_dom->{end}+1, end=>$seq_length, info=>'NODOM');
      if ($new_dom{end} > $new_dom{start}) {push @npf_domains, \%new_dom;}
    }

    # @npf_domains has both old @pf_domains and new neg-domains
    @pf_domains = @npf_domains;
  }

  # now make sure we have useful names: colors

  for my $pf (@pf_domains) {
    $pf->{info} = domain_name($pf->{info}, $acc );
  }

  my @feats = ();
  for my $d_ref (@pf_domains) {
    if ($lav) {
      push @feats, [$d_ref->{start}, $d_ref->{end}, $d_ref->{info}];
    }
    else {
      push @feats, [$d_ref->{start}, '-', $d_ref->{end},  $d_ref->{info} ];
#      push @feats, [$d_ref->{end}, ']', '-', ""];
    }
  }

  return \@feats;
}

# domain name takes a uniprot domain label, removes comments ( ;
# truncated) and numbers and returns a canonical form. Thus:
# Cortactin 6.
# Cortactin 7; truncated.
# becomes "Cortactin"
#

sub domain_name {

  my ($value, $seq_id) = @_;

  unless (defined($value)) {
    warn "missing domain name for $seq_id";
    return "";
  }

  if (!defined($domains{$value})) {
    $domain_cnt++;
    $domains{$value} = $domain_cnt;
  }
  return $value;
}

__END__

=pod

=head1 NAME

ann_feats.pl

=head1 SYNOPSIS

 ann_pfam_www_e.pl --neg-doms  'sp|P09488|GSTM1_NUMAN' | accession.file

=head1 OPTIONS

 -h	short help
 --help include description

 --lav  produce lav2plt.pl annotation format, only show domains/repeats
 --neg-doms : report domains between annotated domains as NODOM
                 (also --neg, --neg_doms)
 --no-over  : generate non-overlapping domains (equivalent to ann_pfam_www.pl)
 --min_nodom=10  : minimum length between domains for NODOM

=head1 DESCRIPTION

C<ann_pfam_www_e.pl> extracts domain information from the Pfam www site
(pfam.xfam.org).  Currently, the program works with database
sequence descriptions in several formats:

 >gi|1705556|sp|P54670.1|CAF1_DICDI
 >sp|P09488|GSTM1_HUMAN
 >sp:CALM_HUMAN 

C<ann_pfam_www_e.pl> uses the Pfam RESTful WWW interface
(C<pfam.xfam.org/help#tabview=10>) to download domain
names/locations/score. C<ann_pfam_www_e.pl> is an alternative to
C<ann_pfam_e.pl> that does not require a MySQL instance with a Pfam
database installation.

If the "--no-over" option is set, overlapping domains are selected and
edited to remove overlaps.  For proteins with multiple overlapping
domains (domains overlap by more than 1/3 of the domain length),
C<auto_pfam_e.pl> selects the domain annotation with the best
C<domain_evalue_score>.  When domains overlap by less than 1/3 of the
domain length, they are shortened to remove the overlap.

C<ann_pfam_www_e.pl> is designed to be used by the B<FASTA> programs
with the C<-V \!ann_pfam_www_e.pl> or C<-V "\!ann_pfam_www_e.pl --neg">
option.

=head1 AUTHOR

William R. Pearson, wrp@virginia.edu

=cut
