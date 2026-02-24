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

## updated 15-April-2024 to use interpro pfam  --- now ann_pfam_www2.pl

## updated 8-Nov-2022 to use pfam-legacy.xfam.org, since pfam has been discontinued

# ann_pfam_www.pl gets an annotation file from fasta36 -V with a line of the form:

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

use warnings;
# use strict;

use Getopt::Long;
use Pod::Usage;
use LWP::Simple;
use LWP::UserAgent;
use JSON qw(decode_json);

# use XML::Twig;
# use Data::Dumper;

my ($auto_reg,$rpd2_fams, $neg_doms, $vdoms, $lav, $no_clans, $pf_acc_flag, $shelp, $help, $no_over, $acc_comment, $bound_comment, $pfamB) =
  (0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0);
my ($show_color) = (1);
my ($min_nodom, $min_vdom) = (10, 10);

my $color_sep_str = " :";
$color_sep_str = '~';

GetOptions(
    "lav" => \$lav,
    "acc_comment" => \$acc_comment,
    "bound_comment" => \$bound_comment,
    "min_nodom=i" => \$min_nodom,
    "neg" => \$neg_doms,
    "neg_doms" => \$neg_doms,
    "neg-doms" => \$neg_doms,
    "no-over" => \$no_over,
    "no_over" => \$no_over,
    "no-clans" => \$no_clans,
    "no_clans" => \$no_clans,
    "color!" => \$show_color,
    "pfamB" => \$pfamB,
    "pfacc" => \$pf_acc_flag,
    "pfam_acc" => \$pf_acc_flag,
    "acc" => \$pf_acc_flag,
    "h|?" => \$shelp,
    "help" => \$help,
    );

pod2usage(1) if $shelp;
pod2usage(exitstatus => 0, verbose => 2) if $help;
pod2usage(1) unless (@ARGV || -p STDIN || -f STDIN);

my %annot_types = ();
my %domains = (NODOM=>0);
my %domain_clan = (NODOM => {clan_id => 'NODOM', clan_acc=>0, domain_cnt=>0});
my @domain_list = (0);
my $domain_cnt = 0;

## my $loc="https://pfam-legacy.xfam.org/";
my $ua = LWP::UserAgent->new(ssl_opts=>{verify_hostname => 0});
my $loc="https://www.ebi.ac.uk/interpro/api/";
my $url;

my @pf_domains = ();
my %pfamA_fams = ('NODOM'=>'NODOM');
my ($pf_seq_length, $pf_model_length)=(0,0);
my ($clan_acc, $clan_id) = ("","");

my $get_annot_sub = \&get_pfam_www;

my ($tmp, $gi, $sdb, $acc, $id, $use_acc);

# get the query
my ($query, $seq_len) = @ARGV;
$seq_len = 0 unless defined($seq_len);

$query =~ s/^>// if ($query);

my @annots = ();

#if it's a file I can open, read and parse it
# unless ($query && ($query =~ m/[\|:]/ ||
# 		   $query =~ m/^[OPQ][0-9][A-Z0-9]{3}[0-9]|[A-NR-Z][0-9]([A-Z][A-Z0-9]{2}[0-9]){1,2}\s/)) {
if (! $query || -r $query) {
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
    if (!$lav && defined($domains{$annot->[-1]})) {
      my ($a_name, $a_num) = domain_num($annot->[-1],$domains{$annot->[-1]});
      $annot->[-1] = $a_name;
      my $tmp_a_num = $a_num;
      $tmp_a_num =~ s/v$//;
      if ($acc_comment) {
	$annot->[-1] .= "{$domain_list[$tmp_a_num]}";
      }
      if ($bound_comment) {
	$annot->[-1] .= $color_sep_str.$annot->[0].":".$annot->[2];
      }
      elsif ($show_color) {
	$annot->[-1] .= $color_sep_str.$a_num;
      }
    }
    print join("\t",@$annot),"\n";
  }
}

exit(0);

sub show_annots {
  my ($query_len, $get_annot_sub) = @_;

  my ($annot_line, $seq_len) = split(/\t/,$query_len);

  my $pfamA_acc;

  my %annot_data = (seq_info=>$annot_line);

  $use_acc = 1;
  if ($annot_line =~ m/^pf26\|/) {
    ($sdb, $gi, $acc, $id) = split(/\|/,$annot_line);
  }
  elsif ($annot_line =~ m/^gi\|/) {
    ($tmp, $gi, $sdb, $acc, $id) = split(/\|/,$annot_line);
  }
  elsif ($annot_line =~ m/^(sp|tr|up)\|/) {
    ($sdb, $acc, $id) = split(/\|/,$annot_line);
    $use_acc = 1;
  }
  elsif ($annot_line =~ m/^(SP|TR):(\w+) (\w+)/) {
    ($sdb, $id, $acc) = (lc($1), $2, $3);
    $use_acc = 1;
  }
  elsif ($annot_line =~ m/^(SP|TR):(\w+)/) {
    ($sdb, $id, $acc) = (lc($1), $2, "");
    $use_acc = 0;
  }
  elsif ($annot_line !~ m/\|/ && $annot_line !~ m/:/) {
    $use_acc = 1;
    ($acc) = split(/\s+/,$annot_line);
  }

  # remove version number
  unless ($use_acc) {
    $annot_data{list} = get_pfam_www($id, $seq_len);
  }
  else {
    $acc =~ s/\.\d+$//;
    $annot_data{list} = get_pfam_www($acc, $seq_len);
  }

  return \%annot_data;
}

sub get_pfam_id_www {
    my ($pf_acc) = @_;

    if ($pf_acc eq 'NODOM') {
	return 'NODOM';
    }

    my $url = "entry/pfam/$pf_acc";

    my $res = get_https($loc . $url);

    $pfam_info= decode_json($res);

    return $pfam_info->{'metadata'}{'name'}{'short'};
}


sub get_pfam_www {
  my ($acc, $seq_length) = @_;

  $url = "entry/pfam/protein/uniprot/$acc";

  my $res = get_https($loc . $url);

  @pf_domains = ();

  if (! $res) {
      return \@pf_domains;
  }

  my $json_info = decode_json($res);

  if (!defined($seq_length) || $seq_length < 1) {
      $seq_length = $json_info->{'results'}[0]{'proteins'}[0]{'protein_length'};
  }

  for my $result ( @{$json_info->{'results'}} ) {
      my $pfam_info = $result->{'metadata'};
      my $pfam_acc = $pfam_info->{'accession'};
      my $pfam_name = $pfam_info->{'name'};

      for my $protein ( @{$result->{'proteins'}} ) {
	  for my $entry ( @{$protein->{'entry_protein_locations'}} ) {
	      for my $frag ( @{$entry->{'fragments'}}) {
		  push @pf_domains, {'pf_acc'=>$pfam_acc, 'start'=>$frag->{'start'}, 'end'=>$frag->{'end'}, 'name'=>$pfam_name};
	      }
	  }
      }
  }

  @pf_domains = sort { $a->{start} <=> $b->{start} } @pf_domains;

  ## put in short name, and clan info, if not seen before

  for my $pf_dom (@pf_domains) {
      my $pf_acc = $pf_dom->{pf_acc};
      if (! exists($pfamA_fams{$pf_acc})) {
	  $pfamA_fams{$pf_acc} = get_pfam_id_www($pf_acc);
      }
      if (! $pf_acc_flag) {
	  $pf_dom->{info} = $pfamA_fams{$pf_acc};
      }
      else {
	  $pf_dom->{info} = $pf_acc;
      }

  }


  # check for domain overlap, and resolve check for domain overlap
  # (possibly more than 2 domains), choosing the domain with the best
  # evalue

  my @raw_pf_domains = @pf_domains;
  @pf_domains = ();

  for my $dom_ref (@raw_pf_domains) {
    next if ($dom_ref->{start} >= $seq_length);
    if ($dom_ref->{end} >= $seq_length) {
	$dom_ref->{end} = $seq_length;
    }
    push @pf_domains, $dom_ref;
  }

  if ($no_over && scalar(@pf_domains) > 1) {

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

  # before checking for domain overlap, check for "split-domains"
  # (self-unbound) by looking for runs of the same domain that are
  # ordered by model_start

  ## cannot find split-domains without model coordinates -- removed

  # $vdoms -- virtual Pfam domains -- the equivalent of $neg_doms,
  # but covering parts of a Pfam model that are not annotated.  split
  # domains have been joined, so simply check beginning and end of
  # each domain (but must also check for bounded-ness)
  # only add when 10% or more is missing and missing length > $min_nodom

  ## vdoms are not available with the InterPro Pfam implementation,
  ## because model coordinates are not provided

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
    $pf->{info} = domain_name($pf->{info}, $acc, $pf->{pf_acc});
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

# in addition, domain_name() looks up each domain name to see if it
# has a clan, and, if different domains share the same clan, they get
# the same colors.

sub domain_name {

  my ($value, $seq_id, $pf_acc) = @_;
  my $is_virtual = 0;

  if ($value =~ m/^@/) {
    $is_virtual = 1;
    $value =~ s/^@//;
  }

  unless (defined($value)) {
    warn "missing domain name for $seq_id";
    return "";
  }

  if ($no_clans) {
    if (! defined($domains{$value})) {
      $domain_clan{$value} = 0;
      $domains{$value} = ++$domain_cnt;
      push @domain_list, $pf_acc;
    }
  }
  elsif (!defined($domain_clan{$value})) {
    ## only do this for new domains, old domains have known mappings

    ## ways to highlight the same domain:
    # (1) for clans, substitute clan name for family name
    # (2) for clans, use the same color for the same clan, but don't change the name
    # (3) for clans, combine family name with clan name, but use colors based on clan

    # return the clan name, identifier if a clan member
    if (!defined($domain_clan{$pf_acc})) {

      my $url = "set/pfam/entry/pfam/$pf_acc";
      my $res = get_https($loc . $url);

      my $clan_info = decode_json($res);

      if (exists($clan_info->{results}[0]{metadata})) {
	  ($clan_acc, $clan_id) = @{$clan_info->{results}[0]{metadata}}{qw(accession name)};
	  $domain_clan{$pf_acc} = { clan_acc=>$clan_acc, clan_id=>$clan_id};
      }
      else {
	  $domain_clan{$pf_acc} = { clan_acc=>"", clan_id=>""};
      }
    }

    ($clan_acc, $clan_id) = @{$domain_clan{$pf_acc}}{qw(clan_acc clan_id)};

    if ($clan_acc) {
      my $c_value = "C." . $clan_id;
      if ($pf_acc_flag) {$c_value = "C." . $clan_acc;}

      $domain_clan{$value} = {clan_id => $clan_id,
			      clan_acc => $clan_acc};

      if ($domains{$c_value}) {
	$domain_clan{$value}->{domain_cnt} =  $domains{$c_value};
	$value = $c_value;
      }
      else {
	$domain_clan{$value}->{domain_cnt} = ++ $domain_cnt;
	$value = $c_value;
	$domains{$value} = $domain_cnt;
	push @domain_list, $pf_acc;
      }
    }
    else {
      $domain_clan{$value} = 0;
      $domains{$value} = ++$domain_cnt;
      push @domain_list, $pf_acc;
    }
  }
  elsif ($domain_clan{$value} && $domain_clan{$value}->{clan_acc}) {
    if ($pf_acc_flag) {$value = "C." . $domain_clan{$value}->{clan_acc};}
    else { $value = "C." . $domain_clan{$value}->{clan_id}; }
  }

  if ($is_virtual) {
    $domains{'@'.$value} = $domains{$value};
    $value = '@'.$value;
  }
  return $value;
}

sub domain_num {
  my ($value, $number) = @_;
  if ($value =~ m/^@/) {
    $value =~ s/^@/v/;
    $number = $number."v";
  }
  return ($value, $number);
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

sub min {
  my ($arg1, $arg2) = @_;

  return ($arg1 <= $arg2 ? $arg1 : $arg2);
}

sub max {
  my ($arg1, $arg2) = @_;

  return ($arg1 >= $arg2 ? $arg1 : $arg2);
}

__END__

=pod

=head1 NAME

ann_feats.pl

=head1 SYNOPSIS

 ann_pfam_www2.pl --neg-doms  'sp|P09488|GSTM1_NUMAN' | accession.file

=head1 OPTIONS

 -h	short help
 --help include description

 --lav  produce lav2plt.pl annotation format, only show domains/repeats
 --neg-doms : report domains between annotated domains as NODOM
                 (also --neg, --neg_doms)
 --no-over  : generate non-overlapping domains
 --no-clans : do not use clans with multiple families from same clan
 --pfam_acc : report Pfam accession
 --min_nodom=10  : minimum length between domains for NODOM

=head1 DESCRIPTION

C<ann_pfam_www2.pl> extracts domain information from the Pfam www site
(pfam.xfam.org).  Currently, the program works with database
sequence descriptions in several formats:

 >gi|1705556|sp|P54670.1|CAF1_DICDI
 >sp|P09488|GSTM1_HUMAN
 >sp:CALM_HUMAN 

C<ann_pfam_www2.pl> uses the Interpro/Pfam RESTful WWW interface
(C<pfam-docs.readthedocs.io/en/latest/api.html>) to download domain
names/locations/score. C<ann_pfam_www2.pl> is an alternative to
C<ann_pfam_sql.pl> that does not require a MySQL instance with a Pfam
database installation.

If the "--no-over" option is set, overlapping domains are selected and
edited to remove overlaps.  For proteins with multiple overlapping
domains (domains overlap by more than 1/3 of the domain length),
C<auto_pfam_e.pl> selects the domain annotation with the best
C<domain_evalue_score>.  When domains overlap by less than 1/3 of the
domain length, they are shortened to remove the overlap.

C<ann_pfam_www2.pl> is designed to be used by the B<FASTA> programs
with the C<-V \!ann_pfam_www_e.pl> or C<-V "\!ann_pfam_www_e.pl --neg">
option.

=head1 AUTHOR

William R. Pearson, wrp@virginia.edu

=cut
