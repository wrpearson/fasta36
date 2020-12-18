#!/usr/bin/env perl

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

# ann_upfeats_pfam_www.pl gets an annotation file from fasta36 -V with a line of the form:

# SP:GSTM1_HUMAN P09488 218
#
# it must:
# (1) read in the line
# (2) parse it to get the up_acc
# (3) return the tab delimited features
#

# this version can read feature2 uniprot features (acc/pos/end/label/value), but returns sorted start/end domains

use warnings;
use strict;

use Getopt::Long;
use Pod::Usage;
use LWP::Simple;
use XML::Twig;
use JSON qw(decode_json);
## use IO::String;

my $up_base = 'http://www.ebi.ac.uk/proteins/api/features';

my %domains = ();
my $domain_cnt = 0;

my ($lav, $neg_doms, $no_doms, $no_feats, $no_over, $shelp, $help, $no_vars) = (0,0,0,0,0,0,0,0);
my ($auto_reg, $vdoms, $no_clans, $pf_acc_flag, $acc_comment, $bound_comment) =  (0, 0, 0, 0, 0, 0);

my $color_sep_str = " :";
$color_sep_str = '~';

my ($min_nodom, $min_vdom) = (10,10);

GetOptions(
    "lav" => \$lav,
    "acc_comment" => \$acc_comment,
    "bound_comment" => \$bound_comment,
    "no-over" => \$no_over,
    "no_doms|no-doms|nodoms" => \$no_doms,
    "neg" => \$neg_doms,
    "neg_doms|neg-doms|negdoms" => \$neg_doms,
    "min_nodom=i" => \$min_nodom,
    "no_feats|no-feats|nofeats" => \$no_feats,
    "no-vars|no_vars" => \$no_vars,
    "vdoms" => \$vdoms,
    "v_doms" => \$vdoms,
    "pfacc" => \$pf_acc_flag,
    "pfam_acc" => \$pf_acc_flag,
    "acc" => \$pf_acc_flag,
    "h|?" => \$shelp,
    "help" => \$help,
    );

pod2usage(1) if $shelp;
pod2usage(exitstatus => 0, verbose => 2) if $help;
pod2usage(1) unless (@ARGV || -f STDIN || -p STDIN);

my @feat_keys = qw( ACT_SITE MOD_RES BINDING METAL SITE );
my @feat_vals = ( '=','*','#','^','@');
my @feat_names = ('Active site', 'Modified', 'Binding', 'Metal binding', 'Site');

unless ($no_vars) {
  push @feat_keys, qw(MUTAGEN VARIANT);
  push @feat_vals, ('V','V',);
  push @feat_names, ("","",);
}

my %feats_text = ();
@feats_text{@feat_keys} = ('Active site', '', 'Substrate binding', 'Metal binding', 'Site', '','');

my %feats_label;
@feats_label{@feat_keys} = ('Active site', 'Modified', 'Substrate binding', 'Metal binding', 'Site', '','');

# my @feat_vals = ( '=','*','#','^','@','V','V');


my @dom_keys = qw( DOMAIN REPEAT );
my @dom_vals = ( [ '[', ']'],[ '[', ']']);

my @ssr_keys = qw(STRAND HELIX);
my @ssr_vals = ( [ '[', ']']);

my %annot_types = ();

# from ann_pfam_www_e.pl 
my %domain_clan = (NODOM => {clan_id => 'NODOM', clan_acc=>0, domain_cnt=>0});
my @domain_list = (0);

my $loc="http://pfam.xfam.org/";
my $pf_url;

my @pf_domains;
my %pfamA_fams = ();
my ($pf_seq_length, $pf_model_length)=(0,0);
my ($clan_acc, $clan_id) = ("","");

my $get_domain_sub = \&get_pfam_annots;

if ($lav) {
  $no_feats = 1;
}

@annot_types{@feat_keys} = @feat_vals unless ($no_feats);

if ($neg_doms) {
  $domains{'NODOM'}=0;
}

my ($tmp, $gi, $sdb, $acc, $id, $use_acc);

unless ($no_feats) {
  for my $i ( 0 .. $#feat_keys) {
    next unless $feats_label{$feat_keys[$i]};
    print "=",$feat_vals[$i],":",$feats_label{$feat_keys[$i]},"\n";
  }
}

# get the query
my ($query, $seq_len) =  @ARGV;
$seq_len = 0 unless defined($seq_len);

$query =~ s/^>// if ($query);

my @annots = ();

#if it's a file I can open, read and parse it

# unless ($query && ($query =~ m/[\|:]/
#		   || $query =~ m/^[OPQ][0-9][A-Z0-9]{3}[0-9]|[A-NR-Z][0-9]([A-Z][A-Z0-9]{2}[0-9]){1,2}\s/
#		   || $query =~ m/^(XN)(MP)_\d+/)) {
if (! $query || -r $query) {
    while (my $a_line = <>) {
	$a_line =~ s/^>//;
	chomp $a_line;
	push @annots, upfeats_pfam_www($a_line, \&up_json_annots, \&get_pfam_www);
    }
} else {
  push @annots, upfeats_pfam_www("$query\t$seq_len", \&up_json_annots, \&get_pfam_www);
}

for my $seq_annot (@annots) {
  print ">",$seq_annot->{seq_info},"\n";
  for my $annot (@{$seq_annot->{list}}) {
    if (!$lav && defined($domains{$annot->[-1]})) {
      my ($a_name, $a_num) = domain_num($annot->[-1],$domains{$annot->[-1]});
      $annot->[-1] = $a_name;
      if ($acc_comment) {
	$annot->[-1] .= "{$domain_list[$a_num]}";
      }
      if ($bound_comment) {
	$annot->[-1] .= $color_sep_str.$annot->[0].":".$annot->[2];
      }
      $annot->[-1] .= $color_sep_str.$a_num;
    }
    print join("\t",@$annot),"\n";
  }
}

exit(0);

sub upfeats_pfam_www {
  my ($query_len, $get_upfeats_sub, $get_pfam_sub) = @_;

  my ($annot_line, $seq_len) = split(/\t/,$query_len);

  my %annot_data = (seq_info=>$annot_line);

  $use_acc = 1;

  if ($annot_line =~ m/^gi\|/) {
    ($tmp, $gi, $sdb, $acc, $id) = split(/\|/,$annot_line);
  } elsif ($annot_line =~ m/^(SP|TR):(\w+)\s(\w+)/) {
    ($sdb, $id, $acc) = (lc($1), $2, $3);
    $use_acc = 1;
  } elsif ($annot_line =~ m/^(SP|TR):(\w+)/) {
    ($sdb, $id) = (lc($1), $2, '');
    warn("*** $0 - accession required: $annot_line");
    $use_acc = 0;
  } elsif ($annot_line =~ m/^(UR\d{3}:UniRef\d{2})_(\w+)/) {
    $sdb = lc($1);
    $id = $2;
#    $acc = $2;
  } elsif ($annot_line =~ m/\|/) {
    ($sdb, $acc, $id) = split(/\|/,$annot_line);
  }
  else {
    ($acc) = split(/\s+/,$annot_line);
  }

  $acc =~ s/\.\d+// if ($acc);
  $annot_data{list} = [];
  my $lwp_features = "";

  if ($acc && ($acc =~ m/^[A-Z][0-9][A-Z0-9]{3}[0-9]/)) {
    $lwp_features = get("$up_base/$acc.json");
  }

  my @annots = ();

  if ($lwp_features && ($lwp_features !~ /ERROR/)) {
    push @annots, up_json_annots(\%annot_types, $lwp_features, $seq_len);
  }

  unless ($use_acc) {
    push @annots, $get_pfam_sub->($id, $seq_len);
  }
  else {
    $acc =~ s/\.\d+$//;
    push @annots, $get_pfam_sub->($acc, $seq_len);
  }

  @annots  = sort { $a->[0] <=> $b->[0] } @annots;
  $annot_data{list} = \@annots;

  return \%annot_data;
}

# parses www.uniprot.org gff feature table
sub up_json_annots {
  my ($annot_types, $annot_data, $seq_len) = @_;

  my $json_ref = decode_json($annot_data);

  my ($acc, $pos, $end, $label, $value, $comment, $len);

  $seq_len = 0;

  my @sites = ();  # sites with one position

  my ($seq_str, $seq_acc, $seq_id)  = @{$json_ref}{qw(sequence accession entryName)};

  for my $feat ( @{$json_ref->{features}} ) {
    next unless ($annot_types->{$feat->{type}});

    my ($label, $pos, $end, $value)  = @{$feat}{qw(type begin end description)};

    if ($label =~ m/VARIANT/ || $label =~ m/MUTAGEN/) {
      push @sites, [$pos, $annot_types->{$label}, $feat->{alternativeSequence}, $value];
    }
    else {
      next unless ($pos == $end);
      if ($feats_text{$label}) {
	my $info = $feats_text{$label};
	if ($value) {
	  $info .= ": $value";
	}
	push @sites, [$pos, $annot_types->{$label}, "-", $info];
      } else {
	push @sites, [$pos, $annot_types->{$label}, "-", $value];
      }
    }
  }

  return @sites;
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

sub get_model_length {
    my ($t, $elt) = @_;
    $pf_model_length = $elt->{att}->{model_length};
}

sub get_clan {
    my ($t, $elt) = @_;
    my $attr_ref = $elt->{att};
#    print Dumper($attr_ref);
    ($clan_acc, $clan_id) = ($attr_ref->{clan_acc},$attr_ref->{clan_id});
}

sub get_pfam_www {
  my ($acc, $seq_length) = @_;

#  if ($acc =~ m/_/) {$url = "protein?id=$acc&output=xml"; }
#  else {$url = "protein/$acc?output=xml"; }

  $pf_url = "protein/$acc?output=xml";

  my $res = get($loc . $pf_url);

  @pf_domains = ();

  my $twig_dom = XML::Twig->new(twig_roots => {matches => 1, sequence => 1},
#			    start_tag_handlers => {
#						   'sequence' => \&get_length,
#						  },
			    twig_handlers => {
					      'match' => \&push_match,
					      'sequence' => \&get_length,
					     },
			    pretty_print => 'indented');
  my $xml = $twig_dom->parse($res);

  if (!$seq_length || $seq_length == 0) {
      $seq_length = $pf_seq_length;
  }

  @pf_domains = sort { $a->{start} <=> $b->{start} } @pf_domains;

  # to look for possible joining, need model_length
  for my $curr_dom (@pf_domains) {
      
      my $acc = $curr_dom->{accession};
      $pf_url = "family/$acc?output=xml";

      my $res = get($loc . $pf_url);

      my $twig_fam = XML::Twig->new(twig_roots => {hmm_details => 1, clan_membership=> 1},
				    twig_handlers => {
					'hmm_details' => \&get_model_length,
					'clan_membership' => \&get_clan,
				    },
				    pretty_print => 'indented');

      ($clan_acc, $clan_id) = ("","");
      my $fam_xml = $twig_fam->parse($res);

      $pfamA_fams{$acc} = { model_length => $pf_model_length, clan_acc=>$clan_acc, clan_id=>$clan_id};
      $curr_dom->{model_length} = $pf_model_length;
  }

  # check for domain overlap, and resolve check for domain overlap
  # (possibly more than 2 domains), choosing the domain with the best
  # evalue

  my @raw_pf_domains = @pf_domains;
  @pf_domains = ();

  for my $dom_ref (@raw_pf_domains) {
    if ($pf_acc_flag) {
      $dom_ref->{info} = $dom_ref->{accession};
    }
    else {
      $dom_ref->{info} = $dom_ref->{id};
    }
    next if ($dom_ref->{start} >= $seq_length);
    if ($dom_ref->{end} >= $seq_length) {
	$dom_ref->{end} = $seq_length;
    }
    push @pf_domains, $dom_ref;
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

  # before checking for domain overlap, check for "split-domains"
  # (self-unbound) by looking for runs of the same domain that are
  # ordered by model_start

  if (scalar(@pf_domains) > 1) {
    my @j_domains;		#joined domains
    my @tmp_domains = @pf_domains;

    my $prev_dom = shift(@tmp_domains);

    for my $curr_dom (@tmp_domains) {
      # to join domains:
      # (1) the domains must be in order by model_start/end coordinates
      # (3) joining the domains cannot make the total combination too long

      # check for model and sequence consistency
      if ($prev_dom->{accession} eq $curr_dom->{accession}) { # same family

	my $prev_dom_len = $prev_dom->{hmm_end}-$prev_dom->{hmm_start}+1;
	my $curr_dom_len = $curr_dom->{hmm_end}-$curr_dom->{hmm_start}+1;
	my $prev_dom_fn = $prev_dom_len/$curr_dom->{model_length};
	my $curr_dom_fn = $curr_dom_len/$curr_dom->{model_length};
	my $missing_dom_fn = max(0,($curr_dom->{hmm_start} - $prev_dom->{hmm_end})/$curr_dom->{model_length});

	if ( $prev_dom->{hmm_start} < $curr_dom->{hmm_start} # model check
	     && $prev_dom->{hmm_end} < $curr_dom->{hmm_end}
	     && ($curr_dom->{hmm_start} > $prev_dom->{hmm_end} * 0.80 # limit overlap
		 || $curr_dom->{hmm_start} <  $prev_dom->{hmm_end} * 1.25)
	     && $prev_dom_fn + $curr_dom_fn < 1.33
	     && $missing_dom_fn < min($prev_dom_fn,$curr_dom_fn)) { # join them by updating $prev_dom
	  $prev_dom->{end} = $curr_dom->{end};
	  $prev_dom->{hmm_end} = $curr_dom->{hmm_end};
	  $prev_dom->{evalue} = ($prev_dom->{evalue} < $curr_dom->{evalue} ? $prev_dom->{evalue} : $curr_dom->{evalue});
	  next;
	}

	push @j_domains, $prev_dom;
	$prev_dom = $curr_dom;
      }
      else {
	  push @j_domains, $prev_dom;
	  $prev_dom = $curr_dom;
      }
    }
    push @j_domains, $prev_dom;
    @pf_domains = @j_domains;
  }

  # $vdoms -- virtual Pfam domains -- the equivalent of $neg_doms,
  # but covering parts of a Pfam model that are not annotated.  split
  # domains have been joined, so simply check beginning and end of
  # each domain (but must also check for bounded-ness)
  # only add when 10% or more is missing and missing length > $min_nodom

  if ($vdoms) {
      my @vpf_domains;

      my $curr_dom = $pf_domains[0];

      my $prev_dom={end=>0, accession=>''};
      my $prev_dom_end = 0;
      my $next_dom_start = $seq_length+1;

      for (my $dom_ix=0; $dom_ix < scalar(@pf_domains); $dom_ix++ ) {
	  $curr_dom = $pf_domains[$dom_ix];

	  my $pfamA =  $curr_dom->{accession};

	  # first, look left, is there a domain there (if there is,
	  # it should be updated right

	  # my $min_vdom = $curr_dom->{model_length} / 10;

	  if ($prev_dom->{accession}) { # look for previous domain
	      $prev_dom_end = $prev_dom->{end};
	  }

	  # there is a domain to the left, how much room is available?
	  my $left_dom_len = min($curr_dom->{start}-$prev_dom_end-1, $curr_dom->{hmm_start}-1);
	  if ( $left_dom_len > $min_vdom) {
	      # there is room for a virtual domain
	      my %new_dom = (start=> $curr_dom->{start}-$left_dom_len,
			     end => $curr_dom->{start}-1,
			     info=>'@'.$curr_dom->{accession},
			     model_length=> $curr_dom->{model_length},
			     hmm_end => $curr_dom->{hmm_start}-1,
			     hmm_start => $left_dom_len,
			     accession=>$pfamA,
		  );
	      push @vpf_domains, \%new_dom;
	  }

	  # save the current domain
	  push @vpf_domains, $curr_dom;
	  $prev_dom = $curr_dom;

	  if ($dom_ix < $#pf_domains) { # there is a domain to the right
	      # first, give all the extra space to the first domain (no splitting)
	      $next_dom_start = $pf_domains[$dom_ix+1]->{start};
	  }
	  else {
	      $next_dom_start = $seq_length;
	  }

	  # is there room for a virtual domain right

	  my $right_dom_len = min($next_dom_start-$curr_dom->{end}-1, # space available 
				  $curr_dom->{model_length}-$curr_dom->{hmm_end} # space needed
	      );
	  if ( $right_dom_len > $min_vdom) {
	      my %new_dom = (start=> $curr_dom->{end}+1,
			     end=> $curr_dom->{end}+$right_dom_len,
			     info=>'@'.$pfamA,
			     model_length => $curr_dom->{model_length},
			     accession=> $pfamA,
		  );
	      push @vpf_domains, \%new_dom;
	      $prev_dom = \%new_dom;
	  }
      }				# all done, check for last one
    # @vpf_domains has both old @pf_domains and new virtural-domains
    @pf_domains = @vpf_domains;
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
    $pf->{info} = domain_name($pf->{info}, $acc, $pf->{accession});
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

  return @feats;
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
    if (!defined($pfamA_fams{$pf_acc})) {

      my $pf_url = "family/$value?output=xml";

      my $res = get($loc . $pf_url);

      my $twig_clan = XML::Twig->new(twig_roots => {'clan_membership'=>1},
				     twig_handlers => {
						       'clan_membership' => \&get_clan,
						      },
				     pretty_print => 'indented');

      # make certain to reinitialize
      ($clan_acc, $clan_id) = ("","");
      my $xml = $twig_clan->parse($res);
    }
    else {
      ($clan_acc, $clan_id) = @{$pfamA_fams{$pf_acc}}{qw(clan_acc clan_id)};
    }

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

ann_feats_up_www2.pl

=head1 SYNOPSIS

 ann_feats_up_www2.pl --no_doms --no_feats --lav 'sp|P09488|GSTM1_NUMAN' | accession.file

=head1 OPTIONS

 -h	short help
 --help include description
 --no-doms  do not show domain boundaries (domains are always shown with --lav)
 --no-feats do not show feature (variants, active sites, phospho-sites)
 --lav  produce lav2plt.pl annotation format, only show domains/repeats

 --neg-doms,  -- report domains between annotated domains as NODOM
                 (also --neg, --neg_doms)
 --min_nodom=10  -- minimum length between domains for NODOM

 --host, --user, --password, --port --db -- info for mysql database

=head1 DESCRIPTION

C<ann_feats_up_www2.pl> extracts feature, domain, and repeat
information from the Uniprot DAS server through an XSLT transation
provided by http://www.ebi.ac.uk/Tools/dbfetch/dbfetch/uniprotkb.
This server provides GFF descriptions of Uniprot entries, with most of
the information provided in UniProt feature tables.

C<ann_feats_up_www2.pl> is an alternative to C<ann_pfam.pl> and
C<ann_pfam.pl> that does not require a local MySQL copy of Pfam.

Given a command line argument that contains a sequence accession
(P09488), the program looks up the features available for that
sequence and returns them in a tab-delimited format:

>sp|P09488|GSTM1_HUMAN
2	[	-	GST N-terminal :1
7	V	F	Mutagen: Reduces catalytic activity 100- fold.
23	*	-	MOD_RES: Phosphotyrosine (By similarity).
33	*	-	MOD_RES: Phosphotyrosine (By similarity).
34	*	-	MOD_RES: Phosphothreonine (By similarity).
88	]	-	
90	[	-	GST C-terminal :2
108	V	Q	Mutagen: Reduces catalytic activity by half.
108	V	S	Mutagen: Changes the properties of the enzyme toward some substrates.
109	V	I	Mutagen: Reduces catalytic activity by half.
116	#	-	BINDING: Substrate.
116	V	A	Mutagen: Reduces catalytic activity 10-fold.
116	V	F	Mutagen: Slight increase of catalytic activity.
173	V	N	in allele GSTM1B; dbSNP:rs1065411.
208	]	-	
210	V	T	in dbSNP:rs449856.

If features are provided, then a legend of feature symbols is provided
as well:

 =*:phosphorylation
 ==:active site
 =@:site
 =^:binding
 =!:metal binding

If the C<--lav> option is specified, domain and repeat features are
presented in a different format for the C<lav2plt.pl> program:

  >sp|P09488|GSTM1_HUMAN
  2	88	GST N-terminal.
  90	208	GST C-terminal.

C<ann_feats_up_www2.pl> is designed to be used by the B<FASTA> programs with
the C<-V \!ann_feats_up_www2.pl> option.  It can also be used with the lav2plt.pl
program with the C<--xA "\!ann_feats_up_www2.pl --lav"> or C<--yA "\!ann_feats_up_www2.pl --lav"> options.

=head1 AUTHOR

William R. Pearson, wrp@virginia.edu

=cut
