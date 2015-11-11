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

# ann_pfam_e.pl gets an annotation file from fasta36 -V with a line of the form:

# gi|62822551|sp|P00502|GSTA1_RAT Glutathione S-transfer\n  (at least from pir1.lseg)
#
# it must:
# (1) read in the line
# (2) parse it to get the up_acc
# (3) return the tab delimited features
#

# this version only annotates sequences known to Pfam:pfamseq:
# >pf26|164|O57809|1A1D_PYRHO
# and only provides domain information

use strict;

use DBI;
use Getopt::Long;
use Pod::Usage;

use vars qw($host $db $port $user $pass);

my $hostname = `/bin/hostname`;

($host, $db, $port, $user, $pass)  = ("wrpxdb.its.virginia.edu", "pfam27", 0, "web_user", "fasta_www");
#$host = 'xdb';

my ($auto_reg,$rpd2_fams, $neg_doms, $lav, $no_doms, $no_clans, $pf_acc, $no_over, $acc_comment, $shelp, $help) = 
  (0, 0, 0, 0,0, 0, 0, 0, 0, 0, 0);
my ($min_nodom) = (10);

my $color_sep_str = " :";
$color_sep_str = '~';

GetOptions(
    "host=s" => \$host,
    "db=s" => \$db,
    "user=s" => \$user,
    "password=s" => \$pass,
    "port=i" => \$port,
    "lav" => \$lav,
    "acc_comment" => \$acc_comment,
    "no-over" => \$no_over,
    "no_over" => \$no_over,
    "no-clans" => \$no_clans,
    "no_clans" => \$no_clans,
    "neg" => \$neg_doms,
    "neg_doms" => \$neg_doms,
    "neg-doms" => \$neg_doms,
    "min_nodom=i" => \$min_nodom,
    "pfacc" => \$pf_acc,
    "RPD2" => \$rpd2_fams,
    "auto_reg" => \$auto_reg,
    "h|?" => \$shelp,
    "help" => \$help,
    );

pod2usage(1) if $shelp;
pod2usage(exitstatus => 0, verbose => 2) if $help;
pod2usage(1) unless (@ARGV || -p STDIN || -f STDIN);

my $connect = "dbi:mysql(AutoCommit=>1,RaiseError=>1):database=$db";
$connect .= ";host=$host" if $host;
$connect .= ";port=$port" if $port;

my $dbh = DBI->connect($connect,
		       $user,
		       $pass
		      ) or die $DBI::errstr;

my %annot_types = ();
my %domains = (NODOM=>0);
my %domain_clan = (NODOM => {clan_id => 'NODOM', clan_acc=>0, domain_cnt=>0});
my @domain_list = (0);
my $domain_cnt = 0;

my $get_annot_sub = \&get_pfam_annots;

my $get_pfam_acc = $dbh->prepare(<<EOSQL);

SELECT seq_start, seq_end, auto_pfamA, pfamA_acc, pfamA_id, auto_pfamA_reg_full, domain_evalue_score as evalue, length
FROM pfamseq
JOIN pfamA_reg_full_significant using(auto_pfamseq)
JOIN pfamA USING (auto_pfamA)
WHERE in_full = 1
AND  pfamseq_acc=?
ORDER BY seq_start

EOSQL

my $get_pfam_refacc = $dbh->prepare(<<EOSQL);

SELECT seq_start, seq_end, auto_pfamA, pfamA_acc, pfamA_id, auto_pfamA_reg_full, domain_evalue_score as evalue, length
FROM pfamseq
JOIN pfamA_reg_full_significant using(auto_pfamseq)
JOIN pfamA USING (auto_pfamA)
JOIN seqdb_demo2.annot as sa1 on(sa1.acc=pfamseq_acc and sa1.db='sp')
JOIN seqdb_demo2.annot as sa2 using(prot_id)
WHERE in_full = 1
AND  sa2.acc=?
AND  sa2.db='ref'
ORDER BY seq_start

EOSQL

my $get_annots_sql = $get_pfam_acc;

my $get_pfam_id = $dbh->prepare(<<EOSQL);

SELECT seq_start, seq_end, auto_pfamA, pfamA_acc, pfamA_id, auto_pfamA_reg_full, domain_evalue_score as evalue, length
FROM pfamseq
JOIN pfamA_reg_full_significant using(auto_pfamseq)
JOIN pfamA USING (auto_pfamA)
WHERE in_full=1
AND  pfamseq_id=?
ORDER BY seq_start

EOSQL

my $get_pfam_clan = $dbh->prepare(<<EOSQL);

SELECT clan_acc, clan_id
FROM clans
JOIN clan_membership using(auto_clan)
WHERE auto_pfamA=?

EOSQL

my $get_rpd2_clans = $dbh->prepare(<<EOSQL);

SELECT auto_pfamA, clan
FROM ljm_db.RPD2_final_fams
WHERE clan is not NULL

EOSQL

# -- LEFT JOIN clan_membership USING (auto_pfamA)
# -- LEFT JOIN clans using(auto_clan)

my ($tmp, $gi, $sdb, $acc, $id, $use_acc);

# get the query
my ($query, $seq_len) = @ARGV;
$seq_len = 0 unless defined($seq_len);

$query =~ s/^>// if ($query);

my @annots = ();

my %rpd2_clan_fams = ();

if ($rpd2_fams) {
  $get_rpd2_clans->execute();
  my ($auto_pfam, $auto_clan);
  while (($auto_pfam, $auto_clan)=$get_rpd2_clans->fetchrow_array()) {
    $rpd2_clan_fams{$auto_pfam} = $auto_clan;
  }
}

#if it's a file I can open, read and parse it
unless ($query && $query =~ m/[\|:]/) {

  while (my $a_line = <>) {
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
      my ($a_name, $a_num) = ($annot->[-1],$domains{$annot->[-1]});
      if ($acc_comment) {
	$annot->[-1] .= "{$domain_list[$a_num]}";
      }
      $annot->[-1] .= $color_sep_str.$a_num;
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
  $get_annots_sql = $get_pfam_acc;

  if ($annot_line =~ m/^pf26\|/) {
    ($sdb, $gi, $acc, $id) = split(/\|/,$annot_line);
    $dbh->do("use RPD2_pfam");
  }
  elsif ($annot_line =~ m/^gi\|/) {
    ($tmp, $gi, $sdb, $acc, $id) = split(/\|/,$annot_line);
    if ($sdb =~ m/ref/) {
	$get_annots_sql = $get_pfam_refacc;
    }
  }
  elsif ($annot_line =~ m/^sp\|/) {
    ($sdb, $acc, $id) = split(/\|/,$annot_line);
  }
  elsif ($annot_line =~ m/^ref\|/) {
    ($sdb, $acc) = split(/\|/,$annot_line);
    $get_annots_sql = $get_pfam_refacc;
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
    $get_annots_sql = $get_pfam_id;
    $get_annots_sql->execute($id);
  }
  else {
    $acc =~ s/\.\d+$//;
    $get_annots_sql->execute($acc);
  }

  $annot_data{list} = $get_annot_sub->($get_annots_sql, $seq_len);

  return \%annot_data;
}

sub get_pfam_annots {
  my ($get_annots, $seq_length) = @_;

  $seq_length = 0 unless $seq_length;

  my @pf_domains = ();

  # get the list of domains, sorted by start
  while ( my $row_href = $get_annots->fetchrow_hashref()) {
    if ($auto_reg) {
      $row_href->{info} = $row_href->{auto_pfamA_reg_full};
    }
    elsif ($pf_acc) {
      $row_href->{info} = $row_href->{pfamA_acc};
    }
    else {
      $row_href->{info} = $row_href->{pfamA_id};
    }

    if ($row_href && $row_href->{length} > $seq_length && $seq_length == 0) { $seq_length = $row_href->{length};}

    next if ($row_href->{seq_start} >= $seq_length);
    if ($row_href->{seq_end} > $seq_length) {
	$row_href->{seq_end} = $seq_length;
    }

    push @pf_domains, $row_href
  }

  # check for domain overlap, and resolve check for domain overlap
  # (possibly more than 2 domains), choosing the domain with the best
  # evalue

  if($no_over && scalar(@pf_domains) > 1) {

    my @tmp_domains = @pf_domains;
    my @save_domains = ();

    my $prev_dom = shift @tmp_domains;

    while (my $curr_dom = shift @tmp_domains) {

      my @overlap_domains = ($prev_dom);

      my $diff = $prev_dom->{seq_end} - $curr_dom->{seq_start};
      # check for overlap > domain_length/3

      my ($prev_len, $cur_len) = ($prev_dom->{seq_end}-$prev_dom->{seq_start}+1, $curr_dom->{seq_end}-$curr_dom->{seq_start}+1);
      my $inclusion = ((($curr_dom->{seq_start} >= $prev_dom->{seq_start}) && ($curr_dom->{seq_end} <= $prev_dom->{seq_end})) ||
		       (($curr_dom->{seq_start} <= $prev_dom->{seq_start}) && ($curr_dom->{seq_end} >= $prev_dom->{seq_end})));

      my $longer_len = ($prev_len > $cur_len) ? $prev_len : $cur_len;

      while ($inclusion || ($diff > 0 && $diff > $longer_len/3)) {
	push @overlap_domains, $curr_dom;
	$curr_dom = shift @tmp_domains;
	last unless $curr_dom;
	$diff = $prev_dom->{seq_end} - $curr_dom->{seq_start};
	($prev_len, $cur_len) = ($prev_dom->{seq_end}-$prev_dom->{seq_start}+1, $curr_dom->{seq_end}-$curr_dom->{seq_start}+1);
	$longer_len = ($prev_len > $cur_len) ? $prev_len : $cur_len;
	$inclusion = ((($curr_dom->{seq_start} >= $prev_dom->{seq_start}) && ($curr_dom->{seq_end} <= $prev_dom->{seq_end})) ||
		      (($curr_dom->{seq_start} <= $prev_dom->{seq_start}) && ($curr_dom->{seq_end} >= $prev_dom->{seq_end})));
      }

      # check for overlapping domains; >1 because $prev_dom is always there
      if (scalar(@overlap_domains) > 1 ) {
	# if $rpd2_fams, check for a chosen one
	if ($rpd2_fams) {
	  for my $dom (@overlap_domains) {
	    if ($rpd2_clan_fams{$dom->{auto_pfamA}}) {
	      $prev_dom = $dom;
	      last;
	    }
	  }
	}
	else {
	  @overlap_domains = sort { $a->{evalue} <=> $b->{evalue} } @overlap_domains;
	  $prev_dom = $overlap_domains[0];
	}
      }

      # $prev_dom should be the best of the overlaps, and we are no longer overlapping > dom_length/3
      push @save_domains, $prev_dom;
      $prev_dom = $curr_dom;
    }
    if ($prev_dom) {push @save_domains, $prev_dom;}

    @pf_domains = @save_domains;

    # now check for smaller overlaps
    for (my $i=1; $i < scalar(@pf_domains); $i++) {
      if ($pf_domains[$i-1]->{seq_end} >= $pf_domains[$i]->{seq_start}) {
	my $overlap = $pf_domains[$i-1]->{seq_end} - $pf_domains[$i]->{seq_start};
	$pf_domains[$i-1]->{seq_end} -= int($overlap/2);
	$pf_domains[$i]->{seq_start} = $pf_domains[$i-1]->{seq_end}+1;
      }
    }
  }

  if ($neg_doms) {
    my @npf_domains;
    my $prev_dom={seq_end=>0};
    for my $curr_dom ( @pf_domains) {
      if ($curr_dom->{seq_start} - $prev_dom->{seq_end} > $min_nodom) {
	my %new_dom = (seq_start=>$prev_dom->{seq_end}+1, seq_end => $curr_dom->{seq_start}-1, info=>'NODOM');
	push @npf_domains, \%new_dom;
      }
      push @npf_domains, $curr_dom;
      $prev_dom = $curr_dom;
    }
    if ($seq_length - $prev_dom->{seq_end} > $min_nodom) {
      my %new_dom = (seq_start=>$prev_dom->{seq_end}+1, seq_end=>$seq_length, info=>'NODOM');
      if ($new_dom{seq_end} > $new_dom{seq_start}) {push @npf_domains, \%new_dom;}
    }

    # @npf_domains has both old @pf_domains and new neg-domains
    @pf_domains = @npf_domains;
  }

  # now make sure we have useful names: colors

  for my $pf (@pf_domains) {
    $pf->{info} = domain_name($pf->{info}, $pf->{auto_pfamA}, $pf->{pfamA_acc});
  }

  my @feats = ();
  for my $d_ref (@pf_domains) {
    if ($lav) {
      push @feats, [$d_ref->{seq_start}, $d_ref->{seq_end}, $d_ref->{info}];
    }
    else {
      push @feats, [$d_ref->{seq_start}, '-', $d_ref->{seq_end},  $d_ref->{info} ];
#      push @feats, [$d_ref->{seq_end}, ']', '-', ""];
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

  my ($value, $auto_pfamA, $pfamA_acc) = @_;

  # check for clan:
  if ($no_clans) {
    if (! defined($domains{$value})) {
      $domain_clan{$value} = 0;
      $domains{$value} = ++$domain_cnt;
      push @domain_list, $pfamA_acc;
    }
  }
  elsif (!defined($domain_clan{$value})) {
    ## only do this for new domains, old domains have known mappings

    ## ways to highlight the same domain:
    # (1) for clans, substitute clan name for family name
    # (2) for clans, use the same color for the same clan, but don't change the name
    # (3) for clans, combine family name with clan name, but use colors based on clan

    # check to see if it's a clan
    $get_pfam_clan->execute($auto_pfamA);

    my $pfam_clan_href=0;

    if ($pfam_clan_href=$get_pfam_clan->fetchrow_hashref()) {  # is a clan
      my ($clan_id, $clan_acc) = @{$pfam_clan_href}{qw(clan_id clan_acc)};

      # now check to see if we have seen this clan before (if so, do not increment $domain_cnt)
      my $c_value = "C." . $clan_id;
      if ($pf_acc) {$c_value = "C." . $clan_acc;}

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
	push @domain_list, $pfamA_acc;
      }
    }
    else {			# not a clan
      $domain_clan{$value} = 0;
      $domains{$value} = ++$domain_cnt;
      push @domain_list, $pfamA_acc;
    }
  }
  elsif ($domain_clan{$value} && $domain_clan{$value}->{clan_acc}) {
    if ($pf_acc) {$value = "C." . $domain_clan{$value}->{clan_acc};}
    else { $value = "C." . $domain_clan{$value}->{clan_id}; }
  }

  return $value;
}

__END__

=pod

=head1 NAME

ann_feats.pl

=head1 SYNOPSIS

 ann_pfam_e.pl --neg-doms  'sp|P09488|GSTM1_NUMAN' | accession.file

=head1 OPTIONS

 -h	short help
 --help include description
 --no-over  : generate non-overlapping domains (equivalent to ann_pfam.pl)
 --no-clans : do not use clans with multiple families from same clan
 --neg-doms : report domains between annotated domains as NODOM
                 (also --neg, --neg_doms)
 --min_nodom=10  : minimum length between domains for NODOM

 --host, --user, --password, --port --db : info for mysql database

=head1 DESCRIPTION

C<ann_pfam_e.pl> extracts domain information from the pfam msyql
database.  Currently, the program works with database sequence
descriptions in one of two formats:

 Currently, the program works with database
sequence descriptions in several formats:

 >gi|1705556|sp|P54670.1|CAF1_DICDI
 >sp|P09488|GSTM1_HUMAN
 >sp:CALM_HUMAN 

C<ann_pfam_e.pl> uses the C<pfamA_reg_full_significant>, C<pfamseq>,
and C<pfamA> tables of the C<pfam> database to extract domain
information on a protein. 

If the "--no-over" option is set, overlapping domains are selected and
edited to remove overlaps.  For proteins with multiple overlapping
domains (domains overlap by more than 1/3 of the domain length),
C<auto_pfam_e.pl> selects the domain annotation with the best
C<domain_evalue_score>.  When domains overlap by less than 1/3 of the
domain length, they are shortened to remove the overlap.

C<ann_pfam_e.pl> is designed to be used by the B<FASTA> programs with
the C<-V \!ann_pfam_e.pl> or C<-V "\!ann_pfam_e.pl --neg"> option.

=head1 AUTHOR

William R. Pearson, wrp@virginia.edu

=cut
