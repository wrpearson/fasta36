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

# ann_pfam.pl gets an annotation file from fasta36 -V with a line of the form:

# gi|62822551|sp|P00502|GSTA1_RAT Glutathione S-transfer\n  (at least from pir1.lseg)
#
# it must:
# (1) read in the line
# (2) parse it to get the up_acc
# (3) return the tab delimited features
#

# this is the first version that works with the new Pfam strategy of
# separating Uniprot reference sequences from the rest of uniprot.  as
# a result, it is possible that 2 SQL queries will be required, one to
# pfamA_reg_full_significant and a second to uniprot_reg_full.

# modified 15-Jan-2017 to reduce the number of calls when the same
# accession is present multiple times.  Accessions are saved in a hash
# than ensures uniqueness.
#
# Uses tmp_annot.temporary table for more rapid joins.  $user must have
# create temporary tables/select permissions for tmp_annot
#

use warnings;
use strict;

use DBI;
use Getopt::Long;
use Pod::Usage;
use File::Temp qw/tempfile/;

use vars qw($host $db $port $user $pass);

my $hostname = `/bin/hostname`;

($host, $db, $port, $user, $pass)  = ("wrpxdb.its.virginia.edu", "pfam30", 0, "web_user", "fasta_www");
#$host = 'xdb';
#$host = 'localhost';
#$db = 'RPD2_pfam28u';

my ($auto_reg,$rpd2_fams, $neg_doms, $vdoms, $lav, $no_doms, $no_clans, $pf_acc, $acc_comment, $bound_comment, $shelp, $help) = 
  (0, 0, 0, 0, 0,0, 0, 0, 0, 0, 0, 0,);
my ($no_over, $split_over, $over_fract) = (0, 0, 3.0);

my ($color_sep_str, $show_color) = (" :",1);
$color_sep_str = '~';

my ($min_nodom, $min_vdom) = (10,10);

GetOptions(
    "host=s" => \$host,
    "db=s" => \$db,
    "user=s" => \$user,
    "password=s" => \$pass,
    "port=i" => \$port,
    "lav" => \$lav,
    "acc_comment" => \$acc_comment,
    "bound_comment" => \$bound_comment,
    "color!" => \$show_color,
    "no-over" => \$no_over,
    "no_over" => \$no_over,
    "split-over" => \$split_over,
    "split_over" => \$split_over,
    "over_fract" => \$over_fract,
    "over-fract" => \$over_fract,
    "no-clans" => \$no_clans,
    "no_clans" => \$no_clans,
    "neg" => \$neg_doms,
    "neg_doms" => \$neg_doms,
    "neg-doms" => \$neg_doms,
    "min_nodom=i" => \$min_nodom,
    "vdoms" => \$vdoms,
    "v_doms" => \$vdoms,
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

my $pfamA_reg_full = 'pfamA_reg_full_significant';
my $uniprot_reg_full = 'uniprot_reg_full';

my @pfam_fields = qw(seq_start seq_end model_start model_end model_length pfamA_acc pfamA_id auto_pfamA_reg_full evalue length);
my @upfam_fields = qw(seq_start seq_end model_start model_end model_length pfamA_acc pfamA_id auto_uniprot_reg_full length);

my $get_pfam_acc = $dbh->prepare(<<EOSQL);
SELECT t_acc, seq_start, seq_end, model_start, model_end, model_length, pfamA_acc, pfamA_id, auto_pfamA_reg_full, domain_evalue_score as evalue, length
FROM pfamseq
JOIN pfamA_reg_full_significant using(pfamseq_acc)
JOIN pfamA USING (pfamA_acc)
JOIN tmp_annot.targets on(pfamseq_acc=t_acc)
WHERE in_full = 1
ORDER BY t_acc,seq_start

EOSQL

my $get_upfam_acc = $dbh->prepare(<<EOSQL);
SELECT t_acc, seq_start, seq_end, model_start, model_end, model_length, pfamA_acc, pfamA_id, auto_uniprot_reg_full as auto_pfamA_reg_full, domain_evalue_score as evalue, length
FROM uniprot
JOIN uniprot_reg_full using(uniprot_acc)
JOIN pfamA USING (pfamA_acc)
JOIN tmp_annot.targets on(uniprot_acc=t_acc)
WHERE in_full = 1
ORDER BY t_acc, seq_start

EOSQL

my $get_pfam_refacc = $dbh->prepare(<<EOSQL);
SELECT t_acc, seq_start, seq_end, model_start, model_end, model_length, pfamA_acc, pfamA_id, auto_pfamA_reg_full, domain_evalue_score as evalue, length
  FROM $pfamA_reg_full
  JOIN pfamseq using(pfamseq_acc)
  JOIN pfamA USING (pfamA_acc)
  JOIN uniprot.refseq2up as rf2up on(rf2up.up_acc=pfamseq_acc)
  JOIN tmp_annot.targets on(rf2up.refseq_acc=t_acc)
 WHERE in_full = 1
 ORDER BY t_acc, seq_start

EOSQL

my $get_upfam_refacc = $dbh->prepare(<<EOSQL);
SELECT t_acc,seq_start, seq_end, model_start, model_end, model_length, pfamA_acc, pfamA_id, auto_uniprot_reg_full as auto_pfamA_reg_full, domain_evalue_score as evalue, length
  FROM uniprot
  JOIN uniprot_reg_full using(uniprot_acc)
  JOIN pfamA USING (pfamA_acc)
  JOIN uniprot.refseq2up as rf2up on(rf2up.up_acc=uniprot_acc)
  JOIN tmp_annot.targets on(rf2up.refseq_acc=t_acc)
WHERE in_full = 1
ORDER BY t_acc, seq_start

EOSQL

my $get_annots_sql = $get_pfam_acc;
my $get_annots_sql_u = $get_upfam_acc;

my $get_pfam_id = $dbh->prepare(<<EOSQL);
SELECT t_acc, seq_start, seq_end, model_start, model_end, model_length, pfamA_acc, pfamA_id, auto_pfamA_reg_full, domain_evalue_score as evalue, length
  FROM pfamseq
  JOIN $pfamA_reg_full using(pfamseq_acc)
  JOIN pfamA USING (pfamA_acc)
  JOIN tmp_annot.targets on(pfamseq_id=t_acc)
WHERE in_full=1
ORDER BY t_acc, seq_start

EOSQL

my $get_upfam_id = $dbh->prepare(<<EOSQL);
SELECT t_acc, seq_start, seq_end, model_start, model_end, model_length, pfamA_acc, pfamA_id, auto_uniprot_reg_full as auto_pfamA_reg_full, domain_evalue_score as evalue, length
  FROM uniprot
  JOIN uniprot_reg_full using(pfamseq_acc)
  JOIN pfamA USING (pfamA_acc)
  JOIN tmp_annot.targets on(uniprot_id=t_acc)
WHERE in_full=1
ORDER BY t_acc, seq_start

EOSQL

my $get_pfam_clan = $dbh->prepare(<<EOSQL);

SELECT clan_acc, clan_id
FROM clan
JOIN clan_membership using(clan_acc)
WHERE pfamA_acc=?

EOSQL

my $get_rpd2_clans = $dbh->prepare(<<EOSQL);

SELECT auto_pfamA, clan
FROM ljm_db.RPD2_final_fams
WHERE clan is not NULL

EOSQL

$dbh->do(<<EOSQL);
create temporary table tmp_annot.targets (t_acc char(10) primary key)
EOSQL

# -- LEFT JOIN clan_membership USING (auto_pfamA)
# -- LEFT JOIN clans using(auto_clan)

my ($tmp, $gi, $sdb, $acc, $id, $use_acc);

my @lav_list = qw(seq_start seq_end info);
my @no_lav_list = qw(seq_start dash seq_end info);
my $out_list_r = \@no_lav_list;
if ($lav) {
  $show_color = 0;
  $out_list_r = \@lav_list;
}

# get the query
my ($query, $seq_len) = @ARGV;
$seq_len = 0 unless defined($seq_len);

$query =~ s/^>// if ($query);

my @annots = ();
my %annot_set = ();

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
    push @annots, show_annots($a_line);
  }
}
else {
  push @annots, show_annots("$query\t$seq_len");
}

# @annots has a list of id's or annotations
# write to temporary file, load data local infile, join to get results

my ($fh, $temp_file) = tempfile(TEMPLATE=>'accannXXXXX');
my @u_annots = keys %annot_set;
print $fh join("\n",@u_annots);
close($fh);

$dbh->do("load data local infile '$temp_file' into table tmp_annot.targets");

unlink($temp_file);

# $get_annots_sql->execute();
# while (my $annot_ar = $get_annots_sql->fetchrow_arrayref()) {
#   my %annot_data = ();
#   @annot_data{@pfam_fields} = @{$annot_ar}[1..10];
#   if (!defined($annot_set{$annot_ar->[0]}->{list})) {
#     $annot_set{$annot_ar->[0]}->{list} = [\%annot_data];
#   }
#   else {
#     push @{$annot_set{$annot_ar->[0]}->{list}}, \%annot_data;
#   }
# }

$get_annots_sql_u->execute();
while (my $annot_hr = $get_annots_sql_u->fetchrow_hashref()) {
  if (!defined($annot_set{$annot_hr->{t_acc}}->{list})) {
    $annot_set{$annot_hr->{t_acc}}->{list} = [$annot_hr];
  }
  else {
    push @{$annot_set{$annot_hr->{t_acc}}->{list}}, $annot_hr;
  }
}

for my $u_acc (@u_annots) {
  map_pfam_annots($annot_set{$u_acc});
}

for my $seq_annot (@annots) {
  next unless $seq_annot;
  my $annot_r = $annot_set{$seq_annot};
  print ">",$annot_r->{seq_info},"\n";
  for my $annot (@{$annot_r->{list}}) {
    $annot->{dash} = '-';
    if (defined($domains{$annot->{info}})) {
      my ($a_name, $a_num) = domain_num($annot->{info},$domains{$annot->{info}});
      $annot->{info} = $a_name;
      my $tmp_a_num = $a_num;
      $tmp_a_num =~ s/v$//;
      if ($acc_comment) {
	$annot->{info} .= "{$domain_list[$tmp_a_num]}";
      }
      if ($bound_comment) {
	$annot->{info} .= $color_sep_str.$annot->{seq_start}.":".$annot->{seq_end};
      }
      elsif ($show_color) {
	  $annot->{info} .= $color_sep_str.$a_num;
      }
    }
    print join("\t",@{$annot}{@{$out_list_r}}),"\n";
  }
}

exit(0);

sub show_annots {
  my ($query_len) = @_;

  my ($annot_line, $seq_len) = split(/\t/,$query_len);

  my $pfamA_acc;

  $use_acc = 1;
  $get_annots_sql = $get_pfam_acc;
  $get_annots_sql_u = $get_upfam_acc;

  if ($annot_line =~ m/^pf\d+\|/) {
    ($sdb, $gi, $pfamA_acc, $acc, $id) = split(/\|/,$annot_line);
#    $dbh->do("use RPD2_pfam");
  }
  elsif ($annot_line =~ m/^gi\|/) {
    ($tmp, $gi, $sdb, $acc, $id) = split(/\|/,$annot_line);
    if ($sdb =~ m/ref/) {
	$get_annots_sql = $get_pfam_refacc;
	$get_annots_sql_u = $get_upfam_refacc;
    }
  }
  elsif ($annot_line =~ m/^(sp|tr)\|/) {
    ($sdb, $acc, $id) = split(/\|/,$annot_line);
  }
  elsif ($annot_line =~ m/^ref\|/) {
    ($sdb, $acc) = split(/\|/,$annot_line);
    $get_annots_sql = $get_pfam_refacc;
    $get_annots_sql_u = $get_upfam_refacc;
  }
  elsif ($annot_line =~ m/^(SP|TR):/i) {
    ($sdb, $id) = split(/:/,$annot_line);
    $use_acc = 0;
  }
  elsif ($annot_line !~ m/\|/ && $annot_line !~ m/:/) {
    $use_acc = 1;
    ($acc) = split(/\s+/,$annot_line);
  }
  # deal with no-database SwissProt/NR
  else {
    ($acc)=($annot_line =~ /^(\S+)/);
  }

  # here we have an $acc or an $id: check to see if we have the data

  my %annot_data = (seq_info=>$annot_line, length=>$seq_len);
  my $annot_key = '';
  unless ($use_acc) {
    next if ($annot_set{$id});
    $annot_set{$id} = \%annot_data;
    $annot_key = $id;

    $get_annots_sql = $get_pfam_id;
  }
  else {
    unless ($acc) {
      warn "missing acc in $annot_line";
      return "";
    }
    else {
      $acc =~ s/\.\d+$//;

      next if ($annot_set{$acc});
      $annot_set{$acc} = \%annot_data;
      $annot_key = $acc;
    }
  }

  return $annot_key;
}

sub map_pfam_annots {
  my ($annot_ref) = @_;

  my $seq_length = $annot_ref->{length};
  my $pf_domains_r = $annot_ref->{list};

  my $row_href=$annot_ref->{list}[0];
  if ($row_href->{length} && $row_href->{length} > $seq_length && $seq_length == 0) {
      $annot_ref->{length} = $seq_length = $row_href->{length};
    }

  # fill in {info} field
  for my $pf_dom (@$pf_domains_r) {
    if ($auto_reg) {
      $pf_dom->{info} = $pf_dom->{auto_pfamA_reg_full};
    } elsif ($pf_acc) {
      $pf_dom->{info} = $pf_dom->{pfamA_acc};
    } else {
      $pf_dom->{info} = $pf_dom->{pfamA_id};
    }
  }

  # before checking for domain overlap, check for "split-domains"
  # (self-unbound) by looking for runs of the same domain that are
  # ordered by model_start

  if (scalar(@{$pf_domains_r}) > 1) {
    my @j_domains;		#joined domains
    my @tmp_domains = @{$pf_domains_r};

    my $prev_dom = shift(@tmp_domains);

    for my $curr_dom (@tmp_domains) {
      # to join domains:
      # (1) the domains must be in order by model_start/end coordinates
      # (3) joining the domains cannot make the total combination too long

      # check for model and sequence consistency
      if (($prev_dom->{pfamA_acc} eq $curr_dom->{pfamA_acc})  # same family
	  && $prev_dom->{model_start} < $curr_dom->{model_start}  # model check
	  && $prev_dom->{model_end} < $curr_dom->{model_end}

	  && ($curr_dom->{model_start} > $prev_dom->{model_end} * 0.80   # limit overlap
	      || $curr_dom->{model_start} <  $prev_dom->{model_end} * 1.25)
	  && ((($curr_dom->{model_end} - $curr_dom->{model_start}+1)/$curr_dom->{model_length} +
	       ($prev_dom->{model_end} - $prev_dom->{model_start}+1)/$prev_dom->{model_length}) < 1.33)
	 ) {			# join them by updating $prev_dom
	$prev_dom->{seq_end} = $curr_dom->{seq_end};
	$prev_dom->{model_end} = $curr_dom->{model_end};
	$prev_dom->{auto_pfamA_reg_full} = $prev_dom->{auto_pfamA_reg_full} . ";". $curr_dom->{auto_pfamA_reg_full};
	$prev_dom->{evalue} = ($prev_dom->{evalue} < $curr_dom->{evalue} ? $prev_dom->{evalue} : $curr_dom->{evalue});
      } else {
	push @j_domains, $prev_dom;
	$prev_dom = $curr_dom;
      }
    }
    push @j_domains, $prev_dom;
    @{$pf_domains_r} = @j_domains;


    if ($no_over) {	# for either $no_over or $split_over, check for overlapping domains and edit/split them

      my @tmp_domains = @{$pf_domains_r};	# allow shifts from copy of @pf_domains
      my @save_domains = ();		# where the new domains go

      my $prev_dom = shift @tmp_domains;

      while (my $curr_dom = shift @tmp_domains) {

	my @overlap_domains = ($prev_dom);

	my $diff = $prev_dom->{seq_end} - $curr_dom->{seq_start};

	my ($prev_len, $cur_len) = ($prev_dom->{seq_end}-$prev_dom->{seq_start}+1,
				    $curr_dom->{seq_end}-$curr_dom->{seq_start}+1);

	my $inclusion = ((($curr_dom->{seq_start} >= $prev_dom->{seq_start})    # start is right && end is left 
			  && ($curr_dom->{seq_end} <= $prev_dom->{seq_end})) || # -- curr inside prev
			 (($curr_dom->{seq_start} <= $prev_dom->{seq_start})    # start is left && end is right
			  && ($curr_dom->{seq_end} >= $prev_dom->{seq_end})));  # -- prev is inside curr

	my $longer_len = ($prev_len > $cur_len) ? $prev_len : $cur_len;

	# check for overlap > domain_length/$over_fract
	while ($inclusion || ($diff > 0 && $diff > $longer_len/$over_fract)) {
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

      if ($prev_dom) {
	push @save_domains, $prev_dom;
      }

      @{$pf_domains_r} = @save_domains;

      # now check for smaller overlaps
      for (my $i=1; $i < scalar(@{$pf_domains_r}); $i++) {
	if ($pf_domains_r->[$i-1]->{seq_end} >= $pf_domains_r->[$i]->{seq_start}) {
	  my $overlap = $pf_domains_r->[$i-1]->{seq_end} - $pf_domains_r->[$i]->{seq_start};
	  $pf_domains_r->[$i-1]->{seq_end} -= int($overlap/2);
	  $pf_domains_r->[$i]->{seq_start} = $pf_domains_r->[$i-1]->{seq_end}+1;
	}
      }
    }
    elsif ($split_over) {   # here, everything that overlaps by > $min_vdom should be split into a separate domain
      my @save_domains = ();		# where the new domains go

      # check to see if one domain is included (or overlapping) more
      # than xx% of the other.  If so, pick the longer one

      my ($prev_dom, $curr_dom) = ($pf_domains_r->[0],0) ;
      for (my $i=1; $i < scalar(@{$pf_domains_r}); $i++) {
	$curr_dom = $pf_domains_r->[$i];

	my ($prev_len, $cur_len) = ($prev_dom->{seq_end}-$prev_dom->{seq_start}+1, $curr_dom->{seq_end}-$curr_dom->{seq_start}+1);
	my $longer_len = ($prev_len > $cur_len) ? $prev_len : $cur_len;

	if (($curr_dom->{seq_start} >= $prev_dom->{seq_start}) && ($curr_dom->{seq_end} <= $prev_dom->{seq_end})
	   && $cur_len / $prev_len > 0.80) {
	  # $prev_dom stays the same, $curr_dom deleted
	  next;
	}
	elsif (($curr_dom->{seq_start} <= $prev_dom->{seq_start}) && ($curr_dom->{seq_end} >= $prev_dom->{seq_end})
	      && $prev_len / $cur_len > 0.80) {
	  $prev_dom = $curr_dom; # this should delete $prev_dom
	  next;
	}

	if ($prev_dom->{seq_end} >= $curr_dom->{seq_start} + $min_vdom) {
	  my ($l_seq_end, $r_seq_start) = ($curr_dom->{seq_start}-1, $prev_dom->{seq_end}+1);

	  $prev_dom->{seq_end} = $l_seq_end;
	  push @save_domains, $prev_dom;
	  my $new_dom = {seq_start => $l_seq_end+1, seq_end=>$r_seq_start-1, 
			 model_length => -1,
			 pfamA_acc=>$prev_dom->{pfamA_acc}."/".$curr_dom->{pfamA_acc},
			 pfamA_id=>$prev_dom->{pfamA_id}."/".$curr_dom->{pfamA_id},
			 };

	  if ($pf_acc) {
	    $new_dom->{info} = $new_dom->{pfamA_acc};
	  }
	  else {
	    $new_dom->{info} = $new_dom->{pfamA_id};
	  }

	  push @save_domains, $new_dom;
	  $curr_dom->{seq_start} = $r_seq_start;
	  $prev_dom = $curr_dom;
	}
	else {
	  push @save_domains, $prev_dom;
	  $prev_dom = $curr_dom;
	}
      }
      push @save_domains, $prev_dom;
      @{$pf_domains_r} = @save_domains;
    }
  }

  # $vdoms -- virtual Pfam domains -- the equivalent of $neg_doms,
  # but covering parts of a Pfam model that are not annotated.  split
  # domains have been joined, so simply check beginning and end of
  # each domain (but must also check for bounded-ness)
  # only add when 10% or more is missing and missing length > $min_nodom

  if ($vdoms && scalar(@{$pf_domains_r})) {
    my @vpf_domains;

    my $curr_dom = $pf_domains_r->[0];
    my $length = $curr_dom->{length};

    my $prev_dom={seq_end=>0, pfamA_acc=>''};
    my $prev_dom_end = 0;
    my $next_dom_start = $length+1;

    for (my $dom_ix=0; $dom_ix < scalar(@{$pf_domains_r}); $dom_ix++ ) {
      $curr_dom = $pf_domains_r->[$dom_ix];

      my $pfamA =  $curr_dom->{pfamA_acc};

      # first, look left, is there a domain there (if there is,
      # it should be updated right

      # my $min_vdom = $curr_dom->{model_length} / 10;

      if ($curr_dom->{model_length} < $min_vdom) {
	push @vpf_domains, $curr_dom;
	next;
      }
      if ($prev_dom->{pfamA_acc}) { # look for previous domain
	$prev_dom_end = $prev_dom->{seq_end};
      }

      # there is a domain to the left, how much room is available?
      my $left_dom_len = min($curr_dom->{seq_start}-$prev_dom_end-1, $curr_dom->{model_start}-1);
      if ( $left_dom_len > $min_vdom) {
	# there is room for a virtual domain
	my %new_dom = (seq_start=> $curr_dom->{seq_start}-$left_dom_len,
	               seq_end => $curr_dom->{seq_start}-1,
		       info=>'@'.$curr_dom->{info},
		       model_length=>$curr_dom->{model_length},
		       model_end => $curr_dom->{model_start}-1,
		       model_start => $left_dom_len,
		       pfamA_acc=>$pfamA,
		      );
	push @vpf_domains, \%new_dom;
      }

      # save the current domain
      push @vpf_domains, $curr_dom;
      $prev_dom = $curr_dom;

      if ($dom_ix < scalar(@$pf_domains_r)-1) { # there is a domain to the right
	# first, give all the extra space to the first domain (no splitting)
	$next_dom_start = $pf_domains_r->[$dom_ix+1]->{seq_start};
      }
      else {
	$next_dom_start = $length;
      }

      # is there room for a virtual domain right
	  
      my $right_dom_len = min($next_dom_start-$curr_dom->{seq_end}-1, # space available 
			      $curr_dom->{model_length}-$curr_dom->{model_end} # space needed
			     );
      if ( $right_dom_len > $min_vdom) {
	my %new_dom = (seq_start=> $curr_dom->{seq_end}+1,
		       seq_end=> $curr_dom->{seq_end}+$right_dom_len,
		       info=>'@'.$curr_dom->{info},
		       model_length => $curr_dom->{model_length},
		       pfamA_acc=> $pfamA,
		      );
	push @vpf_domains, \%new_dom;
	$prev_dom = \%new_dom;
      }
    }				# all done, check for last one

    # $curr_dom=$pf_domains_r->[-1];
    # # my $min_vdom = $curr_dom->{model_length}/10;

    # my $right_dom_len = min($length - $curr_dom->{seq_end}+1,  # space available 
    # 			    $curr_dom->{model_length}-$curr_dom->{model_end} # space needed
    # 			   );
    # if ($right_dom_len > $min_vdom) {
    #   my %new_dom = (seq_start=> $curr_dom->{seq_end}+1,
    # 		     seq_end => $curr_dom->{seq_end}+$right_dom_len,
    # 		     info=>'@'.$curr_dom->{pfamA_acc},
    # 		     model_len=> $curr_dom->{model_len},
    # 		     pfamA_acc => $curr_dom->{pfamA_acc},
    # 		     model_start => $curr_dom->{model_end}+1,
    # 		     model_end => $curr_dom->{model_len},
    # 		     );

    #   push @vpf_domains, \%new_dom;
    # }

    # @vpf_domains has both old @{$pf_domains_r} and new neg-domains
    @{$pf_domains_r} = @vpf_domains;
  }

  if ($neg_doms) {
    my @npf_domains;
    my $prev_dom={seq_end=>0};
    for my $curr_dom ( @{$pf_domains_r}) {
      if ($curr_dom->{seq_start} - $prev_dom->{seq_end} > $min_nodom) {
	my %new_dom = (seq_start=>$prev_dom->{seq_end}+1, seq_end => $curr_dom->{seq_start}-1, info=>'NODOM');
	push @npf_domains, \%new_dom;
      }
      push @npf_domains, $curr_dom;
      $prev_dom = $curr_dom;
    }
    if ($seq_length - $prev_dom->{seq_end} > $min_nodom) {
      my %new_dom = (seq_start=>$prev_dom->{seq_end}+1, seq_end=>$seq_length, info=>'NODOM');
      if ($new_dom{seq_end} > $new_dom{seq_start}) {
	push @npf_domains, \%new_dom;
      }
    }

    # @npf_domains has both old @pf_domains and new neg-domains
    @{$pf_domains_r} = @npf_domains;
  }

  # now make sure we have useful names: colors

  for my $pf (@{$pf_domains_r}) {
    $pf->{info} = domain_name($pf->{info}, $pf->{pfamA_acc});
  }
}

sub min {
  my ($arg1, $arg2) = @_;

  return ($arg1 <= $arg2 ? $arg1 : $arg2);
}

sub max {
  my ($arg1, $arg2) = @_;

  return ($arg1 >= $arg2 ? $arg1 : $arg2);
}

# domain name takes a uniprot domain label, removes comments ( ;
# truncated) and numbers and returns a canonical form. Thus:
# Cortactin 6.
# Cortactin 7; truncated.
# becomes "Cortactin"
#

sub domain_name {

  my ($value, $pfamA_acc) = @_;
  my $is_virtual = 0;

  if ($value =~ m/^@/) {
    $is_virtual = 1;
    $value =~ s/^@//;
  }

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
    $get_pfam_clan->execute($pfamA_acc);

    my $pfam_clan_href=0;

    if ($pfam_clan_href=$get_pfam_clan->fetchrow_hashref()) {  # is a clan
      my ($clan_id, $clan_acc) = @{$pfam_clan_href}{qw(clan_id clan_acc)};

      # now check to see if we have seen this clan before (if so, do not increment $domain_cnt)
      my $c_value = "C." . $clan_id;
      if ($pf_acc) {$c_value = $clan_acc;}

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
    if ($pf_acc) {$value = $domain_clan{$value}->{clan_acc};}
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


__END__

=pod

=head1 NAME

ann_pfam30_tmptbl.pl

=head1 SYNOPSIS

 ann_pfam30_tmptbl.pl --neg-doms  --vdoms 'sp|P09488|GSTM1_NUMAN' | accession.file

=head1 OPTIONS

 -h	short help
 --help include description
 --no-over  : generate non-overlapping domains (equivalent to ann_pfam.pl)
 --split-over  : overlaps of two domains generate a new hybrid domain
 --no-clans : do not use clans with multiple families from same clan
 --neg-doms : report domains between annotated domains as NODOM
                 (also --neg, --neg_doms)
 --vdoms : produce "virtual domains" using model_start, 
           model_end for partial pfam domains
 --min_nodom=10  : minimum length between domains for NODOM

 --host, --user, --password, --port --db : info for mysql database

=head1 DESCRIPTION

C<ann_pfam30_tmptbl.pl> extracts domain information from the pfam msyql
database. Currently, the program works with database
sequence descriptions in several formats:

 >gi|1705556|sp|P54670.1|CAF1_DICDI
 >sp|P09488|GSTM1_HUMAN
 >sp:CALM_HUMAN 

C<ann_pfam30_tmptbl.pl> uses the C<pfamA_reg_full_significant>, C<pfamseq>,
and C<pfamA> tables of the C<pfam> database to extract domain
information on a protein. 

If the C<--no-over> option is set, overlapping domains are selected and
edited to remove overlaps.  For proteins with multiple overlapping
domains (domains overlap by more than 1/3 of the domain length),
C<auto_pfam28.pl> selects the domain annotation with the best
C<domain_evalue_score>.  When domains overlap by less than 1/3 of the
domain length, they are shortened to remove the overlap.

If the C<--split-over> option is set, if two domains overlap, the
overlapping region is split out of the domains and labeled as a new,
virtual-lie, domain.  If one domain is internal to another and spans
80% of the domain, the shorter domain is removed.

C<ann_pfam30_tmptbl.pl> is designed to be used by the B<FASTA> programs with
the C<-V \!ann_pfam30_tmptbl.pl> or C<-V "\!ann_pfam30_tmptbl.pl --neg"> option.

C<ann_pfam30_tmptbl.pl> requires an additional database, C<tmp_annot>,
with C<create temporary tables>, C<insert>, and C<select> privileges
for the default user.

=head1 AUTHOR

William R. Pearson, wrp@virginia.edu

=cut
