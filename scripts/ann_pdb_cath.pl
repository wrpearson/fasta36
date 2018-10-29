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

# ann_feats.pl gets an annotation file from fasta36 -V with a line of the form:

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

use DBI;
use Getopt::Long;
use Pod::Usage;

use vars qw($host $db $port $user $pass);

my $hostname = `/bin/hostname`;

($host, $db, $port, $user, $pass)  = ("wrpxdb.its.virginia.edu", "uniprot", 0, "web_user", "fasta_www");

my ($neg_doms, $lav, $shelp, $help, $class) = (0, 0, 0, 0, 0);
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
    "neg" => \$neg_doms,
    "neg_doms" => \$neg_doms,
    "neg-doms" => \$neg_doms,
    "min_nodom=i" => \$min_nodom,
    "class" => \$class,
    "h|?" => \$shelp,
    "help" => \$help,
    );

pod2usage(1) if $shelp;
pod2usage(exitstatus => 0, verbose => 2) if $help;
pod2usage(1) unless ( @ARGV || -p STDIN || -f STDIN);

my $connect = "dbi:mysql(AutoCommit=>1,RaiseError=>1):database=$db";
$connect .= ";host=$host" if $host;
$connect .= ";port=$port" if $port;

my $dbh = DBI->connect($connect,
		       $user,
		       $pass
		      ) or die $DBI::errstr;

my %domains = (NODOM=>0);
my $domain_cnt = 0;

my $get_offsets_pdb =  $dbh->prepare(<<EOSQL);
SELECT res_beg, pdb_beg, pdb_end, sp_beg, sp_end
FROM pdb_chain_up
WHERE pdb_acc=?
ORDER BY res_beg
EOSQL

my $get_cathdoms_pdb = $dbh->prepare(<<EOSQL);
SELECT s_start, s_stop, p_start, p_stop, cath_class, s_descr as info
FROM cath_doms
JOIN cath_names using(cath_class)
WHERE pdb_acc=?
ORDER BY s_start
EOSQL

my ($tmp, $gi, $sdb, $acc, $id, $use_acc);

# get the query
my ($query, $seq_len) = @ARGV;
$seq_len = 1 unless $seq_len;

$query =~ s/^>// if $query;

my @annots = ();

#if it's a file I can open, read and parse it
unless ($query && $query =~ m/[\|:]/) {

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
  $get_offsets_pdb->execute($pdb_acc);
  my ($res_beg, $pdb_beg, $pdb_end, $sp_beg, $sp_end) = $get_offsets_pdb->fetchrow_array();

  $res_beg = 1 unless defined($res_beg);
  $pdb_beg = 1 unless defined($pdb_beg);
  $sp_beg = 1 unless defined($sp_beg);

  if (defined($sp_end) && $sp_end > $seq_len) {$seq_len = $sp_end;}
  if (defined($pdb_end) && $pdb_end > $seq_len) {$seq_len = $pdb_end;}

  # unless ($seq_len > 1) {
  #   if (defined($sp_end)) {
  #     $seq_len = $sp_end;
  #   }
  #   elsif (defined($pdb_end)) {
  #     $seq_len = $pdb_end;
  #   }
  # }

  $get_cathdoms_pdb->execute($pdb_acc);
  $annot_data{list} = get_cath_annots($lav, $get_cathdoms_pdb, $pdb_beg, $seq_len, $off_flag);

  return \%annot_data;
}

sub get_cath_annots {
  my ($lav, $get_annots, $sp_offset, $seq_length, $is_offset) = @_;

  my @cath_domains = ();

  # get the list of domains, sorted by start
  while ( my $row_href = $get_annots->fetchrow_hashref()) {

    # put in logic to subtract sp_offset when necessary
    if ($is_offset && $row_href->{p_start}) {
      $row_href->{seq_start} = $row_href->{p_start} - $sp_offset;
      $row_href->{seq_end} = $row_href->{p_stop} - $sp_offset;
    }
    else {
      $row_href->{seq_start} = $row_href->{s_start} - $sp_offset;
      $row_href->{seq_end} = $row_href->{s_stop} - $sp_offset;
    }

    if ($seq_length <= 1) {
      $seq_length = $row_href->{seq_end};
    }
    else {
      $row_href->{seq_end} = $seq_length if ($row_href->{seq_end} > $seq_length);
    }
    
    $row_href->{info} =~ s/\s+/_/g;

    push @cath_domains, $row_href
  }

  return unless (scalar(@cath_domains));

  # do a consistency check
  for (my $i=1; $i < scalar(@cath_domains); $i++) {
    if ($cath_domains[$i]->{seq_start} <= $cath_domains[$i-1]->{seq_end}) {
      my $delta = $cath_domains[$i]->{seq_start} - $cath_domains[$i-1]->{seq_end};
      $cath_domains[$i-1]->{seq_end} -= $delta/2;
      $cath_domains[$i]->{seq_start} = $cath_domains[$i-1]->{seq_end}+1;
    }
  }

  if ($neg_doms) {
    my @ncath_domains;
    my $prev_dom={seq_end=>0};
    for my $cur_dom ( @cath_domains) {
      if ($cur_dom->{seq_start} - $prev_dom->{seq_end} > $min_nodom) {
	my %new_dom = (seq_start=>$prev_dom->{seq_end}+1, seq_end => $cur_dom->{seq_start}-1, info=>'NODOM');
	push @ncath_domains, \%new_dom;
      }
      push @ncath_domains, $cur_dom;
      $prev_dom = $cur_dom;
    }
    my %new_dom = (seq_start=>$prev_dom->{seq_end}+1, seq_end=>$seq_length, info=>'NODOM');
    if ($new_dom{seq_end} > $new_dom{seq_start}) {push @ncath_domains, \%new_dom;}

    @cath_domains = @ncath_domains;
  }

  for my $cath (@cath_domains) {
    if ($class && $cath->{cath_class}) {$cath->{info} = $cath->{cath_class};}
    $cath->{info} = domain_name($cath->{info});
  }

  my @feats = ();

  if ($lav) {
    for my $d_ref (@cath_domains) {
      push @feats, [$d_ref->{seq_start}, $d_ref->{seq_end}, $d_ref->{info} ];
    }
  }
  else {
    for my $d_ref (@cath_domains) {
      push @feats, [$d_ref->{seq_start}, '-',  $d_ref->{seq_end}, $d_ref->{info} ];
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

ann_feats.pl

=head1 SYNOPSIS

 ann_pdb_cath.pl --neg 'sp|P09488|GSTM1_NUMAN' | accession.file

=head1 OPTIONS

 -h	short help
 --help include description

 --lav  produce lav2plt.pl annotation format, only show domains/repeats
 --neg-doms,  -- report domains between annotated domains as NODOM
                 (also --neg, --neg_doms)
 --min_nodom=10  -- minimum length between domains for NODOM

 --host, --user, --password, --port --db -- info for mysql database

=head1 DESCRIPTION

C<ann_pfam26.pl> extracts domain information from the pfam26 msyql
database.  Currently, the program works with database sequence
descriptions in one of two formats:

 >pf26|649|O94823|AT10B_HUMAN -- RPD2_seqs

(pf26 databases have auto_pfamseq in the second field) and

 >gi|1705556|sp|P54670.1|CAF1_DICDI

C<ann_pfam26.pl> uses the C<pfamA_reg_full_significant>, C<pfamseq>,
and C<pfamA> tables of the C<pfam26> database to extract domain
information on a protein.  For proteins that have multiple domains
associated with the same overlapping region (domains overlap by more
than 1/3 of the domain length), C<auto_pfam26.pl> selects the domain
annotation with the best C<domain_evalue_score>.  When domains overlap
by less than 1/3 of the domain length, they are shortened to remove
the overlap.

C<ann_pfam26.pl> is designed to be used by the B<FASTA> programs with
the C<-V \!ann_pfam26.pl> or C<-V "\!ann_pfam26.pl --neg"> option.

=head1 AUTHOR

William R. Pearson, wrp@virginia.edu

=cut
