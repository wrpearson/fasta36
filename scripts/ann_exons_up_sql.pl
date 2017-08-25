#!/usr/bin/perl -w

################################################################
# copyright (c) 2014,2015 by William R. Pearson and The Rector &
# Visitors of the University of Virginia
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

# ann_exons_up_sql.pl gets an annotation file from fasta36 -V with a line of the form:

# gi|62822551|sp|P00502|GSTA1_RAT Glutathione S-transfer\n  (at least from pir1.lseg)
#
# it must:
# (1) read in the line
# (2) parse it to get the up_acc
# (3) return the tab delimited features
#

# this version can read feature2 uniprot features (acc/pos/end/label/value), but returns sorted start/end domains
# modified 18-Jan-2016 to produce annotation symbols consistent with ann_exons_up_www2.pl

use strict;

use DBI;
use Getopt::Long;
use Pod::Usage;

use vars qw($host $db $a_table $port $user $pass);

my %domains = ();
my $domain_cnt = 0;

my $hostname = `/bin/hostname`;

unless ($hostname =~ m/ebi/) {
  ($host, $db, $a_table, $port, $user, $pass)  = ("wrpxdb.its.virginia.edu", "uniprot", "annot2", 0, "web_user", "fasta_www");
#  $host = 'xdb';
}
else {
  ($host, $db, $a_table, $port, $user, $pass)  = ("mysql-pearson-prod", "up_db", "annot", 4124, "web_user", "fasta_www");
}

my ($sstr, $lav, $neg_doms, $no_vars, $no_doms, $no_feats, $shelp, $help, $pfam26) = (0,0,0,0,0,0,0,0,0,0);
my ($min_nodom) = (10);

my ($show_color) = (1);
my $color_sep_str = " :";
$color_sep_str = '~';

GetOptions(
    "host=s" => \$host,
    "db=s" => \$db,
    "user=s" => \$user,
    "password=s" => \$pass,
    "port=i" => \$port,
    "lav" => \$lav,
    "no_doms" => \$no_doms,
    "no-doms" => \$no_doms,
    "nodoms" => \$no_doms,
    "no_var" => \$no_vars,
    "no-var" => \$no_vars,
    "novar" => \$no_vars,
    "neg" => \$neg_doms,
    "neg_doms" => \$neg_doms,
    "neg-doms" => \$neg_doms,
    "negdoms" => \$neg_doms,
    "min_nodom=i" => \$min_nodom,
    "min-nodom=i" => \$min_nodom,
    "no_feats" => \$no_feats,
    "no-feats" => \$no_feats,
    "nofeats" => \$no_feats,
    "color!" => \$show_color,
    "sstr" => \$sstr,
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


my $get_annot_sub = \&get_fasta_annots;
if ($lav) {
  $no_feats = 1;
  $get_annot_sub = \&get_lav_annots;
}

my $get_annots_id = $dbh->prepare(qq(select acc, start, end, ix from up_exons join annot2 using(acc) where id=? order by ix));
my $get_annots_acc = $dbh->prepare(qq(select acc, start, end, ix from up_exons where acc=? order by ix));
my $get_annots_refacc = $dbh->prepare(qq(select ref_acc, start, end, ix from up_exons join annot2 using(acc) where ref_acc=? order by ix));

my $get_annots_sql = $get_annots_acc;

my ($tmp, $gi, $sdb, $acc, $id, $use_acc);

# get the query
my ($query, $seq_len) =  @ARGV;
$seq_len = 0 unless defined($seq_len);

$query =~ s/^>// if ($query);

my @annots = ();

#if it's a file I can open, read and parse it
unless ($query && ($query =~ m/[\|:]/ ||
		   $query =~ m/^[NX]P_/ ||
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
    if (!$lav && $show_color && defined($domains{$annot->[-1]})) {
      $annot->[-1] .= $color_sep_str.$domains{$annot->[-1]};
    }
    print join("\t",@$annot),"\n";
  }
}

exit(0);

sub show_annots {
  my ($query_len, $get_annot_sub) = @_;

  my ($annot_line, $seq_len) = split(/\t/,$query_len);

  my %annot_data = (seq_info=>$annot_line);

  if ($annot_line =~ m/^gi\|/) {
    $use_acc = 1;
    ($tmp, $gi, $sdb, $acc, $id) = split(/\|/,$annot_line);
  }
  elsif ($annot_line =~ m/^(SP|TR):(\w+) (\w+)/) {
    ($sdb, $id, $acc) = ($1,$2,$3);
    $use_acc = 1;
    $sdb = lc($sdb)
  }
  elsif ($annot_line =~ m/^(SP|TR):(\w+)/) {
    ($sdb, $id) = ($1,$2);
    $use_acc = 0;
    $sdb = lc($sdb)
  }
  elsif ($annot_line !~ m/\|/) {  # new NCBI swissprot format
    $use_acc =1;
    $sdb = 'sp';
    ($acc) = split(/\s+/,$annot_line);
  }
  else {
    $use_acc = 1;
    ($sdb, $acc, $id) = split(/\|/,$annot_line);
  }

  # remove version number
  unless ($use_acc) {
    $get_annots_sql = $get_annots_id;
    $get_annots_sql->execute($id);
  }
  else {
    unless ($sdb =~ m/ref/) {
      $get_annots_sql = $get_annots_acc;
    } else {
      $get_annots_sql = $get_annots_refacc;
    }
    $acc =~ s/\.\d+$//;
    $get_annots_sql->execute($acc);
  }

  $annot_data{list} = $get_annot_sub->($get_annots_sql, $seq_len);

  return \%annot_data;
}

sub get_fasta_annots {
  my ($get_annots_sql, $seq_len) = @_;

  my ($acc, $start, $end, $ix);
  my @feats = ();

  while (($acc, $start, $end, $ix) = $get_annots_sql->fetchrow_array()) {
      push @feats, [$start, "-", $end, "exon_$ix~$ix"];
  }

  return \@feats;
}

sub get_lav_annots {
  my ($get_annots_sql, $seq_len) = @_;

  my ($pos, $end, $label, $value, $comment);

  my @feats = ();

  my %annot = ();
  while (($acc, $pos, $end, $label, $value) = $get_annots_sql->fetchrow_array()) {
    next unless ($label =~ m/^DOMAIN/ || $label =~ m/^REPEAT/);
    $value =~ s/\s?\{.+\}\.?$//;
    $value = domain_name($label,$value);
    push @feats, [$pos, $end, $value];
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

  my ($label, $value) = @_;

  if ($label =~ /DOMAIN|REPEAT/) {
    $value =~ s/;.*$//;
    $value =~ s/\s+\d+\.?$//;
    $value =~ s/\.\s*$//;
    $value =~ s/\s+\d+\.\s+.*$//;
    $value =~ s/\s+/_/;
    if (!defined($domains{$value})) {
      $domain_cnt++;
      $domains{$value} = $domain_cnt;
    }
    return $value;
  }
  else {
    return $value;
  }
}

__END__

=pod

=head1 NAME

ann_exons_up_sql.pl

=head1 SYNOPSIS

 ann_exons_up_sql.pl --no_doms --no_feats --lav 'sp|P09488|GSTM1_NUMAN' | accession.file

=head1 OPTIONS

 -h	short help
 --help include description
 --no-doms  do not show domain boundaries (domains are always shown with --lav)
 --no-feats do not show features (variants, active sites, phospho-sites)
 --no-var   do not show variant sites (--no_var, --novar)
 --lav  produce lav2plt.pl annotation format, only show domains/repeats
 --neg-doms,  -- report domains between annotated domains as NODOM
                 (also --neg, --neg_doms)
 --min_nodom=10  minimum non-domain length to produce NODOM
 --host, --user, --password, --port --db -- info for mysql database

=head1 DESCRIPTION

C<ann_exons_up_sql.pl> extracts feature, domain, and repeat information from
a msyql database (default name, uniprot) built by parsing the
uniprot_sprot.dat and uniprot_trembl.dat feature tables.  Given a
command line argument that contains a sequence accession (P09488) or
identifier (GSTM1_HUMAN), the program looks up the features available
for that sequence and returns them in a tab-delimited format:

 >sp|P09488|GSTM1_HUMAN
 2	-	88	GST_N-terminal~1
 7	V	F	Mutagen: Reduces catalytic activity 100- fold. {ECO:0000269|PubMed:16548513}.
 34	*	-	MOD_RES: Phosphothreonine. {ECO:0000250|UniProtKB:P10649}.
 90	-	208	GST_C-terminal~2
 108	V	S	Mutagen: Changes the properties of the enzyme toward some substrates. {ECO:0000269|PubMed:16548513, ECO:0000269|PubMed:9930979}.
 108	V	Q	Mutagen: Reduces catalytic activity by half. {ECO:0000269|PubMed:16548513, ECO:0000269|PubMed:9930979}.
 109	V	I	Mutagen: Reduces catalytic activity by half. {ECO:0000269|PubMed:16548513}.
 116	#	-	BINDING: Substrate.
 116	V	A	Mutagen: Reduces catalytic activity 10-fold. {ECO:0000269|PubMed:16548513}.
 116	V	F	Mutagen: Slight increase of catalytic activity. {ECO:0000269|PubMed:16548513}.
 173	V	N	in allele GSTM1B; dbSNP:rs1065411. {ECO:0000269|Ref.3, ECO:0000269|Ref.5}.
 210	*	-	MOD_RES: Phosphoserine. {ECO:0000250|UniProtKB:P04905}.
 210	V	T	in dbSNP:rs449856.

If features are provided, then a legend of feature symbols is provided
as well:

 ==:Active site
 =*:Modified
 =#:Substrate binding
 =^:Site
 =!:Metal binding

If the C<--lav> option is specified, domain and repeat features are
presented in a different format for the C<lav2plt.pl> program:

  >sp|P09488|GSTM1_HUMAN
  2	88	GST N-terminal.
  90	208	GST C-terminal.

C<ann_exons_up_sql.pl> is designed to be used by the B<FASTA> programs
with the C<-V \!ann_exons_up_sql.pl> option, or by the
C<annot_blast_btop.pl> script.  It can also be used with the
lav2plt.pl program with the C<--xA "\!ann_exons_up_sql.pl --lav"> or
C<--yA "\!ann_exons_up_sql.pl --lav"> options.

=head1 AUTHOR

William R. Pearson, wrp@virginia.edu

=cut
