#!/usr/bin/perl -w

# ann_exons_ncbi.pl gets an annotation file from fasta36 -V with a line of the form:

# gi|23065544|ref|NP_000552.2| 
#
# and returns the exons present in the protein from NCBI gff3 tables (human and mouse only)
#
# it must:
# (1) read in the line
# (2) parse it to get the acc
# (3) return the tab delimited exon boundaries


use strict;

use DBI;
use Getopt::Long;
use Pod::Usage;

use vars qw($host $db $port $user $pass);

my $hostname = `/bin/hostname`;

($host, $db, $port, $user, $pass)  = ("wrpxdb.its.virginia.edu", "seqdb_demo2", 0, "web_user", "fasta_www");

my ($auto_reg,$rpd2_fams, $neg_doms, $lav, $no_doms, $pf_acc, $shelp, $help) = (0, 0, 0, 0,0, 0,0,0);
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
my $domain_cnt = 0;

my $get_annot_sub = \&get_refseq_exons;

my $get_exons_acc = $dbh->prepare(<<EOSQL);

SELECT ex_num, ex_p_start as seq_start, ex_p_end as seq_end
FROM ref_exons
WHERE  acc=?
ORDER BY seq_start

EOSQL

my $get_annots_sql = $get_exons_acc;

my ($tmp, $gi, $sdb, $acc, $id, $use_acc);

# get the query
my ($query, $seq_len) = @ARGV;
$seq_len = 0 unless defined($seq_len);

$query =~ s/^>// if ($query);

my @annots = ();

#if it's a file I can open, read and parse it
unless ($query && $query =~ m/[\|:]/ ) {

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
      $annot->[-1] .= $color_sep_str.$domains{$annot->[-1]};
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
  $get_annots_sql = $get_exons_acc;

  if ($annot_line =~ m/^gi\|/) {
    ($tmp, $gi, $sdb, $acc) = split(/\|/,$annot_line);
  }
  elsif ($annot_line =~ m/^ref\|/) {
    ($sdb, $acc) = split(/\|/,$annot_line);
  }

  $acc =~ s/\.\d+$//;
  $get_annots_sql->execute($acc);

  $annot_data{list} = $get_annot_sub->($get_annots_sql, $seq_len);

  return \%annot_data;
}

sub get_refseq_exons {
  my ($get_annots, $seq_length) = @_;

  my @exons = ();

  # get the list of domains, sorted by start
  while ( my $row_href = $get_annots->fetchrow_hashref()) {

    $row_href->{info} = "exon_".$row_href->{ex_num};
    push @exons, $row_href
  }

  # check for domain overlap, and resolve check for domain overlap
  # (possibly more than 2 domains), choosing the domain with the best
  # evalue

  my @feats = ();

  for my $d_ref (@exons) {
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

 ann_pfam.pl --neg-doms  'sp|P09488|GSTM1_NUMAN' | accession.file

=head1 OPTIONS

 -h	short help
 --help include description
 --neg-doms,  -- report domains between annotated domains as NODOM
                 (also --neg, --neg_doms)
 --min_nodom=10  -- minimum length between domains for NODOM

 --host, --user, --password, --port --db -- info for mysql database

=head1 DESCRIPTION

C<ann_pfam.pl> extracts domain information from a msyql
database.  Currently, the program works with database sequence
descriptions in one of two formats:

 >pf26|649|O94823|AT10B_HUMAN -- RPD2_seqs

(pf26 databases have auto_pfamseq in the second field) and

 >gi|1705556|sp|P54670.1|CAF1_DICDI

C<ann_pfam.pl> uses the C<pfamA_reg_full_significant>, C<pfamseq>,
and C<pfamA> tables of the C<pfam> database to extract domain
information on a protein.  For proteins that have multiple domains
associated with the same overlapping region (domains overlap by more
than 1/3 of the domain length), C<auto_pfam.pl> selects the domain
annotation with the best C<domain_evalue_score>.  When domains overlap
by less than 1/3 of the domain length, they are shortened to remove
the overlap.

C<ann_pfam.pl> is designed to be used by the B<FASTA> programs with
the C<-V \!ann_pfam.pl> or C<-V "\!ann_pfam.pl --neg"> option.

=head1 AUTHOR

William R. Pearson, wrp@virginia.edu

=cut
