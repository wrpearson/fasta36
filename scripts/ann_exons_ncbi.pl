#!/usr/bin/env perl

# ann_exons_ncbi.pl gets an annotation file from fasta36 -V with a line of the form:

## modified 17-Dec-2020 to allow [NXYW]P_ accessions
##

# gi|23065544|ref|NP_000552.2|   or
# NP_000552
#
# and returns the exons present in the protein from NCBI gff3 tables (human, mouse, rat, xtrop)
#
# it must:
# (1) read in the line
# (2) parse it to get the acc
# (3) return the tab delimited exon boundaries
#

use warnings;
use strict;

use DBI;
use Getopt::Long;
use Pod::Usage;

use vars qw($host $db $port $user $pass);

my $hostname = `/bin/hostname`;

($host, $db, $port, $user, $pass)  = ("wrpxdb.its.virginia.edu", "seqdb_demo2", 0, "web_user", "fasta_www");

my ($lav, $shelp, $help) = (0, 0, 0);

my $color_sep_str = " :";
$color_sep_str = '~';

GetOptions(
    "host=s" => \$host,
    "db=s" => \$db,
    "user=s" => \$user,
    "password=s" => \$pass,
    "port=i" => \$port,
    "lav" => \$lav,
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

# get the query -- which could be an acession/length OR a filename
my ($query, $seq_len) = @ARGV;
$seq_len = 0 unless defined($seq_len);

$query =~ s/^>// if ($query);

my @annots = ();

# if it's a file I can open, read and parse it
# check to see if it looks like an accession
## unless ($query && ($query =~ m/[\|:]/ || $query =~ m/^[NXYW]P_/)) {
# it would be better to check to see if a file could be opened
if (! $query || -r $query) {
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
  else {
    $acc = $annot_line;
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

    $row_href->{info} = "exon_".$row_href->{ex_num}.$color_sep_str.$row_href->{ex_num};
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

__END__

=pod

=head1 NAME

ann_exons_ncbi.pl

=head1 SYNOPSIS

 ann_exons_ncbi.pl NP_000552

=head1 OPTIONS

 -h	short help
 --help include description
 --lav  produce lav2plt.pl annotation format, only show domains/repeats
 --host, --user, --password, --port --db -- info for mysql database

=head1 DESCRIPTION

C<ann_exons_ncbi.pl> extracts domain information from a msyql
database.  Currently, the program works with database sequence
descriptions in one of two formats:

  >gi|23065544|ref|NP_000552.2|   or
  >NP_000552

C<ann_exons_ncbi.pl> uses the C<ref_exons> table of the C<seqdb2>
database to extract exon position information on a protein.  The
C<seqdb2/ref_exons> table is constructed from refseq gff files using
the C<ncbi_refseq_ex2prot.pl> script.

C<ann_exons_ncbi.pl> is designed to be used by the B<FASTA> programs with
the C<-V \!ann_exons_ncbi.pl> option.

=head1 AUTHOR

William R. Pearson, wrp@virginia.edu

=cut
