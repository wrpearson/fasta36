#!/usr/bin/env perl

################################################################
# copyright (c) 2010, 2014 by William R. Pearson and The Rector &
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

## usage - expand_up_isoforms.pl [--prim_acc] up_hits.file > up_isoforms.file
##
## take a fasta36 -e expand.sh result file of the form:
## sp|P09488_GSTM1_HUMAN|<tab>1.1e-50
##
## and extract the accession number, looking it up from the an SQL
## table $table -- in this case "annot2_iso" to provide Uniprot
## isoforms based on a uniprot accession.
##
## if --prim_acc, then the primary accession (used to find the isoforms) is added to the isoform seq_id, e.g.
## sp|P04988|GSTM1_HUMAN has isoforms:   with --prim_acc, the identifiers become
## >iso|E7EWW9|E7EWW9_HUMAN    >iso|E7EWW9|E7EWW9_HUMAN_P09488
## >iso|H3BRM6|H3BRM6_HUMAN    >iso|H3BRM6|H3BRM6_HUMAN_P09488
## >iso|H3BQT3|H3BQT3_HUMAN    >iso|H3BQT3|H3BQT3_HUMAN_P09488

use warnings;
use strict;
use Getopt::Long;
use Pod::Usage;
use DBI;

my ($host, $db, $port, $user, $pass)  = ("xdb", "uniprot", 0, "web_user", "fasta_www");
$host = 'wrpxdb.its.virginia.edu';
my ($a_table, $i_table) = ("annot2", "annot2_iso");
my ($help, $shelp) = (0,0);
my ($e_thresh, $prim_acc) = (1e-6, 0);

GetOptions(
    "h" => \$shelp,
    "help" => \$help,
    "host=s" => \$host,
    "prim_acc!" => \$prim_acc,
    "db=s" => \$db,
    "expect|evalue|e_thresh=f" => \$e_thresh,
    "user=s" => \$user,
    "password=s" => \$pass,
    "port=i" => \$port,
    "i_table" => \$i_table,
    "a_table" => \$a_table,
    );

pod2usage(1) if $shelp;
pod2usage(exitstatus => 0, verbose => 2) if $help;
pod2usage(1) unless (@ARGV || -p STDIN || -f STDIN);

my $dbh = DBI->connect("dbi:mysql:host=$host:$db",
		       $user, $pass,
                       { RaiseError => 1, AutoCommit => 1}
                      ) or die $DBI::errstr;

my %sth = ( 
    seed2link_acc => "SELECT acc FROM $i_table WHERE prim_acc=?",
    seed2link_id => "SELECT iso_a.acc FROM $i_table as iso_a JOIN $a_table  as an2 on(iso_a.prim_acc=an2.acc) where an2.id=?",
    link2seq => "SELECT db, acc, prim_acc, id, descr, seq FROM annot2_iso JOIN protein_iso USING(acc) WHERE acc=?"
    );

for my $sth (keys(%sth)) {
  $sth{$sth} = $dbh->prepare($sth{$sth});
}

my %acc_uniq = ();

# get the query
my ($query, $eval_arg) =  @ARGV;
$eval_arg = 1e-10 unless $eval_arg;
$query =~ s/^>// if ($query);
my @link_lines = ();

#if it's a file I can open, read and parse it
unless ($query && ($query =~ m/[\|:]/ ||
		   $query =~ m/^[NX]P_/ ||
		   $query =~ m/^[OPQ][0-9][A-Z0-9]{3}[0-9]|[A-NR-Z][0-9]([A-Z][A-Z0-9]{2}[0-9]){1,2}\s/)) {

  while (my $a_line = <>) {
    $a_line =~ s/^>//;
    chomp $a_line;
    push @link_lines, $a_line;
  }
}
else {
  push @link_lines, "$query\t$eval_arg";
}

for my $line ( @link_lines ) {
  my ($hit, $e_val) = split(/\t/,$line);

  if ($e_val <= $e_thresh) {
      process_line($hit,$sth{seed2link_acc},$sth{seed2link_id});
  }
}

for my $acc ( keys %acc_uniq ) {

  $sth{link2seq}->execute($acc);
  while (my $row_href = $sth{link2seq}->fetchrow_hashref ) {
    my $id_str = $row_href->{id};
    if ($prim_acc) {
      $id_str .= "_".$row_href->{prim_acc};
    }

    printf(">%s|%s|%s %s\n","iso",$acc,$id_str,$row_href->{descr});
    my $iso_seq = $row_href->{seq};
    $iso_seq =~ s/(.{60})/$1\n/g;

    print "$iso_seq\n";
  }
  $sth{link2seq}->finish();
}

$dbh->disconnect();

sub process_line{
  my ($seqid,$sth_acc, $sth_id)=@_;

  my $sth = $sth_acc;

  my ($db, $link_acc, $link_id) = ("","","");

  if ($seqid =~ m/\|/) {
      ($db, $link_acc, $link_id) = split('\|',$seqid);
      $link_acc =~ s/\.\d+$//;

      $sth_acc->execute($link_acc);
  }
  elsif ($seqid =~ m/:/) {
      ($db, $link_id) = split(':',$seqid);
      $sth_id->execute($link_id);
      $sth = $sth_id;
  }
  else {
      $link_acc = $seqid;
      $link_acc =~ s/\.\d+$//;
      $sth_acc->execute($link_acc);
  }

  while (my ($acc) = $sth->fetchrow_array()) {
    next if ($acc eq $link_acc);
    $acc_uniq{$acc} = $link_acc unless $acc_uniq{$acc};
  }
  $sth->finish();
}

__END__

=pod

=head1 NAME

 expand_up_isoforms.pl expand_file.tab

=head1 SYNOPSIS

 expand_up_isoforms.pl expand_file.tab

=head1 OPTIONS

 -h	short help
 --help include description
 --evalue E()-value threshold for expansion
 --prim_acc : show primary accession as part of sequence identifier
    >iso|E7EWW9|E7EWW9_HUMAN becomes >iso|E7EWW9|E7EWW9_HUMAN_P09488

 --host, --user, --password, --port --db : info for mysql database
 --a_table, --i_table  -- SQL table names with reference and isoform acc/id/prim_acc mappings.

=head1 DESCRIPTION

C<expand_up_isoforms.pl> uses protein isoform tables in an SQL database to identify and extract
isoforms of proteins in a reference protein sequence database.

C<expand_up_isoforms.pl> takes a file with sequece identifiers and E()-values of the form:

 sp|P09488|GSTM1_HUMAN <tab> 1e-40
 sp:CALM_HUMAN <tab> 1e-40

Lines with E()-values less than --evalue (1E-6 by default) are used to
identify protein isoforms, which are included in the set of sequences to be aligned.

C<expand_up_isoforms.pl> is designed to be used by the B<FASTA> programs with
the C<-e expand_up_isoforms.pl> option.

=head1 AUTHOR

William R. Pearson, wrp@virginia.edu

=cut
