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

## usage - expand_up_isoforms.pl up_hits.file > up_isoforms.file
##
## take a fasta36 -e expand.sh result file of the form:
## sp|P09488_GSTM1_HUMAN|<tab>1.1e-50
##
## and extract the accession number, looking it up from the an SQL
## table $table -- in this case "annot2_iso" to provide Uniprot
## isoforms based on a uniprot accession.

use warnings;
use strict;
use Getopt::Long;
use Pod::Usage;
use DBI;

my ($host, $db, $port, $user, $pass)  = ("xdb", "uniprot", 0, "web_user", "fasta_www");
my ($a_table, $i_table) = ("annot2", "annot2_iso");
my ($help, $shelp) = (0,0);
my ($e_thresh) = 1e-6;


GetOptions(
    "h" => \$shelp,
    "help" => \$help,
    "host=s" => \$host,
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
    link2seq => "SELECT db, acc, id, descr, seq FROM annot2_iso JOIN protein_iso USING(acc) WHERE acc=?"
    );

for my $sth (keys(%sth)) {
  $sth{$sth} = $dbh->prepare($sth{$sth});
}

my %acc_uniq = ();

while (my $line = <>) {
  chomp($line);
  my ($hit, $e_val) = split(/\t/,$line);
  if ($e_val <= $e_thresh) {
      process_line($hit,$sth{seed2link_acc},$sth{seed2link_id});
  }
}

for my $acc ( keys %acc_uniq ) {

  $sth{link2seq}->execute($acc);
  while (my $row_href = $sth{link2seq}->fetchrow_hashref ) {
    printf(">%s|%s|%s %s\n","iso",$acc,$row_href->{id},$row_href->{descr});
    print $row_href->{seq} . "\n";
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
      $sth_acc->execute($link_acc);
  }
  elsif ($seqid =~ m/:/) {
      ($db, $link_id) = split(':',$seqid);
      $sth_id->execute($link_id);
      $sth = $sth_id;
  }
  else {
      $link_acc = $seqid;
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
