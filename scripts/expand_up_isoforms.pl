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
use DBI;

my ($host, $db, $port, $user, $pass, $table)  = ("xdb", "uniprot", 0, "web_user", "fasta_www","annot2_iso");

GetOptions(
    "host=s" => \$host,
    "db=s" => \$db,
    "user=s" => \$user,
    "password=s" => \$pass,
    "port=i" => \$port,
    );

my $dbh = DBI->connect("dbi:mysql:host=$host:$db",
		       $user, $pass,
                       { RaiseError => 1, AutoCommit => 1}
                      ) or die $DBI::errstr;

my %sth = ( 
    seed2link => "SELECT acc FROM $table WHERE prim_acc=?",
    link2seq => "SELECT db, acc, id, descr, seq FROM annot2_iso JOIN protein_iso USING(acc) WHERE acc=?"
    );

for my $sth (keys(%sth)) {
  $sth{$sth} = $dbh->prepare($sth{$sth});
}

my %acc_uniq = ();

while (my $line = <>) {
  chomp($line);
  my ($hit, $e_val) = split(/\t/,$line);
  processLine($hit,$sth{seed2link});
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

sub processLine{
  my ($seqid,$sth)=@_;

  my ($db, $link_acc, $id) = split('\|',$seqid);

  $link_acc = $seqid unless ($link_acc);

  my $result = $sth->execute($link_acc);

  while (my ($acc) = $sth->fetchrow_array()) {
    next if ($acc eq $link_acc);
    $acc_uniq{$acc} = $link_acc unless $acc_uniq{$acc};
  }
  $sth->finish();
}
