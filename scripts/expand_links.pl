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

## usage - expand_link.pl seed_acc.file > linked_fasta.file
##
## take a fasta36 -e expand.sh result file of the form:
## gi|12346|ref|NP_98765.1|<tab>1.1e-50
##
## and extract the accession number, looking it up from the an SQL
## table $table. This script uses the database created by link2sql.pl
## Code is included for linking to UniRef and as well as NCBI refseq
## searches.

## Once the linked accession numbers are found, the sequences are
## extracted from the SQL database seqdb_demo2 (see Mackey and Pearson
## (2004) Current Protocols in Bioinformatics (L. Stein, ed) "Using
## SQL databases for sequence similarity searching and analysis".
## Alternatively, one could use blastdbcmd or fastacmd to extract the
## sequences from an NCBI blast-formatted database.
##

use warnings;
use strict;
use DBI;

my ($host, $db, $port, $user, $pass, $table)  = ("xdb", "wrp_link", 0, "web_user", "fasta_www", "micr_samp_link50");

my $dbh = DBI->connect("dbi:mysql:host=$host:$db",
		       $user, $password,
                       { RaiseError => 1, AutoCommit => 1}
                      ) or die $DBI::errstr;

my %sth = ( 
    seed2link => "SELECT link_acc FROM $table WHERE seed_acc=?",
    link2seq => "SELECT * FROM seqdb_demo2.annot JOIN seqdb_demo2.protein USING(prot_id) WHERE acc=? AND pref=1"
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
    print ">". $row_href->{db} . "|". $row_href->{acc} . " (micr_samp|$acc_uniq{$acc}) " .
      $row_href->{descr}. "\n";
    print $row_href->{seq} . "\n";
  }
  $sth{link2seq}->finish();
}

$dbh->disconnect();

sub processLine{
  my ($id,$sth)=@_;
  my ($link_acc);

  if ($id =~ m/^gi\|/) {
      # $id of the form: gi|12346|ref|NP_98765.1|<tab>1.1e-50
      $link_acc = (split(/\|/,$id))[3];
      $link_acc =~ s/\.\d+$//;
  }
  elsif ($id =~ m/^UniRef/) {
      $link_acc = $id;
      $link_acc =~ s/^UniRef\d+_//;
  }
  else {$link_acc = $id;}

  my $result = $sth->execute($link_acc);

  while (my ($acc) = $sth->fetchrow_array()) {
    next if ($acc eq $link_acc);
    $acc_uniq{$acc} = $link_acc unless $acc_uniq{$acc};
  }
  $sth->finish();
}
