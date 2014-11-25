#!/usr/bin/perl -w

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
## sp|P09488|<tab>1.1e-50
##
## and extract the accession number, looking it up from the an SQL
## table $table. This script uses the database created by link2sql.pl
## Code is included for linking to UniRef and as well as NCBI refseq
## searches.

## Once the linked accession numbers are found, the sequences are
## extracted from the SQL database uniprot (see Mackey and Pearson
## (2004) Current Protocols in Bioinformatics (L. Stein, ed) "Using
## SQL databases for sequence similarity searching and analysis".
## Alternatively, one could use blastdbcmd or fastacmd to extract the
## sequences from an NCBI blast-formatted database.
##

use strict;
use DBI;

my ($host, $port, $db, $table, $user, $pass);

my $hostname = `/bin/hostname`;

unless ($hostname =~ m/ebi/) {
  ($host, $db, $port, $user, $pass)  = ("xdb", "uniprot", 0, "web_user", "fasta_www");
}
else {
  ($host, $db, $port, $user, $pass)  = ("mysql-pearson", "up_db", 4124, "web_user", "fasta_www");
}

my $connect = "dbi:mysql(AutoCommit=>1,RaiseError=>1):database=$db";
$connect .= ";host=$host" if $host;
$connect .= ";port=$port" if $port;

my $dbh = DBI->connect($connect,
		       $user,
		       $pass
		      ) or die $DBI::errstr;

my %sth = ( 
    up2ensg_id => "SELECT acc, ensg FROM ensg JOIN annot2 USING(acc) WHERE id=?",
    up2ensg_acc => "SELECT * FROM ensg WHERE acc=?",
    ensg2seq => "SELECT * FROM ensg JOIN annot2 USING(acc) JOIN protein USING(acc) WHERE ensg=?",
    );

for my $sth (keys(%sth)) {
  $sth{$sth} = $dbh->prepare($sth{$sth});
}

my %acc_uniq = ();
my %ensg_uniq = ();

while (my $line = <>) {
  chomp($line);
  my ($hit, $e_val) = split(/\t/,$line);
  processLine($hit,$sth{up2ensg});
}

for my $ensg_acc ( keys %ensg_uniq ) {

  $sth{ensg2seq}->execute($ensg_acc);
  while (my $row_href = $sth{ensg2seq}->fetchrow_hashref ) {
    next if ($acc_uniq{$row_href->{acc}});
#    print ">". $row_href->{db} . "|". $row_href->{acc} . " (".
#	$ensg_uniq{$acc}->{acc}."|".$ensg_uniq{$acc}->{id}.":$acc) " .
#      $row_href->{descr}. "\n";
#    print ">" . uc($row_href->{db}) . ":$ensg_uniq{$acc}->{acc} $row_href->{acc} $acc $row_href->{descr}\n";
    print ">",join('|', (lc($row_href->{db}),$row_href->{acc},$row_href->{id}))," ($ensg_uniq{$ensg_acc}->{id}|$ensg_acc) $row_href->{descr}\n";
    print $row_href->{seq} . "\n";
  }
  $sth{ensg2seq}->finish();
}

$dbh->disconnect();

sub processLine{
  my ($id)=@_;
  my ($dummy, $link_acc, $link_id);

  my $use_acc = 1;
  my $get_sth = $sth{up2ensg_acc};

  if ($id =~ m/^gi\|/) {
      # $id of the form: gi|12346|ref|NP_98765.1|<tab>1.1e-50
      ($link_acc, $link_id) = (split(/\|/,$id))[3,4];
      $link_acc =~ s/\.\d+$//;
  }
  elsif ($id =~ m/(\w+):(\w+)/) {
    $link_id = $2;
    $link_acc = '';
    $use_acc = 0;
  }
  elsif ($id =~ m/sp\|([\w\-\.]+)/) {
      ($dummy, $link_acc, $link_id) = split(/\|/,$id);
  }
  elsif ($id =~ m/tr\|([\w\-\.]+)/) {
      ($dummy, $link_acc, $link_id) = split(/\|/,$id);
  }
# form: SP:GSTM1_MOUSE P10649
  elsif ($id =~ m/SP\:(\w+)/) {
      ($link_id) = ($1);
      $use_acc = 0;
  } 
  elsif ($id =~ m/TR\:(\w+)/) {
      ($link_id) = ($1);
      $use_acc = 0;
  }
  else {$link_acc = $id;}

  if ($use_acc) {
    return if ($acc_uniq{$link_acc});
    $get_sth->execute($link_acc);
  }
  else {
    $get_sth = $sth{up2ensg_id};
    $get_sth->execute($link_id);
  }

  while (my ($acc, $ensg) = $get_sth->fetchrow_array()) {
    $acc_uniq{$acc} = $acc unless $acc_uniq{$acc};
    $ensg_uniq{$ensg} = {id=>$acc} unless $ensg_uniq{$ensg};
  }
  $get_sth->finish();
}
