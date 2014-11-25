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

## usage - expand_uniref50.pl uref50acc.file > up_fasta.file
#

# (1) take a list of uniref50 accessions and uses uniref50link to get
#     the associated uniprot accessions
# (2) take the uniprot accessions and produce a fasta library file
#     from them

use strict;
use DBI;

my ($host, $db, $port, $user, $pass)  = ("xdb", "uniprot", 0, "web_user", "fasta_www");

my $connect = "dbi:mysql(AutoCommit=>1,RaiseError=>1):database=$db";
$connect .= ";host=$host" if $host;
$connect .= ";port=$port" if $port;

my $dbh = DBI->connect($connect,
		       $user,
		       $pass
		      ) or die $DBI::errstr;

my %up_sth = ( 
    ur50_to_upacc => "SELECT uniprot_acc FROM  uniref50link WHERE uniref50_acc=?",
    upacc_to_seq => "SELECT * FROM annot2 join protein USING(acc) WHERE acc=?",
    );

for my $sth (keys(%up_sth)) {
  $up_sth{$sth} = $dbh->prepare($up_sth{$sth});
}

my %acc_uniq = ();

while (my $line = <>) {
  next if ($line =~ m/^UniRef50_UPI/);	# _UPI accessions are not in sp-trembl
  chomp($line);
  my ($up_acc, $e_val) = split(/\t/,$line);
  processLine($up_acc,$up_sth{ur50_to_upacc});
}

for my $up_acc ( keys %acc_uniq ) {

  $up_sth{upacc_to_seq}->execute($up_acc);
  while (my $row_href = $up_sth{upacc_to_seq}->fetchrow_hashref ) {
    print ">sp|". $row_href->{acc} . "|". $row_href->{id} . " (uref50|$acc_uniq{$up_acc}) " .
      $row_href->{descr}. "\n";
    print $row_href->{seq} . "\n";
  }
  $up_sth{upacc_to_seq}->finish();
}

$dbh->disconnect();

sub processLine{
  my ($id,$sth)=@_;

  $id=~ s/UniRef50_//;
  my $result = $sth->execute($id);

  while (my ($acc) = $sth->fetchrow_array()) {
    $acc_uniq{$acc} = $id unless $acc_uniq{$acc};
  }
  $sth->finish();
}
