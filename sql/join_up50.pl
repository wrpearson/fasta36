#!/usr/bin/perl -w

################################################################
# copyright (c) 2004, 2014 by William R. Pearson and The Rector &
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

# usage - join_up50.pl link_acc_XYZ12 > "link_acc_XYZ12.sql 16" file
#
# this program is designed to read a link_acc file from fasta36 -e "expand.sh" and
# (1) extract the uniref accessions
# (2) build a mySQL fasta_tmp table using those accessions
# (3) write out a link_mysql.sql file that can be used to produce the
#     sequences using format 16
#

use strict;
use DBI;

my $in_file = $ARGV[0];
open (FIN, $in_file) || die "cannot open $in_file\n";

my $db = 'uniprot';
my $db_tmp = 'fasta_tmp';

my $host='host';
my $user = 'user';
my $password = 'password';

my $dbh = DBI->connect("dbi:mysql:host=xdb2:$db",
		       $user, $password,
                       { RaiseError => 1, AutoCommit => 1}
    ) or die $DBI::errstr;

my $dbh_tmp = DBI->connect("dbi:mysql:host=xdb2:$db_tmp",
			   $user, $password,
			   { RaiseError => 1, AutoCommit => 1}
    ) or die $DBI::errstr;

my %up_sth = ( 
    ur50_to_upacc => "SELECT uniprot_acc FROM  uniref50link WHERE uniref50_acc=?",
    upacc_to_seq => "SELECT * FROM  trFull WHERE acc=?",
    );


for my $sth (keys(%up_sth)) {
  $up_sth{$sth} = $dbh->prepare($up_sth{$sth});
}

my %q_acc_uniq = ();

while (my $line = <FIN>) {
  next if ($line =~ m/^UniRef50_UPI/);
  chomp($line);
  my ($descr, $score) = split(/\t/,$line);

  my ($acc) = ($descr =~ m/UniRef\d+_(\w+)/i);
  $q_acc_uniq{$acc} = 1;
}
close FIN;

$dbh->disconnect();

# now we have a hash of unique accessions
# make a table and put them in

$dbh_tmp->do(qq{DROP TABLE IF EXISTS $in_file;}) or die $DBI::errstr;
$dbh_tmp->do(qq{ CREATE TABLE $in_file ( acc CHAR(10) PRIMARY KEY );})
    or die $DBI::errstr;

my @acc_list = sort keys(%q_acc_uniq);

my $sql_insert = "INSERT INTO $in_file (acc) VALUES (\"" . join(q/"),("/,@acc_list) . "\");\n" ;
$dbh_tmp->do($sql_insert);
$dbh_tmp->disconnect();

# now the $in_file table should be full; write out the SQL join to produce the sequences we need.

print "$host $db $user $password;\n";
print qq{SELECT up.acc, protein.seq
 FROM $db_tmp.$in_file AS fa_acc JOIN uniref50link on(fa_acc.acc=uniref50_acc)
 JOIN annot AS up ON(uniprot_acc = up.acc AND fa_acc.acc != up.acc) JOIN protein USING(prot_id);
SELECT acc, concat('up|',acc,'|',name,' ',descr) FROM annot WHERE acc='#';
SELECT acc,protein.seq FROM protein INNER JOIN annot USING(prot_id)
 WHERE annot.acc='#';
DROP TABLE $db_tmp.$in_file;
};
