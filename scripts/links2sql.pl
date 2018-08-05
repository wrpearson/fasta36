#!/usr/bin/env perl

################################################################
# copyright (c) 2014 by William R. Pearson and The Rector &
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

use warnings;
use strict;
use DBI;
use Getopt::Long;

use vars qw($host $db $user $password $table $tab_file);

GetOptions(
    'host=s'=>\$host,
    'db=s'=>\$db,
    'user=s'=>\$user,
    'pass=s'=>\$password,
    'table=s'=>\$table,
    'file=s'=>\$tab_file,
    );

$tab_file ||= "link_tmp.tab";

my $dbh = DBI->connect("dbi:mysql:host=$host:$db",
		       $user, $password,
                       { RaiseError => 1, AutoCommit => 1}
                      ) or die $DBI::errstr;

$dbh->do(qq(drop table if exists $table;));
$dbh->do(qq(create table $table (seed_acc varchar(20) not NULL, link_acc varchar(20) not NULL, key seed_acc (seed_acc), key link_acc (link_acc));));

open(FH, ">$tab_file");

while (my $seed_line = <> ) {
    chomp($seed_line);
    my ($seed, $hit_line) = split(/\s+/,$seed_line);
    my @hits = split(/;/,$hit_line);
    for my $hit (@hits) {
	if ($hit ne $seed) {print FH "$seed\t$hit\n";}
    }
}
close(FH);

$dbh->do(qq(load data local infile '$tab_file' into table $table;));

unlink($tab_file);

$dbh->disconnect();
