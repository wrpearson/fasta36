#!/usr/bin/perl -w

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
