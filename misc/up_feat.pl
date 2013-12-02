#!/usr/bin/perl -w

# up_feat.pl gets an annotation file from fasta36 -V with a line of the form:

# gi|62822551|sp|P00502|GSTA1_RAT Glutathione S-transfer\n  (at least from pir1.lseg)
#
# it must:
# (1) read in the line
# (2) parse it to get the up_acc
# (3) return the tab delimited features
#

use strict;

use DBI;
use Getopt::Long;

my $host = "localhost";
my $db = "uniprot";
my $user = "user-name";
my $pass = "password";
my $port = 0;

my $debug = 0;
my $tag = "";
my $verbose = 0;
GetOptions("debug=i" => \$debug,
	   "verbose=i" => \$verbose,
	   "host=s" => \$host,
	   "db=s" => \$db,
	   "user=s" => \$user,
	   "password=s" => \$pass,
	   "port=i" => \&port,
	  );

my $connect = "dbi:mysql(AutoCommit=>1,RaiseError=>1):database=$db";
$connect .= ";host=$host" if $host;
$connect .= ";port=$port" if $port;

my $dbh = DBI->connect($connect,
		       $user,
		       $pass
		      ) or die $DBI::errstr;

my %h_feature = ( ACT_SITE => '@',
		  MOD_RES => '*',
		  BINDING => '^',
		  SITE => '@',
		  METAL => '!',
		  VARIANT => 'V',
		  MUTAGEN => 'V',
    );

my $get_annots = $dbh->prepare('select * from features where acc=? order by pos');

my ($tmp, $gi, $sdb, $acc, $name);

while (my $annot_line = <>) {
  print ">$annot_line";
  if ($annot_line =~ m/^gi\|/) {
     ($tmp, $gi, $sdb, $acc, $name) = split(/\|/,$annot_line);
  }
  elsif ($annot_line =~ m/(\w+):(\w+)/) {
    $sdb = $1;
    $acc = $2;
  }
  else {
    ($sdb, $acc, $name) = split(/\|/,$annot_line);
  }

  # remove version number
  $acc =~ s/\.\d+$//;

  $get_annots->execute($acc);

  my ($pos, $label, $value, $comment);

  while (($acc, $pos, $label, $value) = $get_annots->fetchrow_array()) {
    if ($h_feature{$label}) {
      if ($label =~ m/VARIANT/) {
	my ($aa_res, $comment) = split(/\(/,$value);
	if ($comment) {	
	    $comment =~ s/\)//;
# remove the  /FTId=VAR_014497 information
	    $comment =~ s/\s+\/FTId=.*$//;
	}
	else {$comment = "";}
	next if ($comment =~ /MISSING/);
	my ($vfrom, $vto) = ($aa_res =~ m/(\w)\s*->\s*(\w)/);
	$comment = '' unless $comment;
	if ($vto) {
	  $value = $vto;
	  print join("\t",($pos, $h_feature{$label}, $value, $comment)),"\n";
	}
      } elsif ($label =~ m/MUTAGEN/) {
	my ($aa_res, $comment) = split(/: /,$value);
	next if ($comment =~ /MISSING/);
	my ($vfrom, $vto) = split(/\->/,$aa_res);
	my @vto_list = ();
	if ($vto) {
	  @vto_list = split(/,/,$vto);
	  $value = $vto;
	  for my $val ( @vto_list) {
	    print join("\t",($pos, $h_feature{$label}, $val, "Mutagen: $comment")),"\n";
	  }
	}
      } else {
# uncomment to remove comments on non-variants
#	print join("\t",($pos, $h_feature{$label})),"\n";

# uncomment to show comments on non-variants
	print join("\t",($pos, $h_feature{$label}, "-", "$label: $value")),"\n";
      }
    }
  }
}

exit(0);
