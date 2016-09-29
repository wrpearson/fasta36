#!/usr/bin/perl -w

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

# ann_feats2ipr.pl gets an annotation file from fasta36 -V with a line of the form:

# gi|62822551|sp|P00502|GSTA1_RAT Glutathione S-transfer\n  (at least from pir1.lseg)
#
# it must:
# (1) read in the line
# (2) parse it to get the up_acc
# (3) return the tab delimited features
#

# this version takes "features:"
#   ACT_SITE, MOD_RES, SITE, METAL, VARIANT, MUTAGEN
# from Uniprot and combines them with domain annotations from my merge of the Interpro database.
#

# ann_feats2ipr.pl is largely identical to ann_feats2l.pl, except that
#   it uses Interpro for domain/repeat information.

use strict;

use DBI;
use Getopt::Long;
use Pod::Usage;

use vars qw($host $db $dom_db $a_table $port $user $pass);

my %domains = ();
my $domain_cnt = 0;

my $hostname = `/bin/hostname`;

unless ($hostname =~ m/ebi/) {
  ($host, $db, $a_table, $port, $user, $pass)  = ("wrpxdb.its.virginia.edu", "uniprot", "annot2", 0, "web_user", "fasta_www");
#  $host = 'localhost';
} else {
  ($host, $db, $a_table, $port, $user, $pass)  = ("mysql-pearson", "up_db", "annot", 4124, "web_user", "fasta_www");
}

my ($lav, $neg_doms, $no_doms, $no_feats, $no_label, $use_ipr, $acc_comment, $shelp, $help, $no_mod, $dom_db, $db_ref_acc) = 
    (0,0,0,0,0,0,0,0,0,0,0,0);

my $color_sep_str = " :";
$color_sep_str = '~';

GetOptions(
	   "host=s" => \$host,
	   "db=s" => \$db,
	   "user=s" => \$user,
	   "password=s" => \$pass,
	   "port=i" => \$port,
	   "lav" => \$lav,
	   "no_mod" => \$no_mod,
	   "no-mod" => \$no_mod,
	   "no-doms" => \$no_doms,
	   "nodoms" => \$no_doms,
           "dom_db=s" => \$dom_db,
	   "dom_acc" => \$db_ref_acc,
	   "dom-acc" => \$db_ref_acc,
	   "neg" => \$neg_doms,
	   "neg_doms" => \$neg_doms,
	   "neg-doms" => \$neg_doms,
	   "negdoms" => \$neg_doms,
	   "no_feats" => \$no_feats,
	   "no-feats" => \$no_feats,
	   "nofeats" => \$no_feats,
	   "no_label" => \$no_label,
	   "no-label" => \$no_label,
	   "nolabel" => \$no_label,
	   "ipr" => \$use_ipr,
	   "acc_comment" => \$acc_comment,
	   "h|?" => \$shelp,
	   "help" => \$help,
	  );

pod2usage(1) if $shelp;
pod2usage(exitstatus => 0, verbose => 2) if $help;
pod2usage(1) unless (-p STDIN || -f STDIN || @ARGV);

my $connect = "dbi:mysql(AutoCommit=>1,RaiseError=>1):database=$db";
$connect .= ";host=$host" if $host;
$connect .= ";port=$port" if $port;

my $dbh = DBI->connect($connect,
		       $user,
		       $pass
		      ) or die $DBI::errstr;

my @feat_keys = qw(ACT_SITE MOD_RES BINDING SITE METAL VARIANT MUTAGEN);
my @feat_vals = ( '=','*','#','^','!','V','V');
my @feat_text = ( "active site", "phosphorylation", "binding site", "site", "metal binding");

my @dom_vals = ( [ '[', ']'],[ '[', ']']);

my %annot_types = ();
@annot_types{@feat_keys} = @feat_vals;

my $get_annot_sub = \&get_fasta_annots;
if ($lav) {
  $no_feats = 1;
  $get_annot_sub = \&get_lav_annots;
}

if ($neg_doms) {
  $domains{'NODOM'}=0;
}

if ($no_mod) {
  @feat_keys = qw(ACT_SITE BINDING SITE METAL);
  @feat_text = ( "active site", "binding site", "site", "metal binding");
  @feat_vals = ( '=','#','^','!');
  delete($annot_types{'MOD_RES'});
  delete($annot_types{'MUTAGEN'});
  delete($annot_types{'VARIANT'});
}

my $get_ft2_sites_id = $dbh->prepare(qq(select acc, pos, end, label, value, len from features2 join $a_table using(acc) where id=? and label in ('ACT_SITE','MOD_RES','BINDING','SITE','METAL','VARIANT','MUTAGEN') order by pos));

my $get_ft2_sites_acc = $dbh->prepare(qq(select acc, pos, end, label, value, len from features2 join $a_table using(acc) where acc=? and label in ('ACT_SITE','MOD_RES','BINDING','SITE','METAL','VARIANT','MUTAGEN') order by pos));

my $get_ft2_sites_refacc= $dbh->prepare(qq(select ref_acc, pos, end, label, value, len from features2 join $a_table using(acc) where ref_acc=? and label in ('ACT_SITE','MOD_RES','BINDING','SITE','METAL','VARIANT','MUTAGEN') order by pos));

my $get_ipr_doms_id = $dbh->prepare(qq(select acc, start, stop, ipr_acc, db_ref, s_descr, len from prot2ipr_s join $a_table using(acc) join ipr_annot using(ipr_acc) where id=? order by start));

my $get_ipr_domdb_id = $dbh->prepare(qq(select acc, start, stop, ipr_acc, db_ref, s_descr, len from prot2ipr join $a_table using(acc) join ipr_annot using(ipr_acc) where dom_db='$dom_db' AND id=? order by start));

my $get_ipr_doms_acc = $dbh->prepare(qq(select acc, start, stop, ipr_acc, db_ref, s_descr, len from prot2ipr_s join $a_table using(acc) join ipr_annot using(ipr_acc) where acc=? order by start));

my $get_ipr_doms_refacc = $dbh->prepare(qq(select ref_acc, start, stop, ipr_acc, db_ref, s_descr, len from prot2ipr_s join $a_table using(acc) join ipr_annot using(ipr_acc) where ref_acc=? order by start));

my $get_ipr_domdb_acc = $dbh->prepare(qq(select acc, start, stop, ipr_acc, db_ref, s_descr, len from prot2ipr join $a_table using(acc) join ipr_annot using(ipr_acc) where dom_db='$dom_db' AND acc=? order by start));

my $get_sites_sql = $get_ft2_sites_id;
my $get_doms_sql = $get_ipr_doms_id;

my ($tmp, $gi, $sdb, $acc, $id, $use_acc);

unless ($no_feats || $no_label) {
  for my $i ( 0 .. $#feat_text ) {
    print "=",$feat_vals[$i],":",$feat_text[$i],"\n";
  }
  # print "=*:phosphorylation\n";
  # print "==".":active site\n";
  # print "=@".":site\n";
  # print "=^:binding\n";
  # print "=!:metal binding\n";
}

# get the query
my ($query, $seq_len) = @ARGV;

$seq_len = 0 unless defined($seq_len);

$query =~ s/^>// if $query;

my @annots = ();

#if it's a file I can open, read and parse it
unless ($query && $query =~ m/[\|:]/ ) {

  while (my $a_line = <>) {
    $a_line =~ s/^>//;
    chomp $a_line;
    my $annots_ref = show_annots($a_line, $get_annot_sub);
    push @annots, $annots_ref if ($annots_ref);
  }
} else {
  my $annots_ref = show_annots("$query $seq_len", $get_annot_sub);
  push @annots, $annots_ref if ($annots_ref);
}

for my $seq_annot (@annots) {
  print ">",$seq_annot->{seq_info},"\n";
  for my $annot (@{$seq_annot->{list}}) {
    if (!$lav && defined($domains{$annot->[4]})) {
      $annot->[-2] .= $color_sep_str.$domains{$annot->[4]};
    }
    if ($lav) {
      print join("\t",@$annot[0 .. 2]),"\n";
    }
    else {
      print join("\t",@$annot[0 .. 3]),"\n";
    }
  }
}

exit(0);

sub show_annots {
  my ($query_length, $get_annot_sub) = @_;

  my ($annot_line, $seq_length) = split(/\t/,$query_length);

  my %annot_data = (seq_info=>$annot_line);

  if ($annot_line =~ m/^gi\|/ && $annot_line =~ m/\|[sp|ref]\|/) {
    $use_acc = 1;
    ($tmp, $gi, $sdb, $acc, $id) = split(/\|/,$annot_line);
    if ($sdb !~ m/sp/ && $annot_line =~ m/\|sp\|(\w+)/) {
      ($acc) = ($annot_line =~ m/\|sp\|(\w+)/);
    }
  }
  elsif ($annot_line =~ m/^[sp|tr|ref]\|/ ) {
    $use_acc = 1;
    ($sdb, $acc) = split(/\|/,$annot_line);
  }
  elsif ($annot_line =~ m/^gi\|/) {
    $use_acc = 1;
    ($tmp, $gi, $sdb, $acc, $id) = split(/\|/,$annot_line);
  }
  elsif ($annot_line =~ m/SP:(\w+)/) {
    $use_acc = 0;
    $sdb = 'sp';
    $id = $1;
  } elsif ($annot_line =~ m/TR:(\w+)/) {
    $use_acc = 0;
    $sdb = 'tr';
    $id = $1;
  }
  elsif ($annot_line =~ m/\|/) {  # new NCBI swissprot format
    $use_acc = 1;
    ($sdb, $acc, $id) = split(/\|/,$annot_line);
  } else {
    $use_acc =1;
    $sdb = 'sp';
    ($acc) = split(/\s+/,$annot_line);
  }

  # remove version number
  unless ($use_acc) {
    unless ($no_feats) {
      $get_sites_sql = $get_ft2_sites_id;
      $get_sites_sql->execute($id);
    }
    unless ($no_doms) {
      if ($dom_db) {
	$get_doms_sql = $get_ipr_domdb_id;
      }
      else {
	$get_doms_sql = $get_ipr_doms_id;
      }

      $get_doms_sql->execute($id);
    }
  } else {
    unless ($acc) {
      print STDERR "ann_feats2ipr.pl no acc: $annot_line\n";
      return 0;
    }
    $acc =~ s/\.\d+$//;
    if ($sdb eq 'ref') {
      unless ($no_feats) {
	$get_sites_sql = $get_ft2_sites_refacc;
	$get_sites_sql->execute($acc);
      }
      unless ($no_doms) {
	$get_doms_sql = $get_ipr_doms_refacc;
	$get_doms_sql->execute($acc);
      }
    }
    else {
      unless ($no_feats) {
	$get_sites_sql = $get_ft2_sites_acc;
	$get_sites_sql->execute($acc);
      }
      unless ($no_doms) {
	if ($dom_db) {
	  $get_doms_sql = $get_ipr_domdb_acc;
	}
	else {
	  $get_doms_sql = $get_ipr_doms_acc;
	}
	$get_doms_sql->execute($acc);
      }
    }
  }

  $annot_data{list} = $get_annot_sub->($seq_length, $get_sites_sql, $get_doms_sql);

  return \%annot_data;
}

sub get_fasta_annots {
  my ($seq_len, $get_sites_sql, $get_doms_sql) = @_;

  my ($acc, $pos, $end, $label, $value, $comment, $len);

  $seq_len = 0;

  my @feats2 = (); # features with start/stop, for checking overlap, adding negative
  my @sites = ();  # sites with one position

  # get sites
  unless ($no_feats) {
    while (($acc, $pos, $end, $label, $value, $len) = $get_sites_sql->fetchrow_array()) {
      $seq_len = $len if ($len > $seq_len);
      next unless $annot_types{$label};
      if ($label =~ m/VARIANT/) {
	my ($aa_res, $comment) = split(/\(/,$value);
	if ($comment) {	
	  $comment =~ s/\)//;
	  # remove the  /FTId=VAR_014497 information
	  $comment =~ s/\s+\/FTId=.*$//;
	} else {
	  $comment = "";
	}
	next if ($comment =~ /MISSING/);
	my ($vfrom, $vto) = ($aa_res =~ m/(\w)\s*->\s*(\w)/);
	if ($vto) {
	  $comment = '' unless $comment;
	  $value = $vto;
	  push @sites, [$pos, $annot_types{$label}, $value, $comment, ""];
	}
      } elsif ($label =~ m/MUTAGEN/) {
	my ($aa_res, $comment) = split(/: /,$value);
	next if ($comment =~ /MISSING/);
	my ($vfrom, $vto) = split(/\->/,$aa_res);
	if ($vto) {
	  my @vto_list = split(/,/,$vto);
	  $value = $vto;
	  for my $val ( @vto_list) {
	    push @sites, [$pos, $annot_types{$label}, $val, "Mutagen: $comment", ""];
	  }
	}
      } else {
	push @sites, [$pos, $annot_types{$label}, "-", "$label: $value", ""];
      }
    }
  }

  unless ($no_doms) {
    my ($ipr_acc, $db_ref, $s_descr) = ("","","");
    while (($acc, $pos, $end, $ipr_acc, $db_ref, $s_descr, $len) = $get_doms_sql->fetchrow_array()) {
      $db_ref =~ s/G3DSA://;
      $seq_len = $len unless ($seq_len > $len);

      $value = domain_name($ipr_acc, $s_descr);
      if ($acc_comment) {
	$value .= "{$ipr_acc}";
      }
      if ($db_ref_acc) {
	$value = $db_ref;
      }
      elsif ($use_ipr) {
	$value = $ipr_acc;
      }

      push @feats2, [$pos, "-", $end, $value, $ipr_acc];
    }
  }

  # ensure that domains do not overlap
  for (my $i=1; $i < scalar(@feats2); $i++) {
    my $diff = $feats2[$i-1]->[2] - $feats2[$i]->[0];
    if ($diff >= 0) {
      $feats2[$i-1]->[2] = $feats2[$i]->[0]+ int($diff/2);
      $feats2[$i]->[0] = $feats2[$i-1]->[2] + 1;
    }
  }

  my @n_feats2 = ();

  if ($neg_doms) {
    my $last_end = 0;
    for my $feat ( @feats2 ) {
      if ($feat->[0] - $last_end > 10) {
	push @n_feats2, [$last_end+1, "-", $feat->[0]-1, "NODOM", "NODOM"];
      }
      $last_end = $feat->[2];
    }
    if ($seq_len - $last_end > 10) {
      push @n_feats2, [$last_end+1, "-", $seq_len, "NODOM", "NODOM"];
    }
  }

  my @feats = ();
  for my $feat (@feats2, @n_feats2) {
    push @feats, [$feat->[0], '[', '-', $feat->[-2], $feat->[-1] ];
    push @feats, [$feat->[2], ']', '-', "", ""];
  }

  @feats = sort { $a->[0] <=> $b->[0] } (@sites, @feats);

  return \@feats;
}

sub get_lav_annots {
  my ($seq_len, $get_sites_sql, $get_doms_sql) = @_;

  my @feats = ();

  my %annot = ();
  while (my ($acc, $pos, $end, $ipr_acc, $db_ref, $s_descr, $len) = $get_doms_sql->fetchrow_array()) {
    #    $value = domain_name($label,$value);
    my $value = domain_name($ipr_acc,$s_descr);
    push @feats, [$pos, $end, $value];
  }

  return \@feats;
}

# domain name takes a uniprot domain label, removes comments ( ;
# truncated) and numbers and returns a canonical form. Thus:
# Cortactin 6.
# Cortactin 7; truncated.
# becomes "Cortactin"
#

sub domain_name {

  my ($ipr_acc, $s_descr) = @_;

  $s_descr =~ s/[\-_]domain//;
  $s_descr =~ s/[\-_]homology//;

  $s_descr =~ s/^(.{20})/$1/;

  if (!defined($domains{$ipr_acc})) {
      $domain_cnt++;
      $domains{$ipr_acc} = $domain_cnt;
  }
  return $s_descr;
}


__END__

=pod

=head1 NAME

ann_feats2.pl

=head1 SYNOPSIS

 ann_feats2.pl --no_doms --no_feats --lav 'sp|P09488|GSTM1_NUMAN' | accession.file

=head1 OPTIONS

 -h	short help
 --help include description

 --acc_comment provide the InterPro accession in {IPR00123} brackets for links
 --dom_db=G3DSA use a single domain database (e.g. PF, G3DSA, PS5) from InterPro
 --dom_acc provide the domain accession, not the description, as the domain label
 --neg, --neg_doms, --neg-doms label non-domain regions > 10 residues as "NODOM"
 --ipr proide InterPro accession as label
 --no-doms  do not show domain boundaries (domains are always shown with --lav)
 --no-feats do not show feature (variants, active sites, phospho-sites)
 --no-label do show feature key (==*phosphorylation, etc)

 --lav  produce lav2plt.pl annotation format, only show domains/repeats

 --host, --user, --password, --port --db -- info for mysql database

=head1 DESCRIPTION

C<ann_feats2ipr.pl> extracts feature, domain, and repeat information from
two msyql databases (default names: uniprot/ipr2) built by parsing the
uniprot_sprot.dat and uniprot_trembl.dat feature tables.  Given a
command line argument that contains a sequence accession (P09488) or
identifier (GSTM1_HUMAN), the program looks up the features available
for that sequence and returns them in a tab-delimited format:

 >sp|P09488
 2	-	88	DOMAIN: GST N-terminal.
 7	V	F	Mutagen: Reduces catalytic activity 100- fold.
 23	*	-	MOD_RES: Phosphotyrosine (By similarity).
 33	*	-	MOD_RES: Phosphotyrosine (By similarity).
 34	*	-	MOD_RES: Phosphothreonine (By similarity).
 90	-	208	DOMAIN: GST C-terminal.
 108	V	S	Mutagen: Changes the properties of the enzyme toward some substrates.
 108	V	Q	Mutagen: Reduces catalytic activity by half.
 109	V	I	Mutagen: Reduces catalytic activity by half.
 116	#	-	BINDING: Substrate.
 116	V	A	Mutagen: Reduces catalytic activity 10-fold.
 116	V	F	Mutagen: Slight increase of catalytic activity.
 173	V	N	in allele GSTM1B; dbSNP:rs1065411.
 210	V	T	in dbSNP:rs449856.

If features are provided, then a legend of feature symbols is provided
as well (disabled with C<--no-label>):

 =*:phosphorylation
 ==:active site
 =@:site
 =^:binding
 =!:metal binding

If the C<--lav> option is specified, domain and repeat features are
presented in a different format for the C<lav2plt.pl> program:

  >sp|P09488|GSTM1_HUMAN
  2	88	GST N-terminal.
  90	208	GST C-terminal.

C<ann_feats2.pl> is designed to be used by the B<FASTA> programs with
the C<-V \!ann_feats2.pl> option.  It can also be used with the lav2plt.pl
program with the C<--xA "\!ann_feats2.pl --lav"> or C<--yA "\!ann_feats2.pl --lav"> options.

=head1 AUTHOR

William R. Pearson, wrp@virginia.edu

=cut
