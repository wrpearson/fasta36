#!/usr/bin/perl -w
#
# parse_m9.pl -- a simple script to parse fasta/fastx/ssearch -m 9 output and produce a simple set of results:
# >query_id<tab>len
# hit_acc<tab>len<tab>score<tab>bits<tab>expect<tab>f_id<tab>f_sim<tab>...
#
use strict;
use Getopt::Long;
use vars qw($e_cutoff $p_cutoff $head);

die "usage -- parse_m9.pl [--head] [--expect e_cut] [--percid p_cut] m9_out.file\n" unless @ARGV;

$e_cutoff = 10.0;
$p_cutoff = 0.0;
$head = 0;

GetOptions("expect=s" => \$e_cutoff,
	   "percid=s" => \$p_cutoff,
	   "head" => \$head,
	  );

my @hit_fields = ();
my @m9_fields = ();
my $first_hit = 1;

my $res_handle;
for my $s_res_file ( @ARGV ) {

  next unless open($res_handle, $s_res_file);

  while (my ($q_num, $query_descr, $query_len, $best_yes) = skip_to_results($res_handle)) {
    last unless $query_descr;

    unless ($best_yes) {
      <$res_handle>;	# skip >>><<<
# uncomment for queries with no hits
#    print ">$query_descr\t$query_len\n";
      next;
    }

    print ">$query_descr\t$query_len\n" unless $head;

    while (my $line = <$res_handle>) { # for each result
      last if $line =~ m/>>><<</;
      next if $line =~ m/^\+\-/;	# skip over HSPs
      chomp ($line);
      my ($left, $right) = split(/\t/,$line);
      my @fields = split(/\s+/,$left);
      my @afields = split(/\s+/,$right);

      my $evalue = $fields[-1];
      my $percid = $afields[0]*100.0;
      last if ($evalue > $e_cutoff && $percid < $p_cutoff);

      my $frame = "";
      my $l_len = $fields[-4];
      if ($fields[-4] =~ m/\[f|r\]/) {
	$l_len = $fields[-5];
	$frame = $fields[-4];
	if ($head) {unshift @hit_fields, "[fr]";}
      }

      if ($head && $first_hit) {
	unshift @hit_fields, qw(acc llen); 
	print "#" . join("\t",(@hit_fields, @m9_fields)) . "\n";
	print ">$query_descr\t$query_len\n";
	$first_hit = 0;
	$head = 0;
      }

      $l_len =~ s/\(//;
      $l_len =~ s/\)//;
      my ($l_db,$l_acc) = parse_descr($fields[0]);

      my @out_fields = ($l_acc, $l_len);
      if ($l_db) { unshift @out_fields, $l_db;}
      if ($frame) { push @out_fields, $frame;}
      print join("\t",(@out_fields, @fields[-3,-2,-1], @afields)) . "\n";
    }
  }
}

sub skip_to_results {
  my ($res_handle) = @_;
  my ($q_num, $query_desc, $query_len, $best_yes);

  while (my $line = <$res_handle>) {
    if ($line =~ m/^\s*(\d+)>>>(\S+)\s/) {
      ($q_num,$query_desc) = ($1,$2);
      ($query_len) = ($line =~ m/\s(\d+)\s\w+$/);
      goto have_query;
    }
    elsif ($line =~ m/>>>\/\/\//) {goto done;}
  }
  warn "EOF - no query\n";
 done:
  return "";

 have_query:
  while (my $line = <$res_handle>) {
    $best_yes = 0;
    if ($line =~ m/^The best scores are:/) {
      my ($left, $right) = split(/\t/,$line);
      if ($head) {
	my @afields = split(/\s+/,$left);
	@m9_fields = split(/\s+/,$right);
	@hit_fields = @afields[-3,-2,-1];
      }
      $best_yes = 1;
      last;
    }
    last if ($line =~ m/^!! No sequences/);
  }
  return ($q_num, $query_desc, $query_len, $best_yes);
}

sub parse_descr {
  my ($descr) = @_;

  my ($dummy, $gi, $db, $acc);

  if ($descr !~ m/\|/) {
    $db="";
    $acc=$descr;
  }
  elsif ($descr =~ m/gi\|/) {
    # has std gi|12345|ref|acc
    ($dummy, $gi, $db, $acc) = split(/\|/,$descr);
  }
  elsif ($descr =~ m/\d+\|\w+/) {
    # has std 12345|ref|acc from libtype=10
    ($gi, $db, $acc) = split(/\|/,$descr);
  }

  # remove version number
  $acc =~ s/\.\d+$//;

  return ($db, $acc);
}
