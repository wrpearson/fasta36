#!/usr/bin/env perl

################################################################
# copyright (c) 2014,2015 by William R. Pearson and The Rector &
# Visitors of the University of Virginia
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

# ann_exons_up_sql.pl gets an annotation file from fasta36 -V with a line of the form:

# gi|62822551|sp|P00502|GSTA1_RAT Glutathione S-transfer\n  (at least from pir1.lseg)
#
# it must:
# (1) read in the line
# (2) parse it to get the up_acc
# (3) return the tab delimited features

# this version can read feature2 uniprot features (acc/pos/end/label/value), but returns sorted start/end domains
# modified 18-Jan-2016 to produce annotation symbols consistent with ann_exons_up_www2.pl
# modified Dec 2018 to generate genomic coordinates with --gen_coord
# modified 3-Jan-2019 to merge sql and www (--www) access to exon coordinates

use warnings;
use strict;

use DBI;
use Getopt::Long;
use Pod::Usage;
use LWP::Simple;
use LWP::UserAgent;
use JSON qw(decode_json);

use vars qw($host $db $a_table $port $user $pass);

my %domains = ();
my $domain_cnt = 0;

my $hostname = `/bin/hostname`;

unless ($hostname =~ m/ebi/) {
  ($host, $db, $a_table, $port, $user, $pass)  = ("wrpxdb.its.virginia.edu", "uniprot", "annot2", 0, "web_user", "fasta_www");
#  $host = 'xdb';
}
else {
  ($host, $db, $a_table, $port, $user, $pass)  = ("mysql-pearson-prod", "up_db", "annot", 4124, "web_user", "fasta_www");
}

my ($lav, $gen_coord, $exon_label, $use_www, $shelp, $help) = (0,0,0,0,0,0);

my ($show_color) = (1);
my $color_sep_str = " :";
$color_sep_str = '~';

GetOptions(
    "gen_coord|gene_coord!" => \$gen_coord,
    "exon_label|label_exons!" => \$exon_label,
    "www!" => \$use_www,
    "host=s" => \$host,
    "db=s" => \$db,
    "user=s" => \$user,
    "password=s" => \$pass,
    "port=i" => \$port,
    "lav" => \$lav,
    "h|?" => \$shelp,
    "help" => \$help,
    );

pod2usage(1) if $shelp;
pod2usage(exitstatus => 0, verbose => 2) if $help;
pod2usage(1) unless (@ARGV || -p STDIN || -f STDIN);

my $connect = "dbi:mysql(AutoCommit=>1,RaiseError=>1):database=$db";
$connect .= ";host=$host" if $host;
$connect .= ";port=$port" if $port;

my $dbh = DBI->connect($connect,
		       $user,
		       $pass
		      ) or die $DBI::errstr;


my $get_annot_sub = \&get_annots;

my $ua = LWP::UserAgent->new(ssl_opts=>{verify_hostname => 0});
my $uniprot_url = 'https://www.ebi.ac.uk/proteins/api/coordinates/';
my $uniprot_suff = ".json";

if ($use_www) {
  $get_annot_sub = \&get_annots_up_www;
}


my $get_annots_id = $dbh->prepare(qq(select up_exons.* from up_exons join annot2 using(acc) where id=? order by ix));
my $get_annots_acc = $dbh->prepare(qq(select up_exons.* from up_exons where acc=? order by ix));
my $get_annots_refacc = $dbh->prepare(qq(select ref_acc, start, end, ix from up_exons join annot2 using(acc) where ref_acc=? order by ix));

my $get_annots_sql = $get_annots_acc;

my ($tmp, $gi, $sdb, $acc, $id, $use_acc);

# get the query
my ($query, $seq_len) =  @ARGV;
$seq_len = 0 unless defined($seq_len);

$query =~ s/^>// if ($query);

my @annots = ();

#if it's a file I can open, read and parse it
unless ($query && ($query =~ m/[\|:]/ ||
		   $query =~ m/^[NX]P_/ ||
		   $query =~ m/^[OPQ][0-9][A-Z0-9]{3}[0-9]|[A-NR-Z][0-9]([A-Z][A-Z0-9]{2}[0-9]){1,2}\s/)) {

  while (my $a_line = <>) {
    $a_line =~ s/^>//;
    chomp $a_line;
    push @annots, show_annots($a_line, $get_annot_sub, $use_www);
  }
}
else {
  push @annots, show_annots("$query\t$seq_len", $get_annot_sub, $use_www);
}

for my $seq_annot (@annots) {
  print ">",$seq_annot->{seq_info},"\n";
  for my $annot (@{$seq_annot->{list}}) {
    if (!$lav && $show_color && defined($domains{$annot->[-1]})) {
      $annot->[-1] .= $color_sep_str.$domains{$annot->[-1]};
    }
    print join("\t",@$annot),"\n";
  }
}

exit(0);

sub show_annots {
  my ($query_len, $get_annot_sub, $use_www) = @_;

  my ($annot_line, $seq_len) = split(/\t/,$query_len);

  my %annot_data = (seq_info=>$annot_line);

  if ($annot_line =~ m/^gi\|/) {
    $use_acc = 1;
    ($tmp, $gi, $sdb, $acc, $id) = split(/\|/,$annot_line);
  }
  elsif ($annot_line =~ m/^(SP|TR):(\w+) (\w+)/) {
    ($sdb, $id, $acc) = ($1,$2,$3);
    $use_acc = 1;
    $sdb = lc($sdb)
  }
  elsif ($annot_line =~ m/^(SP|TR):(\w+)/) {
    ($sdb, $id) = ($1,$2);
    $use_acc = 0;
    $sdb = lc($sdb)
  }
  elsif ($annot_line !~ m/\|/) {  # new NCBI swissprot format
    $use_acc =1;
    $sdb = 'sp';
    ($acc) = split(/\s+/,$annot_line);
  }
  else {
    $use_acc = 1;
    ($sdb, $acc, $id) = split(/\|/,$annot_line);
  }

  unless ($use_acc) {
    $get_annots_sql = $get_annots_id;
    $get_annots_sql->execute($id);
  }
  else {
    unless ($sdb =~ m/ref/) {
      $get_annots_sql = $get_annots_acc;
    } else {
      $get_annots_sql = $get_annots_refacc;
    }
    $acc =~ s/\.\d+$//;

    unless ($use_www) {
      $get_annots_sql->execute($acc);
    }
    else {
      $get_annots_sql = $acc;
    }
  }

  $annot_data{list} = $get_annot_sub->($get_annots_sql, $seq_len);

  return \%annot_data;
}

sub get_annots {
  my ($get_annots_sql, $seq_len) = @_;

  my @feats = ();

  while (my $exon_hr = $get_annots_sql->fetchrow_hashref()) {
    my $ix = $exon_hr->{ix};
    if ($lav) {
      push @feats, [$exon_hr->{start}, $exon_hr->{end}, "exon_$ix~$ix"];
    } else {
      my ($exon_info,$ex_info_start, $ex_info_end) = ("exon_$ix~$ix","","");
      if ($gen_coord) {
	if (defined($exon_hr->{g_start})) {
	  my $chr=$exon_hr->{chrom};
	  $chr = "unk" unless $chr;
	  if ($chr =~ m/^\d+$/ || $chr =~m/^[XYZ]+$/) {
	    $chr = "chr$chr";
	  }
	  $ex_info_start = sprintf("exon_%d::%s:%d",$ix, $chr, $exon_hr->{g_start});
	  $ex_info_end   = sprintf("exon_%d::%s:%d",$ix, $chr, $exon_hr->{g_end});
	  if ($exon_label) {
	    $exon_info = sprintf("exon_%d{%s:%d-%d}~%d",$ix, $chr, $exon_hr->{g_start}, $exon_hr->{g_end}, $ix);
	    push @feats, [$exon_hr->{start}, "-", $exon_hr->{end}, $exon_info];
	  } else {
	    push @feats, [$exon_hr->{start}, "-", $exon_hr->{end}, $exon_info];
	    push @feats, [$exon_hr->{start},'<','-',$ex_info_start];
	    push @feats, [$exon_hr->{end},'>','-',$ex_info_end];
	  }
	}
      } else {
	push @feats, [$exon_hr->{start}, "-", $exon_hr->{end}, $exon_info];
      }
    }
  }

  return \@feats;
}

sub get_annots_up_www {
  my ($acc, $seq_len) = @_;

  my @feats = ();

#  my  $exon_json = get_https($uniprot_url.$acc.$uniprot_suff);
  my  $exon_json = get($uniprot_url.$acc.$uniprot_suff);

  unless (!$exon_json || $exon_json =~ m/errorMessage/ || $exon_json =~ m/Can not find/) {
    return parse_json_up_exons($exon_json);
  }
  else {
    return ();
  }
}

sub parse_json_up_exons {
  my ($exon_json) = @_;

  my @exons = ();
  my @ex_coords = ();

  my $acc_exons = decode_json($exon_json);

  my $exon_num = 1;
  my $last_end = 0;
  my $last_phase = 0;

  my $chrom = $acc_exons->{'gnCoordinate'}[0]{'genomicLocation'}{'chromosome'};
  my $rev_strand = $acc_exons->{'gnCoordinate'}[0]{'genomicLocation'}{'reverseStrand'};

  for my $exon ( @{$acc_exons->{'gnCoordinate'}[0]{'genomicLocation'}{'exon'}} ) {
    my ($p_begin, $p_end) = ($exon->{'proteinLocation'}{'begin'}{'position'},$exon->{'proteinLocation'}{'end'}{'position'});
    my ($g_begin, $g_end) = ($exon->{'genomeLocation'}{'begin'}{'position'},$exon->{'genomeLocation'}{'end'}{'position'});

    my $this_phase = 0;
    if (defined($g_begin) && defined($g_end)) {
      $this_phase = ($g_end - $g_begin + 1) % 3;
    }

    if (!defined($p_begin) || !defined($p_end)) {
      $exon_num++;
      $last_phase = 0;
      next;
    }

    if ($p_end >= $p_begin) {
      if ($p_begin == $last_end) {
	if ($last_phase==2) {
	  $p_begin += 1;
	}
	elsif ($last_phase==1) {
	  $last_end -= 1;
	  $exons[-1]->{seq_end} -= 1;
	}
      }

      if ($p_begin <= $last_end && $p_end > $last_end) {
	$p_begin = $last_end+1;
      }
      $last_end = $p_end;
      $last_phase = $this_phase;

      my ($gs_begin, $gs_end) = ($g_begin, $g_end);
      if ($rev_strand) {
	($gs_begin, $gs_end) = ($g_end, $g_begin);
      }

      push @exons, {
		    ix=>$exon_num,
		    start=>$p_begin,
		    end=>$p_end,
		    g_start=>$gs_begin,
		    g_end=>$gs_end,
		    chrom=>$chrom,
		   };

      $exon_num++;
    }
  }

  # check for domain overlap, and resolve check for domain overlap
  # (possibly more than 2 domains), choosing the domain with the best
  # evalue

  my @ex_feats = ();

  for my $exon_hr (@exons) {
    my $ix = $exon_hr->{ix};
    if ($lav) {
      push @ex_feats, [$exon_hr->{start}, $exon_hr->{end}, "exon_$ix~$ix" ];
    }
    else {
      my ($exon_info,$ex_info_start, $ex_info_end) = ("exon_$ix~$ix","","");
      if ($gen_coord) {
	if (defined($exon_hr->{g_start})) {
	  my $chr=$exon_hr->{chrom};
	  $chr = "unk" unless $chr;
	  if ($chr =~ m/^\d+$/ || $chr =~m/^[XYZ]+$/) {
	    $chr = "chr$chr";
	  }
	  $ex_info_start = sprintf("exon_%d::%s:%d",$ix, $chr, $exon_hr->{g_start});
	  $ex_info_end   = sprintf("exon_%d::%s:%d",$ix, $chr, $exon_hr->{g_end});
	  if ($exon_label) {
	    $exon_info = sprintf("exon_%d{%s:%d-%d}~%d",$ix, $chr, $exon_hr->{g_start}, $exon_hr->{g_end},$ix);
	    push @ex_feats, [$exon_hr->{start}, "-", $exon_hr->{end}, $exon_info];
	  } else {
	    push @ex_feats, [$exon_hr->{start}, "-", $exon_hr->{end}, $exon_info];
	    push @ex_feats, [$exon_hr->{start},'<','-',$ex_info_start];
	    push @ex_feats, [$exon_hr->{end},'>','-',$ex_info_end];
	  }
	}
      } else {
	push @ex_feats, [$exon_hr->{start}, "-", $exon_hr->{end}, $exon_info];
      }
    }
  }
  return \@ex_feats;
}

sub get_https {
  my ($url) = @_;

  my $result = "";
  my $response = $ua->get($url);

  if ($response->is_success) {
    $result = $response->decoded_content;
  } else {
    $result = '';
  }
  return $result;
}



__END__

=pod

=head1 NAME

ann_exons_up_sql.pl

=head1 SYNOPSIS

 ann_exons_up_sql.pl --lav 'sp|P09488|GSTM1_NUMAN' | accession.file

=head1 OPTIONS

 -h	short help
 --help include description
 --gen_coord  -- provide genomic exon start/stop coordinates as features
 --lav  produce lav2plt.pl annotation format, only show domains/repeats
 --host, --user, --password, --port --db -- info for mysql database

=head1 DESCRIPTION

C<ann_exons_up_sql.pl> extracts exon location information from
a msyql database (default name, uniprot) built from EBI/proteins API data.

Given a command line argument that contains a sequence accession
(P09488) or identifier (GSTM1_HUMAN), the program looks up the
features available for that sequence and returns them in a
tab-delimited format:

 >sp|P09488|GSTM1_HUMAN
 1	-	12	exon_1~1
 13	-	38	exon_2~2
 39	-	59	exon_3~3
 60	-	87	exon_4~4
 88	-	120	exon_5~5
 121	-	152	exon_6~6
 153	-	189	exon_7~7
 190	-	218	exon_8~8

C<ann_exons_up_sql.pl --gen_coord 'sp|P09488|GSTM1_HUMAN'>also provides genomic coordinates:

>sp|P09488|GSTM1_HUMAN
1	-	12	exon_1~1
1	<	-	exon_1::chr1:109687874
12	>	-	exon_1::chr1:109687909
13	-	37	exon_2~2
13	<	-	exon_2::chr1:109688170
37	>	-	exon_2::chr1:109688245
38	-	59	exon_3~3
38	<	-	exon_3::chr1:109688673
59	>	-	exon_3::chr1:109688737
...
190	-	218	exon_8~8
190	<	-	exon_8::chr1:109693206
218	>	-	exon_8::chr1:109693292

C<ann_exons_up_sql.pl> is designed to be used by the B<FASTA> programs
with the C<-V \!ann_exons_up_sql.pl> option, or by the
C<annot_blast_btop.pl> script.  It can also be used with the
lav2plt.pl program with the C<--xA "\!ann_exons_up_sql.pl --lav"> or
C<--yA "\!ann_exons_up_sql.pl --lav"> options.

=head1 AUTHOR

William R. Pearson, wrp@virginia.edu

=cut
