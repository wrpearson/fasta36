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

# ann_ipr_www.pl gets an annotation file from fasta36 -V with a line of the form:


# gi|62822551|sp|P00502|GSTA1_RAT Glutathione S-transfer\n  (at least from pir1.lseg)

# this version only annotates sequences known to InterPro
# and only provides domain information

# This script uses the dbfetch iprmc database, which REQUIRES a
# Uniprot Acc (not ID).  If an Acc is not provided, we must get an ACC
# first from the ID.

# SP:GSTM1_HUMAN P09488 218
#
# it must:
# (1) read in the line
# (2) parse it to get the up_acc
# (3) return the tab delimited domains
#

use strict;

use Getopt::Long;
use Pod::Usage;
use LWP::Simple;
## use IO::String;

# use dbfetch and IPRMC to get Interpro domain coordinates
# http://www.ebi.ac.uk/Tools/dbfetch/dbfetch?db=iprmc&id=gstm1_human&format=gff2&style=default&Retrieve=Retrieve

my $ipr_base = 'http://www.ebi.ac.uk/Tools/dbfetch/dbfetch?db=iprmc&id=';
my $gff_post = '&format=gff2&style=default&Retrieve=Retrieve';

################################################################
#
##gff-version 2
##Type Protein
# InterPro Matches for UniProtKB entries 
##source-version InterProMatches 49.0
##date 20-NOV-14

##sequence-region P09488 1 218
# P09488	InterProScan	region	99	190	1.1999999998684077E-49	.	.	Signature GENE3D G3DSA:1.20.1050.10 "G3DSA:1.20.1050.10" T ; InterPro IPR010987 "Glutathione S-transferase, C-terminal-like"
# P09488	InterProScan	region	2	98	5.800000000494973E-51	.	.	Signature GENE3D G3DSA:3.40.30.10 "G3DSA:3.40.30.10" T ; InterPro IPR012336 "Thioredoxin-like fold"
# P09488	InterProScan	region	105	189	3.900000000000007E-16	.	.	Signature PFAM PF00043 "GST_C" T ; InterPro IPR004046 "Glutathione S-transferase, C-terminal"
# P09488	InterProScan	region	4	82	7.299999999999985E-21	.	.	Signature PFAM PF02798 "GST_N" T ; InterPro IPR004045 "Glutathione S-transferase, N-terminal"
# P09488	InterProScan	region	31	43	1.1000015067164208E-25	.	.	Signature PRINTS PR01267 "GSTRNSFRASEM" T ; InterPro IPR003081 "Glutathione S-transferase, Mu class"
# P09488	InterProScan	region	44	56	1.1000015067164208E-25	.	.	Signature PRINTS PR01267 "GSTRNSFRASEM" T ; InterPro IPR003081 "Glutathione S-transferase, Mu class"
# P09488	InterProScan	region	87	98	1.1000015067164208E-25	.	.	Signature PRINTS PR01267 "GSTRNSFRASEM" T ; InterPro IPR003081 "Glutathione S-transferase, Mu class"
# P09488	InterProScan	region	139	152	1.1000015067164208E-25	.	.	Signature PRINTS PR01267 "GSTRNSFRASEM" T ; InterPro IPR003081 "Glutathione S-transferase, Mu class"
# P09488	InterProScan	region	1	88	0.0	.	.	Signature PROFILE PS50404 "GST_NTER" T ; InterPro IPR004045 "Glutathione S-transferase, N-terminal"
# P09488	InterProScan	region	90	208	0.0	.	.	Signature PROFILE PS50405 "GST_CTER" T ; InterPro IPR010987 "Glutathione S-transferase, C-terminal-like"
# P09488	InterProScan	region	1	217	0.0	.	.	Signature PANTHER PTHR11571 "PTHR11571" T
# P09488	InterProScan	region	1	217	0.0	.	.	Signature PANTHER PTHR11571:SF117 "PTHR11571:SF117" T
# P09488	InterProScan	region	86	217	8.190000000746436E-47	.	.	Signature SSF SSF47616 "SSF47616" T ; InterPro IPR010987 "Glutathione S-transferase, C-terminal-like"
# P09488	InterProScan	region	3	85	3.339999999911062E-23	.	.	Signature SSF SSF52833 "SSF52833" T ; InterPro IPR012336 "Thioredoxin-like fold"
###

my %domains = ();
my $domain_cnt = 0;

my $hostname = `/bin/hostname`;

my ($sstr, $lav, $neg_doms, $no_doms, $no_feats, $no_over, $data_file, $shelp, $help) = (0,0,0,0,1,0,0,0,0);
my $dom_dbs = "PFAM+PROFILE+GENE3D";

my ($min_nodom) = (10);

my $color_sep_str = " :";
$color_sep_str = '~';

GetOptions(
    "lav" => \$lav,
    "no-over" => \$no_over,
	   "no_doms" => \$no_doms,
	   "no-doms" => \$no_doms,
	   "nodoms" => \$no_doms,
    	   "dom_dbs:s" => \$dom_dbs,	# PF, PS, 
    	   "dbs:s" => \$dom_dbs,
	   "neg" => \$neg_doms,
	   "neg_doms" => \$neg_doms,
	   "neg-doms" => \$neg_doms,
	   "negdoms" => \$neg_doms,
	   "min_nodom=i" => \$min_nodom,
	   "no_feats" => \$no_feats,
	   "no-feats" => \$no_feats,
	   "nofeats" => \$no_feats,
	   "data:s" => \$data_file,
	   "sstr" => \$sstr,
	   "h|?" => \$shelp,
	   "help" => \$help,
	  );

pod2usage(1) if $shelp;
pod2usage(exitstatus => 0, verbose => 2) if $help;
pod2usage(1) unless (-p STDIN || -f STDIN || @ARGV || $data_file);

my @feat_keys = qw(catalytic_residue posttranslation_modification binding_motif metal_contact
		   polypeptide_region mutated_variant_site natural_variant_site);

my %feats_text = ();
@feats_text{@feat_keys} = ('Active site', '', 'Substrate binding', 'Metal binding', 'Site', '','');

my %feats_label;
@feats_label{@feat_keys} = ('Active site', 'Modified', 'Substrate binding', 'Metal binding', 'Site', '','');

my @feat_vals = ( '=','*','#','^','@','V','V');

my %annot_types = ();

my $get_annot_sub = \&iprmc_annots;
if ($lav) {
  $no_feats = 1;
}

if ($dom_dbs) {
  my @dom_db_list = split(/\+/,$dom_dbs);

  for my $dom_db (@dom_db_list) {
    $annot_types{$dom_db} = $dom_db;
  }
}

if ($neg_doms) {
  $domains{'NODOM'}=0;
}

my ($tmp, $gi, $sdb, $acc, $id, $use_acc);

unless ($no_feats || $sstr) {
  for my $i ( 0 .. $#feat_keys) {
    next unless $feats_label{$feat_keys[$i]};
    print "=",$feat_vals[$i],":",$feats_label{$feat_keys[$i]},"\n";
  }
}

# get the query
my ($query, $seq_len) =  @ARGV;
$seq_len = 0 unless defined($seq_len);

$query =~ s/^>// if ($query);

my @annots = ();

#if it's a file I can open, read and parse it

unless ($data_file) {
  unless ($query && $query =~ m/[\|:]/) {

    while (my $a_line = <>) {
      $a_line =~ s/^>//;
      chomp $a_line;
      push @annots, lwp_annots($a_line, $get_annot_sub);
    }
  } else {
    push @annots, lwp_annots("$query\t$seq_len", $get_annot_sub);
  }
} else {   # just read the data from a file, give to $get_annot_sub().
  my %annot_data = (seq_info => ">$data_file DATA");

  open(DATA_IN, $data_file) || die "Cannot read $data_file";

  my $lwp_data = "";
  while (<DATA_IN>) {
    $lwp_data .= $_;
  }

  $annot_data{list} = $get_annot_sub->(\%annot_types, $lwp_data,0);

  push @annots, \%annot_data;
}


for my $seq_annot (@annots) {
  print ">",$seq_annot->{seq_info},"\n";
  for my $annot (@{$seq_annot->{list}}) {
    if (!$lav && defined($domains{$annot->[-1]})) {
      $annot->[-2] .= $color_sep_str.$domains{$annot->[-1]};
    }
    print join("\t",@{$annot}[0..3]),"\n";
  }
}

exit(0);

sub lwp_annots {
  my ($query_len, $get_annot_sub) = @_;

  my ($annot_line, $seq_len) = split(/\t/,$query_len);

  my %annot_data = (seq_info=>$annot_line);

  if ($annot_line =~ m/^gi\|/) {
    ($tmp, $gi, $sdb, $acc, $id) = split(/\|/,$annot_line);
  } elsif ($annot_line =~ m/^(SP|TR):(\w+)/) {
    $sdb = lc($1);
    $id = $2;
#     $acc = $2;
  } elsif ($annot_line =~ m/^(UR\d{3}:UniRef\d{2})_(\w+)/) {
    $sdb = lc($1);
    $id = $2;
#    $acc = $2;
  } else {
    ($sdb, $acc, $id) = split(/\|/,$annot_line);
  }

  $acc =~ s/\.\d+// if ($acc);

  $annot_data{list} = [];
  my $lwp_domains = "";

  if ($acc && ($acc =~ m/^[A-Z][0-9][A-Z0-9]{3}[0-9]/)) {
    $lwp_domains = get($ipr_base . $acc . $gff_post);
  } elsif ($id && ($id =~ m/^\w+$/)) {
    $lwp_domains = get($ipr_base . $id . $gff_post);
  }

  if ($lwp_domains && ($lwp_domains !~ /ERROR/)) {
    $annot_data{list} = $get_annot_sub->(\%annot_types, $lwp_domains, $seq_len);
  }

  return \%annot_data;
}

# parses www.uniprot.org gff feature table
sub iprmc_annots {
  my ($annot_types, $annot_data, $seq_len) = @_;

  my ($acc, $pos, $end, $label, $value, $comment, $len);
  my ($seq_acc, $seq_start, $seq_end, $tmp);

  $seq_len = 0;

  my @feats2 = (); # domains with start/stop, for checking overlap, adding negative
  my @sites = ();  # sites with one position

  my @gff_lines = split(/\n/m,$annot_data);

  while (my $gff_line = shift @gff_lines) {
    chomp $gff_line;
    if ($gff_line =~ m/^#sequence-region/) {
      my @fields = split($gff_line, /\s+/);
      $seq_end = $fields[-1];
      last;
    }
  }

  while (my $gff_line = shift(@gff_lines)) {
    next if ($gff_line =~ m/^#/);
    chomp($gff_line);

    my @gff_line_arr = split(/\t/,$gff_line);
    ($acc, $pos, $end, $comment) = @gff_line_arr[(0,3,4,-1)];

    # parse the comment to get signature (domain_db), domain_db_acc, interpro_acc, description
    my ($domain_info, $dom_acc) = parse_ipr_comment($comment);

    next unless $domain_info;

    push @feats2, [$pos, "-", $end, $domain_info, $dom_acc];

    $value = '' unless $value;
    #	print join("\t",($pos, $annot_types->{$label})),"\n";
    #	print join("\t",($pos, $annot_types->{$label}, "-", "$label: $value")),"\n";
  }

  @feats2 = sort { $a->[0] <=> $b->[0] } @feats2;

  if ($no_over) {
    # check for containment
    my $have_contained = 0;
    my $last_container = 0;
    for (my $i=1; $i < scalar(@feats2); $i++) {
      if ($feats2[$i]->[0] >= $feats2[$last_container]->[0] && $feats2[$i]->[2] <= $feats2[$last_container]->[2]) {
	$feats2[$i]->[1] = 'Delete';
	$have_contained = 1;
      } else {
	$last_container=$i;
      }
    }

    if ($have_contained) {
      @feats2 = grep { $_->[1] !~ /Delete/ } @feats2;
    }

    # ensure that domains do not overlap
    for (my $i=1; $i < scalar(@feats2); $i++) {
      my $diff = $feats2[$i-1]->[2] - $feats2[$i]->[0];
      if ($diff >= 0) {
	$feats2[$i-1]->[2] = $feats2[$i]->[0]+ int($diff/2);
	$feats2[$i]->[0] = $feats2[$i-1]->[2] + 1;
      }
    }
  }

  my @n_feats2 = ();

  if ($neg_doms) {
    my $last_end = 0;
    for my $feat ( @feats2 ) {
      if ($feat->[0] - $last_end > $min_nodom) {
	push @n_feats2, [$last_end+1, "-", $feat->[0]-1, "NODOM", ""];
      }
      $last_end = $feat->[2];
    }
    if ($seq_len - $last_end > $min_nodom) {
      push @n_feats2, [$last_end+1, "-", $seq_len, "NODOM", ""];
    }
  }

  my @feats = ();

  for my $feat (@feats2, @n_feats2) {
    if (!$lav)  {
      push @feats, [$feat->[0], '-', $feat->[2], $feat->[-2], $feat->[-1] ];
#      push @feats, [$feat->[2], ']', '-', ""];
    }
    else {
      push @feats, [$feat->[0], $feat->[2], $feat->[-1]];
    }
  }

  @feats = sort { $a->[0] <=> $b->[0] } (@sites, @feats);

  # now that domains are sorted, give them names
  for my $feat ( @feats ) {
    $feat->[-2] = domain_name($feat->[-2],$feat->[-1]);
  }

  return \@feats;
}

sub parse_ipr_comment {
  my ($comment_str) = @_;

  my @comments = split(/\s+;\s+/,$comment_str);
  my @comment_info = ();
  my $ipr_info = "";
  $comments[0] =~ s/^Signature\s+//;

  return ("","") unless @comments;

  for my $comment (@comments) {
    my %ipr_data = ();
    @ipr_data{qw(db acc descr)} = ($comment =~ m/(\S+)\s+(\S+)\s+"([^"]+)"/);
    return ("","") if $ipr_data{db} =~ m/(PRINTS|PROSITE)/i;
    $ipr_data{descr} =~ s/\s+/_/g;
    push @comment_info, \%ipr_data;
  }

  my $primary_acc = $comment_info[0]->{acc};

  for my $comment (@comment_info) {
    if ($comment->{db} =~ m/InterPro/) {
      return ("$primary_acc:".$comment->{descr},$comment->{acc});
    }
  }
  return ("","");
}

# domain name takes a uniprot domain label, removes comments ( ;
# truncated) and numbers and returns a canonical form. Thus:
# Cortactin 6.
# Cortactin 7; truncated.
# becomes "Cortactin"
#

sub domain_name {

  my ($value, $ipr_acc) = @_;

  $value = 'UnDef' unless $value;

  $value =~ s/;.*$//;
  $value =~ s/\.\s*$//;
  $value =~ s/\s+\d+$//;
  if (!defined($domains{$ipr_acc})) {
    $domain_cnt++;
    $domains{$ipr_acc} = $domain_cnt;
  }
  return $value;
}

__END__

=pod

=head1 NAME

ann_ipr_www.pl

=head1 SYNOPSIS

 ann_ipr_www.pl --no_doms --no_feats --lav 'sp|P09488|GSTM1_NUMAN' | accession.file

=head1 OPTIONS

 -h	short help
 --help include description
 --no-doms  do not show domain boundaries (domains are always shown with --lav)
 --no-feats do not show feature (variants, active sites, phospho-sites)
 --lav  produce lav2plt.pl annotation format, only show domains/repeats

 --neg-doms,  -- report domains between annotated domains as NODOM
                 (also --neg, --neg_doms)
 --min_nodom=10  -- minimum length between domains for NODOM

 --host, --user, --password, --port --db -- info for mysql database

=head1 DESCRIPTION

C<ann_ipr_www.pl> extracts feature, domain, and repeat
information from the Uniprot DAS server through an XSLT transation
provided by http://www.ebi.ac.uk/Tools/dbfetch/dbfetch/uniprotkb.
This server provides GFF descriptions of Uniprot entries, with most of
the information provided in UniProt feature tables.

C<ann_ipr_www.pl> is an alternative to C<ann_pfam.pl> and
C<ann_pfam.pl> that does not require a local MySQL copy of Pfam.

Given a command line argument that contains a sequence accession
(P09488), the program looks up the domains available for that
sequence and returns them in a tab-delimited format:

>sp|P09488|GSTM1_HUMAN
2	-	88	GST N-terminal :1
90	-	208	GST C-terminal :2

If the C<--lav> option is specified, domain and repeat domains are
presented in a different format for the C<lav2plt.pl> program:

  >sp|P09488|GSTM1_HUMAN
  2	88	GST N-terminal.
  90	208	GST C-terminal.

C<ann_ipr_www.pl> is designed to be used by the B<FASTA> programs with
the C<-V \!ann_ipr_www.pl> option.  It can also be used with the lav2plt.pl
program with the C<--xA "\!ann_ipr_www.pl --lav"> or C<--yA "\!ann_ipr_www.pl --lav"> options.

=head1 AUTHOR

William R. Pearson, wrp@virginia.edu

=cut
