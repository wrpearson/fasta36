#!/usr/bin/perl -w

################################################################
# copyright (c) 2014,2015 by William R. Pearson and The Rector &
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

## modified 29-Sept-2016 to use EBI/proteins JSON URL:
## http://www.ebi.ac.uk/proteins/api/features/p12345

# ann_feats_up_www2.pl gets an annotation file from fasta36 -V with a line of the form:

# SP:GSTM1_HUMAN P09488 218
#
# it must:
# (1) read in the line
# (2) parse it to get the up_acc
# (3) return the tab delimited features
#

# this version can read feature2 uniprot features (acc/pos/end/label/value), but returns sorted start/end domains

use strict;

use Getopt::Long;
use Pod::Usage;
use LWP::Simple;
use JSON qw(decode_json);

## use IO::String;

my $up_base = 'http://www.ebi.ac.uk/proteins/api/features';

my %domains = ();
my $domain_cnt = 0;

my $hostname = `/bin/hostname`;

my ($sstr, $lav, $neg_doms, $no_doms, $no_feats, $no_vars, $no_over, $data_file, $shelp, $help) = (0,0,0,0,0,0,0,0,0,0);
my ($min_nodom) = (10);

my ($show_color) = (1);
my $color_sep_str = " :";
$color_sep_str = '~';

GetOptions(
    "lav" => \$lav,
    "no-over" => \$no_over,
	   "no_doms" => \$no_doms,
	   "no-doms" => \$no_doms,
	   "nodoms" => \$no_doms,
	   "no_vars" => \$no_vars,
	   "no-vars" => \$no_vars,
	   "novars" => \$no_vars,
	   "neg" => \$neg_doms,
	   "neg_doms" => \$neg_doms,
	   "neg-doms" => \$neg_doms,
	   "negdoms" => \$neg_doms,
	   "min_nodom=i" => \$min_nodom,
	   "no_feats" => \$no_feats,
	   "no-feats" => \$no_feats,
	   "nofeats" => \$no_feats,
	   "data:s" => \$data_file,
    	   "color!" => \$show_color,
	   "sstr" => \$sstr,
	   "h|?" => \$shelp,
	   "help" => \$help,
	  );

pod2usage(1) if $shelp;
pod2usage(exitstatus => 0, verbose => 2) if $help;
pod2usage(1) unless @ARGV || $data_file || -p STDIN || -f STDIN;

my @feat_keys = qw( ACT_SITE MOD_RES BINDING METAL SITE );
my @feat_vals = ( '=','*','#','^','@');
my @feat_names = ('Active site', 'Modified', 'Binding', 'Metal binding', 'Site');

unless ($no_vars) {
  push @feat_keys, qw(MUTAGEN VARIANT);
  push @feat_vals, ('V','V',);
  push @feat_names, ("","",);
}

my %feats_text = ();
@feats_text{@feat_keys} = @feat_names;

my %feats_label;
@feats_label{@feat_keys} = @feat_names;

my @dom_keys = qw( DOMAIN REPEAT );
my @dom_vals = ( [ '[', ']'],[ '[', ']']);

my @ssr_keys = qw(STRAND HELIX);
my @ssr_vals = ( [ '[', ']'], [ '[', ']']);

my %annot_types = ();

my $get_annot_sub = \&json_annots;
if ($lav) {
  $no_feats = 1;
}

if ($sstr) {
  @annot_types{@ssr_keys} = @ssr_vals;
} else {
  @annot_types{@feat_keys} = @feat_vals unless ($no_feats);
  @annot_types{@dom_keys} = @dom_vals unless ($no_doms);
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
  unless ($query && $query =~ m/[\|:]/ ) {

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
    if (!$lav && $show_color && defined($domains{$annot->[-1]})) {
      $annot->[-1] .= $color_sep_str.$domains{$annot->[-1]};
    }
    print join("\t",@$annot),"\n";
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
  } elsif ($annot_line =~ m/\|/) {
    ($sdb, $acc, $id) = split(/\|/,$annot_line);
  }
  else {
    ($acc) = ($annot_line =~ m/^(\S+)/);
  }

  $acc =~ s/\.\d+// if ($acc);

  $annot_data{list} = [];
  my $lwp_features = "";

  if ($acc && ($acc =~ m/^[A-Z][0-9][A-Z0-9]{3}[0-9]/)) {
    $lwp_features = get("$up_base/$acc");
  }
#  elsif ($id && ($id =~ m/^\w+$/)) {
#    $lwp_features = get("$up_base/$id/$gff_post");
#  }

  if ($lwp_features && ($lwp_features !~ /ERROR/)) {
    my $annot_json = decode_json($lwp_features);
    $annot_data{list} = $get_annot_sub->(\%annot_types, $annot_json, $seq_len);
  }

  return \%annot_data;
}

####
# parses www.ebi.ac.uk/uniprot/api json
#
sub json_annots {
  my ($annot_types, $json_ref, $seq_len) = @_;

  my ($acc, $pos, $end, $label, $value, $comment, $len);

  $seq_len = 0;

  my @feats2 = (); # features with start/stop, for checking overlap, adding negative
  my @sites = ();  # sites with one position

  my ($seq_str, $seq_acc, $seq_id)  = @{$json_ref}{qw(sequence accession entryName)};
  $seq_len = length($seq_str);

  for my $feat ( @{$json_ref->{features}} ) {
    if ($annot_types->{$feat->{type}}) {

      my ($label, $pos, $end, $value)  = @{$feat}{qw(type begin end description)};

      $pos =~ s/[<>]//g;
      $end =~ s/[<>]//g;

      if ($label =~ m/DOMAIN/ || $label =~ m/REPEAT/) {
	$value = domain_name($label,$value);
	push @feats2, [$pos, "-", $end, $value];
      } elsif ($label =~ m/HELIX/) {
	push @feats2, [$pos, "-", $end, $label];
      } elsif ($label =~ m/STRAND/) {
	push @feats2, [$pos, "-", $end, $label];
      } elsif ($label =~ m/VARIANT/ || $label =~ m/MUTAGEN/) {
	push @sites, [$pos, $annot_types->{$label}, $feat->{alternativeSequence}, $value];
      } 
      else {
	next unless ($pos == $end);
	if ($feats_text{$label}) {
	  my $info = $feats_text{$label};
	  if ($value) {
	    $info .= ": $value";
	  }
	  push @sites, [$pos, $annot_types->{$label}, "-", $info];
	} else {
	  push @sites, [$pos, $annot_types->{$label}, "-", $value];
	}
      }
    }
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
	push @n_feats2, [$last_end+1, "-", $feat->[0]-1, "NODOM"];
      }
      $last_end = $feat->[2];
    }
    if ($seq_len - $last_end > $min_nodom) {
      push @n_feats2, [$last_end+1, "-", $seq_len, "NODOM"];
    }
  }

  my @feats = ();
  for my $feat (@feats2, @n_feats2) {
    if (!$lav)  {
      push @feats, [$feat->[0], '-', $feat->[2], $feat->[-1] ];
#      push @feats, [$feat->[2], ']', '-', ""];
    }
    else {
      push @feats, [$feat->[0], $feat->[2], $feat->[-1]];
    }
  }

  @feats = sort { $a->[0] <=> $b->[0] } (@sites, @feats);

  return \@feats;
}

# domain name takes a uniprot domain label, removes comments ( ;
# truncated) and numbers and returns a canonical form. Thus:
# Cortactin 6.
# Cortactin 7; truncated.
# becomes "Cortactin"
#

sub domain_name {

  my ($label, $value) = @_;

  $value = 'UnDef' unless $value;

  $value =~ s/ /_/g;

  if ($label =~ /Domain|Repeat/i) {
    $value =~ s/;.*$//;
    $value =~ s/\.\s*$//;
    $value =~ s/\s+\d+$//;
    if (!defined($domains{$value})) {
      $domain_cnt++;
      $domains{$value} = $domain_cnt;
    }
    return $value;
  }
  else {
    return $value;
  }
}



__END__

=pod

=head1 NAME

ann_feats_up_www2.pl

=head1 SYNOPSIS

 ann_feats_up_www2.pl --no_doms --no_feats --lav 'sp|P09488|GSTM1_NUMAN' | accession.file

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

C<ann_feats_up_www2.pl> extracts feature, domain, and repeat
information from the Uniprot DAS server through an XSLT transation
provided by http://www.ebi.ac.uk/Tools/dbfetch/dbfetch/uniprotkb.
This server provides GFF descriptions of Uniprot entries, with most of
the information provided in UniProt feature tables.

C<ann_feats_up_www2.pl> is an alternative to C<ann_pfam.pl> and
C<ann_pfam.pl> that does not require a local MySQL copy of Pfam.

Given a command line argument that contains a sequence accession
(P09488), the program looks up the features available for that
sequence and returns them in a tab-delimited format:

>sp|P09488|GSTM1_HUMAN
2	[	-	GST N-terminal :1
7	V	F	Mutagen: Reduces catalytic activity 100- fold.
23	*	-	MOD_RES: Phosphotyrosine (By similarity).
33	*	-	MOD_RES: Phosphotyrosine (By similarity).
34	*	-	MOD_RES: Phosphothreonine (By similarity).
88	]	-	
90	[	-	GST C-terminal :2
108	V	Q	Mutagen: Reduces catalytic activity by half.
108	V	S	Mutagen: Changes the properties of the enzyme toward some substrates.
109	V	I	Mutagen: Reduces catalytic activity by half.
116	#	-	BINDING: Substrate.
116	V	A	Mutagen: Reduces catalytic activity 10-fold.
116	V	F	Mutagen: Slight increase of catalytic activity.
173	V	N	in allele GSTM1B; dbSNP:rs1065411.
208	]	-	
210	V	T	in dbSNP:rs449856.

If features are provided, then a legend of feature symbols is provided
as well:

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

C<ann_feats_up_www2.pl> is designed to be used by the B<FASTA> programs with
the C<-V \!ann_feats_up_www2.pl> option.  It can also be used with the lav2plt.pl
program with the C<--xA "\!ann_feats_up_www2.pl --lav"> or C<--yA "\!ann_feats_up_www2.pl --lav"> options.

=head1 AUTHOR

William R. Pearson, wrp@virginia.edu

=cut
