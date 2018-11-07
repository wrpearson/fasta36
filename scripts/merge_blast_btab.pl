#!/usr/bin/env perl

################################################################
# copyright (c) 2018 by William R. Pearson and The Rector &
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

################################################################
# merge_blast_btab.pl --btab .btab file html_file
################################################################

use warnings;
use strict;
use Getopt::Long;
use Pod::Usage;
use URI::Encode qw(uri_encode);
use URI::Escape qw(uri_escape);

my ($btab_file, $have_qslen, $help, $shelp, $dom_info) = ("", 0, 0, 0, 0);
my ($plot_url) = ("");

GetOptions(
    "btab_file|btab=s" => \$btab_file,
    "have_qslen|have_sqlen!" => \$have_qslen,
    "domain_info|dom_info!" => \$dom_info,
    "plot_url=s"=> \$plot_url,
    "h|?" => \$shelp,
    "help" => \$help,
     );

pod2usage(1) if $shelp;
pod2usage(exitstatus => 0, verbose => 2) if $help;
unless (-f STDIN || -p STDIN || @ARGV) {
  pod2usage(1);
}

# require a btab file

# read it in, save structure as list/hash on accession (list more robust)
# what happens with multiple hits for same library -- need to add code
#

my @bl_fields = qw(q_seqid s_seqid percid alen mismatch gopen q_start q_end s_start s_end evalue bits score annot);

if ($have_qslen) {
  @bl_fields = qw(q_seqid q_len s_seqid s_len percid alen mismatch gopen q_start q_end s_start s_end evalue bits score annot);
}

if ($dom_info) {
  push @bl_fields, "dom_info";
}

my %tab_data = ();
my @sseq_ids = ();

unless ($btab_file) {
  die "--btab_file required"
}
else {
  # read in btab file
  open(my $fd, $btab_file) || die "cannot open $btab_file";

  while (my $line = <$fd>) {
    next if ($line =~ m/^#/);  # ignore comments
    chomp($line);
    my %a_data = ();
    @a_data{@bl_fields} = split(/\t/,$line);

    # here we should confirm that the sseqid is new.  If it is not, then add to a list.
    my $sseqid = $a_data{'s_seqid'};

    if (defined($tab_data{$sseqid})) {
      push @{$tab_data{$sseqid}}, \%a_data
    }
    else {
      $tab_data{$sseqid} = [ \%a_data ];
      push @sseq_ids, $sseqid;
    }
  }
}

# have the annotation data in %tab_data{} and @seq_ids
# read in the blastp html file and annotate it

my ($in_best, $in_align) = (0,0);
my ($best_ix, $align_ix, $hsp_ix) = (0,0,0);

while (my $line = <>) {
  chomp($line);
  unless ($line) {
    print "\n";
    next;
  }
  if ($line =~ m/^Sequences producing/) {
    $in_best = 1;
    $best_ix = 0;
    print "$line\n";
    next;
  }

  if ($in_best) {
    if ($line =~ /^>/) {
      $in_best = 0;
      $in_align = 1;
      $align_ix = 0;
      $hsp_ix = 0;
      # print out the first line
      print "$line\n";
      next;
    }
    else {
      $line = add_best($line, $tab_data{$sseq_ids[$best_ix]}->[0]);
      $best_ix++;
    }
  }

  if ($in_align) {
    if ($line =~ m/^\s+Score = \d+/) { # have Length= match, put out annotations if available
      my $regions_str = regions_to_str($tab_data{$sseq_ids[$align_ix]}->[$hsp_ix]);
      print $regions_str;

      if ($plot_url) {
	my $raw_dom_str = "";
	if ($dom_info) {
	  $raw_dom_str = dom_info_str($tab_data{$sseq_ids[$align_ix]}->[$hsp_ix]{'dom_info'});
	}

	my $plot_tag = plot_tag_str($plot_url, $tab_data{$sseq_ids[$align_ix]}->[$hsp_ix], $regions_str, $raw_dom_str);
	if ($plot_tag) {print $plot_tag,"\n";}
      }

      $hsp_ix++;

    }
    elsif ($line =~ m/^>/) {
      $align_ix++;
      $hsp_ix = 0;
    }
  }

  print "$line\n";
}

sub parse_annots {
  my ($annot_str) = @_;

  my @annot_list = ();

  unless ($annot_str && $annot_str =~ m/^\|/) {
    return \@annot_list;
  }

  my @annots = split('\|',$annot_str);
  shift @annots;

  for my $annot ( @annots ) {
    my %annot_data = ();
    next unless ($annot =~ m/^[XR][RX]/);
    my @a_fields = split(/;/,$annot);
    for my $f (@a_fields) {
      if ($f =~ m/^[XR][XR]/) {
	my @a2_f = split(':',$f);
	if ($a2_f[0] =~ m/^XR/) {
	  $annot_data{target} = 'subj';
	}
	else {
	  $annot_data{target} = 'query';
	}
	$annot_data{coord} = "$a2_f[1]:$a2_f[2]";
	$annot_data{score} = (split('=',$a2_f[3]))[1]
      }
      elsif ($f =~ m/(\w)=(.+)/) {
	$annot_data{$1} = $2;
      }
    }
    $annot_data{name} = $a_fields[-1];
    $annot_data{name} =~ s/^C=//;
    push @annot_list, \%annot_data;
  }
  return \@annot_list;
}

sub regions_to_str {
  my ($a_data_r) = @_;

  my $annot_ref = parse_annots($a_data_r->{annot});

  my $region_str = "";
  my $annot_str = "";

  for my $annot ( @{$annot_ref}) {
    if ($annot->{target} =~ m/^q/) {
      $region_str = "qRegion";
    }
    else {
      $region_str = " Region";
    }

    $annot_str .= sprintf "%s: %s : score=%d; bits=%.1f; Id=%.3f; Q=%.1f : %s\n", $region_str,
      @{$annot}{qw(coord score b I Q name)};
  }
  return $annot_str;
}

sub add_best {
  my ($line, $a_data) = @_;

  my $annot_str = '';

  my $annot_refs = parse_annots($a_data->{annot});

  for my $annot ( @$annot_refs) {
    if ($annot->{target} !~ m/^q/) {
      $annot_str .= $annot->{name} . ";"
    }
  }

  if ($annot_str) {
    return "$line  $annot_str";
  }
  else {
    return $line;
  }
}

sub plot_tag_str {

  my ($plot_script, $align_data_r, $regions_str, $doms_str) = @_;

  my $svg_pref =  q(<object type="image/svg+xml" );
  my $svg_post =  q( width="660" height="76" ></object>);

  #build argument string
  my %plt_args = ();
  @plt_args{qw(q_cstart l_cstart)} = (1, 1);
  @plt_args{qw(q_name q_cstop q_astart q_astop l_name l_cstop l_astart l_astop)} =
    @{$align_data_r}{qw(q_seqid q_len q_start q_end s_seqid s_len s_start s_end)};
  $plt_args{'regions'}= uri_escape(uri_encode($regions_str));
  if ($doms_str) {
    $plt_args{'doms'} = uri_encode($doms_str);
  }

  my $dom_info = ();

  my @args = map {"$_=$plt_args{$_}"} keys(%plt_args);

  return $svg_pref . qq( data="$plot_url?) . join('&amp;',@args) . '"' . $svg_post;
}

sub dom_info_str {
  my ($raw_dom_info) = @_;

  my $dom_str = "";
  
  unless ($raw_dom_info) { return "";}

  my @raw_doms = split('\|',$raw_dom_info);
  shift(@raw_doms);

  for my $dom ( @raw_doms ) {
    my $tmp_dom = $dom;
    $tmp_dom =~ s/^DX:/qDomain:\t/g;
    $tmp_dom =~ s/^XD:/lDomain:\t/g;
    $tmp_dom =~ s/;C=/\t/g;

    $dom_str .= "$tmp_dom\n";
  }

  return $dom_str;
}


__END__

=pod

=head1 NAME

merge_blast_btab.pl

=head1 SYNOPSIS

merge_blast_btab.pl --btab_file=result.b_tab result.html

=head1 OPTIONS

 -h	short help
 --help include description

 --btab_file|--btab file_name -- blast tabular output file with
   sub-alignment scoring

=head1 DESCRIPTION

C<merge_blast_btab.pl> merges the domain annotations and sub-alignment scoring from C<annot_blast_btop2.pl> blast tabular output file with a conventional blast result file.

The tab file is read and parsed, and then the subject/query seqid is used to
capture domain locations in the subject/query sequence.  If the domains
overlap the aligned region, the domain names are appended to the output.

=head1 AUTHOR

William R. Pearson, wrp@virginia.edu

=cut
