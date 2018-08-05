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

# summ_domain_ident.pl takes a -m 8CC result file with query-annotated
# domains, and produces a tab-delimited summary of identity across the domains
# parse:
# sp|P09488|GSTM1_HUMAN	gi|121735|sp|P09488.3|GSTM1_HUMAN	100.00	218	0	0	1	218	1	218	2.9e-113	408.2	218M	|RX:1-12:1-12:s=64;b=25.0;I=1.000;Q=47.5;C=exon_1|RX:13-37:13-37:s=128;b=49.9;I=1.000;Q=121.4;C=exon_2|RX:38-59:38-59:s=125;b=48.7;I=1.000;Q=117.9;C=exon_3|RX:60-86:60-86:s=145;b=56.5;I=1.000;Q=141.0;C=exon_4|RX:87-120:87-120:s=185;b=72.1;I=1.000;Q=187.2;C=exon_5|RX:121-152:121-152:s=174;b=67.8;I=1.000;Q=174.5;C=exon_6|RX:153-189:153-189:s=197;b=76.8;I=1.000;Q=201.0;C=exon_7|RX:190-218:190-218:s=151;b=58.9;I=1.000;Q=147.9;C=exon_8

use warnings;
use strict;
use Getopt::Long;
use Pod::Usage;

my ($shelp, $help) = (0, 0);

GetOptions(
    "h|?" => \$shelp,
    "help" => \$help,
    );

pod2usage(1) if $shelp;
pod2usage(exitstatus => 0, verbose => 2) if $help;
pod2usage(1) unless @ARGV;

my $first_line =1;

my @a_field_names = qw( score_field bits id qval comment );
my @domain_names = ();

while (my $line = <>) {
  next if $line =~ m/^#/;
  chomp($line);
  
  # get last (annotation) field

  my @fields = split(/\t/,$line);
  $fields[1] =~ s/^gi\|\d+\|//;
  $fields[1] =~ s/\.\d+\|/\|/;

  my @annots = split(/\|/,$fields[-1]);
  shift @annots; # first is blank

  # $annots[...]:
  # RX:1-12:1-12:s=64;b=25.0;I=1.000;Q=47.5;C=exon_1

  my %dom_ids = ();

  for my $annot (@annots) {
    my %a_fields = ();
    @a_fields{@a_field_names} = split(/;/,$annot);
    $a_fields{'id'} =~ s/^I=//;
    $a_fields{'comment'} =~ s/^C=//;

    $dom_ids{$a_fields{'comment'}} = $a_fields{'id'};

    if ($first_line) {
      push @domain_names, $a_fields{'comment'};
    }
  }
  if ($first_line) {
    print join("\t",("subj_acc          ","ident",@domain_names)),"\n";
    $first_line = 0;
  }

  for my $dom ( @domain_names ) {
    if (defined($dom_ids{$dom})) {
      if (100.0 * $dom_ids{$dom} >= $fields[2]) {
	$dom_ids{$dom} .= '+';
      }
      else {
	$dom_ids{$dom} .= '-';
      }
    }
    else {
      $dom_ids{$dom} = '';
    }
  }

  print join("\t",($fields[1],$fields[2],@dom_ids{@domain_names})),"\n";
}
