#!/usr/bin/env perl

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

################################################################
# clustal2fasta.pl 
################################################################
# clustal2fasta.pl takes a standard clustal format alignment file
# and produces the corresponding FASTA file.
#
################################################################

use warnings;
use strict;
use Pod::Usage;
use Getopt::Long;

my ($shelp, $help, $trim) = (0, 0);

GetOptions(
    "h|?" => \$shelp,
    "help" => \$help,
    );

pod2usage(1) if $shelp;
pod2usage(exitstatus => 0, verbose => 2) if $help;
unless (-f STDIN || -p STDIN || @ARGV) {
 pod2usage(1);
}

my @seq_ids = ();
my %msa = ();
    
# read the first line, first should not be blank
my $title = <>;

while (my $line = <>) {
  chomp $line;
  next unless ($line);
  next if ($line =~ m/^[\s:\*\+\.]+$/);   # skip conservation line

  my ($seq_id, $align) = split(/\s+/,$line);

  if (defined($msa{$seq_id})) {
    $msa{$seq_id} .= $align;
  }
  else {
    $msa{$seq_id} = $align;
    push @seq_ids, $seq_id;
  }
}

for my $seq_id ( @seq_ids ) {
  my $fmt_seq = $msa{$seq_id};
  $fmt_seq =~ s/(.{0,60})/$1\n/g;
  print ">$seq_id\n$fmt_seq";
}

__END__

=pod

=head1 NAME

 clustal2fasta.pl

=head1 SYNOPSIS

 clustal2fasta.pl clustal.msa

=head1 OPTIONS

 -h	short help
 --help include description


=head1 DESCRIPTION

C<clustal2fasta.pl> takes a Clustal format interleaved multiple
sequence alignment and produces the corresponding fasta format library.

=head1 AUTHOR

William R. Pearson, wrp@virginia.edu

=cut
