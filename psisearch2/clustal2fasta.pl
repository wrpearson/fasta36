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

################################################################
# clustal2fasta.pl 
################################################################
# clustal2fasta.pl takes a standard clustal format alignment file
# and produces the corresponding FASTA file.
#
# if --end_mask or --int_mask are set, then end or internal '-'s are converted to the query (first) sequence
# if --trim is set, then alignments beyond the beginning/end of the query sequence are trimmed
#
################################################################

use strict;
use Pod::Usage;
use Getopt::Long;

my ($shelp, $help, $trim) = (0, 0, 0);
my ($query_file, $sel_file, $bound_file_in, $bound_file_only, $bound_file_out, $masked_lib_out,$mask_type_end, $mask_type_int) = ("","","","","","","","");

my $query_lib_r = 0;

GetOptions(
    "query=s" => \$query_file,
    "query_file=s" => \$query_file,
    "masked_library_out=s" => \$masked_lib_out,
    "masked_lib_out=s" => \$masked_lib_out,
    "mask_lib_out=s" => \$masked_lib_out,
    "mask_out=s" => \$masked_lib_out,
    "end_mask_type=s" => \$mask_type_end,
    "end_mask=s" => \$mask_type_end,
    "int_mask_type=s" => \$mask_type_int,
    "int_mask=s" => \$mask_type_int,
    "h|?" => \$shelp,
    "help" => \$help,
    );

pod2usage(1) if $shelp;
pod2usage(exitstatus => 0, verbose => 2) if $help;
unless (-f STDIN || -p STDIN || @ARGV) {
 pod2usage(1);
}

################
# initialization
my @random_res = ();
if ($mask_type_end =~ m/^rand/i) {
  @random_res = init_random_res();
}

my @hit_list = ();

my %seq_bound = ();  # boundaries for each accession
my %acc_names = ();  # generate uniq s_seq_id names
my %multi_align = ();
my @multi_names = ();


################
# get sequence boundaries if available
#
my $seq_bound_hr = 0;
my @seq_bound_accs = ();

$seq_bound_hr = \%seq_bound;

if ($bound_file_in) {
  $seq_bound_hr = parse_bound_file($bound_file_in);
}
elsif ($bound_file_only) {
  $seq_bound_hr = parse_bound_file($bound_file_only);
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

sub parse_bound_file {
  my ($bound_file) = @_;

  open(my $qfd, $bound_file) || return 0;

  while (my $line = <$qfd>) {
    next if ($line =~ m/^#/);
    chomp $line;
    my @data = split(/\t/,$line);
    if (!defined($seq_bound{$data[0]})) {
      $seq_bound{$data[0]} = {start=>$data[1], end=>$data[2]};
      push @seq_bound_accs, $data[0];
    }
    else {
      warn "multiple boundaries for $data[0]";
    }
  }

  return \%seq_bound;
}

################
# init_random_res initializes a 1000 element array of amino acid residues with Robinson/Robinson frequencies

sub init_random_res {
  my @rr_res = qw(A R N D C Q E G H I L K M F P S T W Y V);
  my @rr_counts = (35155, 23105, 20212, 24161,  8669,  19208,
		   28354, 33229,  9906, 23161, 40625,  25872,
		   10101, 17367, 23435, 32070, 26311,   5990,
		   14488, 29012);
  my $rr_total = 450431;

  my $rr_seq = "";
  for (my $i=0; $i < 20; $i++) {
    $rr_seq = $rr_seq . $rr_res[$i] x int(1000.0 *$rr_counts[$i]/$rr_total + 0.5);
  }
  return split(//,$rr_seq);
}

__END__

=pod

=head1 NAME

 m89_btop_msa2.pl

=head1 SYNOPSIS

 m89_btop_msa2.pl --query_file query.fasta fasta_m8CB_output.file
 m89_btop_msa2.pl --query_file query.fasta blast_outfmt7_BTOP.output.file
 m89_btop_msa2.pl --query_file query.fasta [--m_format m9] fasta_m9C_output.file

=head1 OPTIONS

 -h	short help
 --help include description

 --query_file -- query sequence file
 --query      -- same as --query_file
 (only one sequence per file)

 --eval2 : "": use E()-value, "eval2": use E2()/eval2, "ave": use geom. mean

 --bound_file_in -- tab delimited accession<tab>start<tab>end that
                    specifies MSA boundaries WITHIN alignment.
                    Additional hits use alignment (or domain)
                    boundaries.

 --bound_file_only -- tab delimited accession<tab>start<tab>end that
                      specifies MSA boundaries WITHIN alignment.
                      Only sequences in --bound_file_only will be in the MSA.

 --bound_file_out -- "--bound_file" for next iteration of psisearch2

 --domain_bound  parse domain annotations (-V) from m9B file
 --domain

 --masked_lib_out -- FASTA format library of MSA sequences

 --int_mask_type = "query", "rand", "X", "none"
 --end_mask_type = "query", "rand", "X", "none"
specify the residues to be inserted into output library

=head1 DESCRIPTION

C<m89_btop_msa2.pl> takes a fasta36/ssearch36 -m 9B ouput file, which
includes a BTOP encoded alignment string, and produces the multiple
sequence alignment (MSA) implied by the query sequence, alignment
boundaries, and pairwise alignments.  The alignment does not allow
gaps in the query sequence, only in the subject sequences.

The C<--query_file> must be specified, and the query sequence is
provided as the first sequence in the MSA.

If a C<--bound_file> is provided, then the ends of the alignments are
reduced to the coordinates specified by the C<bound_file>.  In
addition, only sequences included in the C<bound_file> are included in
the MSA.

Output: A clustal-like interleaved multiple sequence alignment that
can be used as input (using the C<-in_msa> option) to C<psiblast>.

If an C<--masked_lib_out> filename.fasta is specified, a version of the MSA
in FASTA format is written to filename.fasta.  This file can be
converted to BLAST format (C<makeblastdb --in filename.fasta>) and
and the converted blast library can be used to rebuild a PSSM with
C<psiblast -num_iterations 2 -db filename.fasta -in_msa filename.msa
-out_pssm filename.asn_txt>.

The sequences in the C<--masked_lib_out> fasta file can be modifed
where gaps are present in the MSA as specified by the C<--int_mask_type> and C<--end_mask_type>
options. If no --int_mask_type/--end_mask_type is specified
("none"), then the subject sequence in the output library matches the
aligned part of the subject sequenc (gaps characters are deleted).  If
the --end/int_mask_type is "query", "rand", or "X", then either the
aligned query residue, a random residue, or an "X" substituted at
each gap, producing a library of subject sequences that may differ
from the original subject sequences.  These different options can be
used to force C<psiblast> to build a PSSM that more accurately
reflects the original C<ssearch36> alignment.

=head1 AUTHOR

William R. Pearson, wrp@virginia.edu

=cut
