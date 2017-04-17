#!/usr/bin/perl -w

use strict;
use Getopt::Long;
use Pod::Usage;

my ($window, $insert, $shelp, $help, $n_shuff) = (20, 1, 0, 0,1);

GetOptions("window=i" => \$window,
	   "insert=i" => \$insert,
	   "n=i" => \$n_shuff,
	   "h|?" => \$shelp,
	   "help" => \$help,

    );

pod2usage(1) if $shelp;
pod2usage(exitstatus => 0, verbose => 2) if $help;

my ($seq, $header) = ("","");

while (my $line = <>) {
  chomp($line);
  if ($line =~ /^>/) {
    if ($seq) { process_seq($header, $seq, $window, $insert, $n_shuff);}
    $header = $line;
    $seq = "";
  }
  else {
    $seq .= $line;
  }
}

if ($seq) { process_seq($header, $seq, $window, $insert, $n_shuff);}

exit(0);

sub process_seq {
  my ($header, $seq, $window, $insert, $n_shuff) = @_;

  # remove non amino-acids
  $seq =~ s/[^A-Za-z]//g;
  my $seq_len = length($seq);

  for (my $shuff_cnt = 0; $shuff_cnt < $n_shuff; $shuff_cnt++) {

    my $shuff_seq = win_shuffle($seq, $window);
    my $left_sseq = substr($shuff_seq, 0, ($seq_len+1)/2);
    my $right_sseq = substr($shuff_seq, ($seq_len+1)/2, $seq_len - ($seq_len+1)/2 +1);

    my $embed_seq = $left_sseq;
    if ($insert) {
      $embed_seq .= $seq;
    }
    $embed_seq .= $right_sseq;

    my ($acc, $descr) = ($header =~ m/^>(\S+)\s*(.*)$/);

    if ($insert) {
      if ($n_shuff > 1) {
	printf ">%s_%d e:%d-%d %s\n",$acc,$shuff_cnt,length($left_sseq)+1,length($left_sseq) + $seq_len,$descr;
      }
      else {
	printf ">%s e:%d-%d %s\n",$acc,length($left_sseq)+1,length($left_sseq) + $seq_len,$descr;
      }
    } else {
      if ($n_shuff > 1) {
	printf ">%s_shuff_%d %s\n",$acc,$shuff_cnt, $descr;
      }
      else {
	printf ">%s_shuff %s\n",$acc, $descr;
      }
    }

    $embed_seq =~ s/(.{60})/$1\n/g;

    print "$embed_seq\n";
  }
}

my $random_seed;

sub win_shuffle {
  my ($seq, $win) = @_;

  # break sequence into $win len pieces
  my @seq_arr = ($seq =~ m/(.{1,$win})/g);

  # shuffle the subsets
  for (my $j = 0; $j < @seq_arr; $j++) {
    my @subs_arr = split(//,$seq_arr[$j]);
    fy_shuffle(\@subs_arr);
    $seq_arr[$j] = join("",@subs_arr);
  }

  # now shuffle the window order
  fy_shuffle(\@seq_arr);
  # and put it back together
  return join("",@seq_arr);
}

# fy_shuffle array_ref
sub fy_shuffle {
  my $arr = shift;

  die "fy_shuffle (array_ref)" unless (ref($arr) eq 'ARRAY');

  return unless @$arr;
  my $i = scalar(@$arr)-1;

  while ($i > 0) {
    my $is = int(rand($i));
    ($arr->[$i],$arr->[$is]) = ($arr->[$is],$arr->[$i]);
    $i--;
  }
}

__END__

=pod

=head1 NAME

shuffle_embed.pl

=head1 SYNOPSIS

shuffle_embed.pl --n=1 --insert=1 --window=20 file.seq > file.shuff_emb

=head1 OPTIONS

 -h	    short help
 --help     include description
 --insert=0 shuffle only, do not insert unshuffled
 --n=1      number of shuffles
 --window   size of shuffle window

=head1 DESCRIPTION

shuffle_embed.pl takes a fasta formatted protein or DNA sequence file,
reads the sequence, shuffles it, splits the shuffled sequence in the
middle, and embeds the unshuffled sequence between the two halves of
the shuffled sequence.

With --insert 0, the sequences produced are random, no unshuffled
sequence is embedded.

=cut
