#!/usr/bin/perl -w

# convert FASTA .res file to fields for 'R'

# lose the first line
<>;

print "len\tscore\n";

my @line;

while(<>) {
    last if (m/\/\*\*/);

    @line = split(/\s+/);

#fields are:
# [0] ACC; [1] 0; [2] len; [3] frame; [4] comp; [5] H; [6-8] score[0-2];
# [9] rst.escore [10] segnum; [11] seglen; [12]lseek

    print "$line[2]\t$line[6]\n";
}
