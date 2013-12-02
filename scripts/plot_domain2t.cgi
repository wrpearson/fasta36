#!/usr/bin/perl -w

# plot_domain2.pl - produce SVG plot for aligned domains
# version2 plots both n0 and n1 sequences, with 2 axes
#
# args:
#  q_cstop - n0 - query length
#  l_cstop - n1 - lib length
#  q_name - query_acc
#  l_name - library_acc
#  l_astart= lib start  (need q_start, q_astop, l_astart, l_astop)
#  l_astop= lib stop 
#  pgm = program used 
#  regions -- same as annotations on alignment
#  doms -- domains on (library) sequence
#
#  l_annot - script to run to annotate library domain (separate from alignment)
#  q_annot - script to run to annotate library domain (separate from alignment)
#

# 9-May-2013 -- modify to accomodate reverse-complement coordinates

use strict;
use Getopt::Long;
use Pod::Usage;

use CGI qw(header param end_html);
#use URI::Escape;

use vars qw($pminx $pmaxx $pminy $pmaxy $lvstr $max_x $max_y
	    $fxscal $fyscal $fxoff $fyoff $x_rev $y_rev
	    @block_colors
	    $annot_color %annot_names %color_names);

@block_colors = qw( slategrey lightgreen lightblue pink cyan tan gold plum mediumplum );

$annot_color = 1;
%annot_names = ();

my @annot_scripts = ("", "",
		     "ann_feats2ipr.pl",
		     "ann_feats2l.pl",
		     "ann_feats2ipr.pl",
		     "ann_pfam26.pl",
		     "ann_pfam26.pl --pfacc",
		     "ann_pdb_cath.pl",
		     "ann_pdb_cath.pl --class",
		    );

# $max_x, $max_y define the maximum plotting area
# the actual bounding box/view area will be larger if annotation comments are available
($max_x,$max_y)=(540,24);

my @xax_len = (200, 300, 500, 1000, 2000, 3000, 5000, 10000, 20000, 30000, 50000);
my $max_xax = -1;
my $comb_xax = 0;  # comb_xax captures the length of the two sequences minus the offset
my ($x0c_off, $x1c_off, $xdc_off) = (0,0,0);

my ($x0f3, $x1f3) = (1,1);	# set to 3 for fastx (x0f3) or tfastx (x1f3)

# tick array - values */
my @tarr = (50, 100,200,500,1000,2000,5000,10000,20000);
my $MAX_INTERVAL=20000;

my $q = new CGI;

my %dom_colors = ();
my $max_color = 0;

my @arg_names = $q->param();

my %args = map { $_ => $q->param($_) } @arg_names;

if ($args{pgm} =~ m/^f[xy]$/) { $x0f3 = 3;}
elsif ($args{pgm} =~ m/^tf[xy]/) { $x1f3 = 3;}

#unless ($ENV{DOCUMENT_ROOT}) {
#  %args = map { $_ => uri_unescape($args{$_}) } keys %args;
#}

my ($region_info_r, $q_dom_info_r, $l_dom_info_r);

if ($args{regions}) {
  $region_info_r = parse_regions($args{regions});
}
else {$region_info_r = [];}

if ($args{doms}) {
  ($q_dom_info_r, $l_dom_info_r) = parse_domains($args{doms});
} else {
  $q_dom_info_r = [];
  $l_dom_info_r = [];
}

my @q_annots = ();
# unless (scalar(@{$q_dom_info_r})) {
#   my $q_annot_script = "";
#   if (defined($args{q_annot}) && $args{q_annot}) {
#     $q_annot_script = $annot_scripts[$args{q_annot}];
#   }
#   if (!$q_annot_script && defined($args{l_annot}) && $args{l_annot}) {
#     $q_annot_script = $annot_scripts[$args{l_annot}];
#   }

#   if ($q_annot_script) {
#     open(S_IN, '-|',"./$q_annot_script --lav \'sp|$args{q_name}'") || die "cannot open $args{q_name}\n";
#     while (my $s_line = <S_IN>) {
#       next if ($s_line =~ m/^#/);
#       next if ($s_line =~ m/^>/);
#       chomp($s_line);
#       my %q_data = ();
#       @q_data{qw(beg end descr)} = split(/\t/,$s_line);
#       if ($dom_colors{$q_data{descr}}) {
#   	$q_data{color} = $dom_colors{$q_data{descr}}
#       } else {
#   	$q_data{color} = ++$max_color;
#   	$dom_colors{$q_data{descr}} = $max_color;
#       }
#       push @q_annots, \%q_data;
#     }
#   }
# }


openplt(($args{q_cstop}-$args{q_cstart})+1, ($args{l_cstop}-$args{l_cstart})+1, $q_dom_info_r, $l_dom_info_r);
draw_align(\%args);
if (scalar(@{$region_info_r})) {
  draw_regions($region_info_r, $args{l_cstop});
}

my $q_annot_script = "";
if ($args{doms}) {
  if (scalar(@{$l_dom_info_r})) {
    draw_doms($l_dom_info_r, $x1c_off, -12, $args{l_cstart}, $args{l_cstop});
  }
  if (scalar(@{$q_dom_info_r})) {
    draw_doms($q_dom_info_r, $x0c_off, 48, $args{q_cstart}, $args{q_cstop});
  }
  elsif (scalar(@q_annots)) {
    draw_doms(\@q_annots, $x0c_off, 48, $args{q_cstart}, $args{q_cstop});
  }
}

closeplt($args{l_cstop});

exit(0);

# have all the data (and length of sequence), scale it and color it

#define SX(x) (int)((double)(x)*fxscal+fxoff+6)
sub SX {
  my $xx = shift;
  return int($xx*$fxscal+$fxoff+18);
}

sub SY {
  my $yy = shift;
  return $max_y - int($yy*$fyscal+$fyoff);
}

my $y_delta = 0;

#void openplt(long n0, long n1, int sq0off, int sq1off, char *xtitle, char *ytitle)
sub openplt
{
  my ($n0, $n1, $q_dom_info_r, $l_dom_info_r) = @_;

  my ($xbound, $ybound) = ($max_x + 24, 48);
  my ($x0_rev, $x1_rev) = (0,0);

  if (scalar(@{$q_dom_info_r})) {
    $ybound += 14;
  }
  elsif (scalar(@q_annots)) {
    $ybound += 14;
  }
  $ybound += 14 if (scalar(@{$l_dom_info_r}));

  if ($n0 < 0) {$x0_rev=1; $n0 = 2 - $n0;}
  if ($n1 < 0) {$x1_rev=1; $n1 = 2 - $n1;}

  print $q->header('image/svg+xml') if ($ENV{DOCUMENT_ROOT});
  print("<?xml version=\"1.0\" standalone=\"no\"?>\n");
  print("<!DOCTYPE svg PUBLIC \"-//W3C//DTD SVG 1.1//EN\" \n");
  print("\"http://www.w3.org/Graphics/SVG/1.1/DTD/svg11.dtd\">\n\n");

  print(qq(<!-- l_name=$args{l_name} -->\n));
  print("<svg width=\"$xbound\" height=\"$ybound\" version=\"1.1\"\n");
#  print("<svg version=\"1.1\"\n");
  print("xmlns=\"http://www.w3.org/2000/svg\" xmlns:xlink=\"http://www.w3.org/1999/xlink\">\n\n");

# simple things first, if query is shorter and inside, or library is
# shorter and inside, then just use the longer sequence.

  ($x0c_off, $x1c_off, $xdc_off) = (0,0,0);

  ($comb_xax, $xdc_off, $x0c_off, $x1c_off) = calc_offsets($n0, $x0f3, $x0_rev, $n1, $x1f3, $x1_rev, \%args);

#  if ($args{pgm} =~ m/lsw/) {
#    $max_xax = $comb_xax;
#  } else {
    for (my $i=0; $i < scalar(@xax_len); $i++) {
      if ($comb_xax <= $xax_len[$i]) {
	$max_xax = $xax_len[$i];
	last;
      }
#    }
    $max_xax = $xax_len[$#xax_len] if ($max_xax <= 0);
  }

  $fxscal = ($max_x-1)/$max_xax;
  $fyscal = 1;

  $fxscal *= 0.9; $fxoff = 24;
  $fyoff = -14;
  if (scalar(@{$q_dom_info_r}) || scalar(@q_annots)) {$fyoff -= 12;}


  xaxis2($n0/$x0f3,$args{q_cstart}/$x0f3, $x0_rev, $n1/$x1f3, $args{l_cstart}/$x1f3, $x1_rev);

  newline(qq(stroke="black" stroke-width="1.5"));
  move(SX($x0c_off), SY(15));
  draw(SX($x0c_off+($n0/$x0f3)),SY(15));
  clsline($n0,0);

  newline(qq(stroke="black" stroke-width="1.5"));
  move(SX($x1c_off), SY(9));
  draw(SX($x1c_off+($n1/$x1f3)),SY(9));

  clsline($n1,0);
}

sub closeplt {
  print "</svg>\n";
}

sub calc_offsets {
  my ($n0, $x0f3, $x0_rev, $n1, $x1f3, $x1_rev, $args_r) = @_;

  my ($comb_xax, $x0c_off, $x1c_off, $xdc_off) = (0,0,0,0);

  my ($n0f, $n1f, $q_start_d, $q_stop_d, $l_start_d, $l_stop_d) =
    (
     $n0/$x0f3, # $n0f
     $n1/$x0f3, # $n1f
     abs($args_r->{q_astart} - $args_r->{q_cstart})/$x0f3, # $q_start_d
     abs($args_r->{q_astop} - $args_r->{q_cstop})/$x0f3,   # $q_stop_d
     abs($args_r->{l_astart} - $args_r->{l_cstart})/$x1f3, # $l_start_d
     abs($args_r->{l_astop} - $args_r->{l_cstop})/$x1f3    # $l_stop_d
    );

  if (($n1f >= $n0f) && ($l_start_d >= $q_start_d) && ($l_stop_d >= $q_stop_d)) {
    # n1 is longer and n0 is contained
    $comb_xax = $n1f; # $comb_xax : combined x-axis, in amino-acids if translated
    $x0c_off = $l_start_d - $q_start_d;
  }
  elsif (($n1f < $n0f) && ($l_start_d <= $q_start_d) && ($l_stop_d <= $q_stop_d)) {
    # n0 is longer and n1 is contained
    $comb_xax = $n0f;
    $xdc_off = $x1c_off = $q_start_d - $l_start_d;
  }
    # some kind of extension is necessary
  elsif ($l_start_d >= $q_start_d) {
    $x0c_off = $l_start_d - $q_start_d;
    $comb_xax = $n0f + $x0c_off;
  }
  else {
    $xdc_off = $x1c_off = $q_start_d - $l_start_d;
    $comb_xax = $n1f + $x1c_off;
  }

  return ($comb_xax, $xdc_off, $x0c_off, $x1c_off);
}

sub draw_trapz {
  my ($start0, $stop0, $start1, $stop1, $color, $text) = @_;

  $color = ($color % scalar(@block_colors));
  my $svg_color = $block_colors[$color];
  my $tx = $start1 + int(($stop1-$start1+1)/2);
  my $ty = 10 + 9;

  $text = substr($text,0,10);

  newline(qq(stroke="black" stroke-width="1.5"));

  move(SX($start0+$x0c_off), SY(20));
  draw(SX($stop0+$x0c_off),SY(20));
  draw(SX($stop1+$x1c_off),SY(10));
  draw(SX($start0+$x1c_off),SY(10));
  move(SX($start0+$x0c_off), SY(20));

  print(qq(" fill="$svg_color" />\n));


#  print (qq(<rect x="$x" y="$y" width="$w" height="$h" fill="$svg_color" stroke="white" stroke-width="1" />\n));
#  print (qq(<text x="$tx" y="$ty" font-size="9" font-family="sans-serif" fill="white" text-anchor="middle">$text</text>\n));
}

# draws a colored solid block, and labels it, to indicate domain
sub draw_block {
  my ($x, $y, $w, $h, $color, $text, $Q) = @_;

  $color = ($color % scalar(@block_colors));
  my $svg_color = $block_colors[$color];
  my $tx = $x + int($w/2);
  my $ty = $y + 9;

  $text = substr($text,0,10);

  my $stroke_width = 0.5;
  if ($Q < 30.0) {$stroke_width = 2;}

  print (qq(<rect x="$x" y="$y" width="$w" height="$h" fill="$svg_color" stroke="white" stroke-width="$stroke_width" />\n));
  print (qq(<text x="$tx" y="$ty" font-size="9" font-family="sans-serif" fill="white" text-anchor="middle">$text</text>\n));
}

sub draw_regions {
  my ($annot_arr_r, $n1) = @_;

  for my $annot ( @$annot_arr_r) {
    draw_block(SX($annot->{beg1} - $args{l_cstart}+$xdc_off), SY(18), SX($annot->{end1}+$xdc_off)-SX($annot->{beg1}+$xdc_off),
	       12, $annot->{color}, $annot->{descr}, $annot->{Q});
  }
}

sub draw_doms {
  my ($annot_arr_r, $xc_off, $y_off, $xc_start, $xc_stop) = @_;

  for my $annot ( @$annot_arr_r) {
    draw_block(SX($annot->{beg}+$xc_off), SY($y_off), SX($annot->{end}+$xc_off)-SX($annot->{beg}+$xc_off),
	       12, $annot->{color}, $annot->{descr}, 100.0);
  }
}

sub draw_align {
  my $arg_r = shift;

  my ($x, $y, $w, $h) = (SX($args{l_astart} - $args{l_cstart} + $xdc_off), SY(21), SX($args{l_astop}+$xdc_off) - SX($args{l_astart}+$xdc_off), 18);

  print (qq(<rect x="$x" y="$y" width="$w" height="$h" stroke="black" fill-opacity="0" stroke-width="1" />\n));
}

# void newline(char *options)
sub newline
{
  my $options = shift;

  if ($options)  {
      printf("<path %s d=\"",$options);
  }
  else {print("<path stroke=\"black\" d=\"");}
}

# void clsline(long x, long y, int s)
sub clsline
{
  my ($x, $y, $s) = @_;

  print("\" fill=\"none\" />\n");
}

#void move(int x, int y)
sub move
{
  my ($x, $y) = @_;
  printf(" M %d %d",$x,$y);
}

# void draw(int x, int y)
sub draw
{
  my ($x, $y) = @_;
  printf(" L %d %d",$x,$y);
}

# void xaxis(long n, int offset, char *title)
# coordinates in amino acids - modify for final axes
sub xaxis2 {
  my ($n0, $offset0, $x0_rev, $n1, $offset1,$x1_rev) = @_;

  my ($v_offset0, $v_offset1) = ($offset0, $offset1);

  my ($i, $jm, $tick_length, $max_ticks, $tick_inc);
  my ($js, $jl0, $jl1);
  my ($sgn0, $sgn1) = (1,1);

  my $num_len;
  my $numstr;

  if ($x0_rev) {
    $sgn0 = -1;
    $v_offset0 = $offset0 - $n0;
  }

  if ($x1_rev) {
    $sgn1 = -1;
    $v_offset1 = $offset1 - $n1;
  }

  # for translated-DNA/protein searches, both $n0 and $n1 are in amino-acids
  my $n_max = $n1;
  $n_max = $n0 if ($n0 > $n1);
  my $offset = 0;

  $tick_length = 2;

  # search for the correct increment for the tick array */
  # @tarr[] has the list of axis increments we might use
  # we want a tick increment that gives < 6 ticks
  for ($i=0; $i< @tarr; $i++) {
    # seek to divide into 10 or fewer divisions */
    if (($max_ticks = $n_max/$tarr[$i]) <= 6) {goto found;}
  }

  # these happen only if no tick increment was found
  # point $i to the last element
  $i = scalar(@tarr)-1;

  # $max_ticks is the number of increments for longest sequence

  $max_ticks = ($n_max)/$tarr[-1];
  $i = -1;

 found:
  $max_ticks += max($offset0, $offset1)/$tarr[$i];
  $tick_inc = $tarr[$i];

  # jo is the offset for the coordinate system, e.g. if we are
  # plotting an alignment from 101 - 300 rather than 1 - 400 we may
  # show partial sequences in alignments, it should be kept, but is is
  # different from the axis shift ($x0c_off, $x1c_off)

  my ($xx0c_off, $xx1c_off) = ($x0c_off, $x1c_off);

  unless ($x0_rev) {
#    $xx0c_off -= $offset0;
# draw up-tick
    newline("stroke=\"black\" stroke-width=\"1.5\"");
    for ($i=1; $i<=$max_ticks; $i++) {
      last if ($i*$tick_inc > $n0 + $v_offset0);
      next if ($i*$tick_inc < $v_offset0);

      move(SX($i*$tick_inc + $xx0c_off - $v_offset0),SY(26));
      draw(SX($i*$tick_inc + $xx0c_off - $v_offset0),SY(26)+$tick_length);
    }
    clsline($n_max,$n_max,10000);

    for ($i=1; $i<=$max_ticks; $i++) {
      last if ($i*$tick_inc > $n0 + $v_offset0);
      next if ($i*$tick_inc < $v_offset0);
      $numstr = sprintf("%ld",$i*$tick_inc*$x0f3 );
      printf(qq(<text x="%d" y="%d" font-family="sans-serif" font-size='9' text-anchor="middle">%s</text>\n),
	     SX($i*$tick_inc+$xx0c_off - $v_offset0),SY(28)+$tick_length-1,$numstr);
    }
  }
  else {
    $xx0c_off += $offset0;
    # if $x0_rev need to know $x0_max_ticks
    my $x0_max_ticks = int(($n0+$offset0)/$tick_inc);

    newline("stroke=\"black\" stroke-width=\"1.5\"");
    for ($i=$x0_max_ticks; $i>0; $i--) {
      move(SX($xx0c_off - $i*$tick_inc),SY(26));
      draw(SX($xx0c_off - $i*$tick_inc),SY(26)+$tick_length);
    }
    clsline($n_max,$n_max,10000);

    # now put in the numbers, using the same reverse counting
    for ($i=$x0_max_ticks; $i>0; $i--) {
      $numstr = sprintf("%ld",$i*$tick_inc*$x0f3);
      printf(qq(<text x="%d" y="%d" font-family="sans-serif" font-size='9' text-anchor="middle">%s</text>\n),
	     SX($xx0c_off - $i*$tick_inc - $v_offset0),SY(28)+$tick_length-1,$numstr);
    }
  }

  unless ($x1_rev) { 
#    $xx1c_off -= $offset1;
    # draw down-tick
    newline("stroke=\"black\" stroke-width=\"1.5\"");
    for ($i=1; $i<=$max_ticks; $i++) {
      last if ($i*$tick_inc > $n1 + $v_offset1);
      next if ($i*$tick_inc < $v_offset1);
      move(SX($i*$tick_inc + $xx1c_off - $v_offset1),SY(0));
      draw(SX($i*$tick_inc + $xx1c_off - $v_offset1),SY(0)+$tick_length);
    }
    clsline($n_max,$n_max,10000);

    $numstr = sprintf("%ld",$tick_inc*$x0f3);
    $num_len = length($numstr);

    for ($i=1; $i<=$max_ticks; $i++) {
      last if ($i*$tick_inc > $n1 + $offset1);
      next if ($i*$tick_inc < $offset1);
      $numstr = sprintf("%ld",$i*$tick_inc*$x1f3);
      printf(qq(<text x="%d" y="%d" font-family="sans-serif" font-size='9' text-anchor="middle">%s</text>\n),
	     SX($i*$tick_inc+$xx1c_off-$v_offset1),SY(0)+$tick_length+8,$numstr);
    }
  }
  else {
    my $x1_max_ticks = int($n1/$tick_inc);
    $xx1c_off += $offset1;

    # draw down-tick
    newline("stroke=\"black\" stroke-width=\"1.5\"");
    for ($i=$x1_max_ticks; $i>0; $i--) {
      move(SX($xx1c_off - $i*$tick_inc),SY(0));
      draw(SX($xx1c_off - $i*$tick_inc),SY(0)+$tick_length);
    }
    clsline($n_max,$n_max,10000);

    $numstr = sprintf("%ld",$tick_inc*$x0f3);
    $num_len = length($numstr);

    for ($i=$x1_max_ticks; $i>0; $i--) {
      $numstr = sprintf("%ld",$i*$tick_inc*$x1f3);
      printf(qq(<text x="%d" y="%d" font-family="sans-serif" font-size='9' text-anchor="middle">%s</text>\n),
	     SX($xx1c_off - $i*$tick_inc),SY(0)+$tick_length+8,$numstr);
    }
  }
}

# void xaxis(long n, int offset, char *title)
sub xaxis_d {
  my ($n, $offset) = @_;

  my ($i, $jm, $tick);
  my ($js, $jo, $jl);
  my $num_len;
  my $numstr;

  $tick = 2;

  # search for the correct increment for the tick array */
  for ($i=0; $i< @tarr; $i++) {
    # seek to divide into 10 or fewer divisions */
    if (($jm = $n/$tarr[$i])<6) {goto found;}
  }
  $i= scalar(@tarr)-1;
  $jm = $n/$tarr[$i];
 found:
  # js is the start of the value - modify to accomodate offset */
  $js = $tarr[$i];

  # jo is the offset */
  $jo = $offset % $tarr[$i];	# figure out offset in tarr[i] increments */

  # jl is the label */
  $jl = $offset/$tarr[$i];	# figure out offset in tarr[i] increments */
  $jl *= $tarr[$i];

  newline("stroke=\"black\" stroke-width=\"1.5\"");
  for ($i=1; $i<=$jm; $i++) {
    move(SX($i*$js - $jo),SY(0));
    draw(SX($i*$js - $jo),SY(0)+$tick);
  }
  clsline($n,$n,10000);

  $numstr = sprintf("%ld",$js + $jl );
  $num_len = length($numstr);
  if ($num_len > 4) {

    printf(qq(<text x="%d" y="%d" font-family="sans-serif" font-size='9' text-anchor="middle">%s</text>\n),SX($js-$jo),SY(0)+$tick+8,$numstr);

    $numstr = sprintf("%ld",$jm*$js+$jl);
    printf(qq(<text x="%d" y="%d" font-family="sans-serif" font-size='9' text-anchor="middle">%s</text>\n),SX($jm*$js-$jo),SY(0)+$tick+8,$numstr);
  }
  else {
    for ($i=1; $i<=$jm; $i++) {
      $numstr = sprintf("%ld",$i*$js+$jl);
      printf(qq(<text x="%d" y="%d" font-family="sans-serif" font-size='9' text-anchor="middle">%s</text>\n),SX($i*$js-$jo),SY(0)+$tick+8,$numstr);
    }
  }
}

sub parse_regions {
  my $region_str = shift;

  my @regions = split(/\n\s*/,$region_str);

  my @region_info = ();

  for my $region ( @regions) {
    $region =~ s/^\s+//;
    next unless ($region =~ m/^Region/);

    my @fields = split(/\s+:\s*/,$region);

    my %data = ();

    @data{qw(descr color)} = @fields[-2,-1];

    $dom_colors{$data{descr}}=$data{color} unless defined($dom_colors{$data{descr}});
    $max_color = $data{color} if ($data{color} > $max_color);

    my @scores = split(/;\s*/,$fields[1]);

    for my $score (@scores) {
      my ($key, $value) = split(/=/,$score);
      $data{$key} = $value;
    }

    # this line hides low-score NODOMs
    next if ($data{color}==0 && $data{Q} < 30.0);

    @data{qw(beg0 end0 beg1 end1)} = ($fields[0] =~ m/^Region:\s*(\d+)-(\d+):(\d+)-(\d+)$/);

    push @region_info, \%data;
  }

  return \@region_info;
}

sub parse_domains {
  my $domain_str = shift;

  my @domains = split(/\n\s*/,$domain_str);

  my @q_domain_info = ();
  my @l_domain_info = ();

  for my $domain ( @domains) {
    $domain =~ s/^\s+//;
    next unless ($domain =~ m/^[ql]Domain/);

    my @fields = split(/\t/,$domain);

    next if ($fields[-1] =~ m/NODOM/);

    my %data = ();

    @data{qw(beg end)}  = ($fields[1]) =~ m/(\-?\d+)\-(\-?\d+)/;
    @data{qw(descr color)} = split(/ :/,$fields[-1]);

    $dom_colors{$data{descr}}=$data{color} unless defined($dom_colors{$data{descr}});
    $max_color = $data{color} if ($data{color} > $max_color);

    if ($fields[0] =~ m/^qDomain/) {
	$data{beg} -= $args{q_cstart} + 1;
	$data{end} -= $args{q_cstart} + 1;
	next if $data{end} < 1;
	$data{beg} = 1 if $data{beg} < 1;
	next if $data{beg} > $args{q_cstop};
	$data{end} = $args{q_cstop} if $data{end} > $args{q_cstop};

	push @q_domain_info, \%data;
    }
    elsif ($fields[0] =~ m/^lDomain/) {
	$data{beg} -= $args{l_cstart} + 1;
	$data{end} -= $args{l_cstart} + 1;
	next if $data{end} < 1;
	$data{beg} = 1 if $data{beg} < 1;
	next if $data{beg} > $args{l_cstop};
	$data{end} = $args{l_cstop} if $data{end} > $args{l_cstop};
	push @l_domain_info, \%data;
    }
  }

  return (\@q_domain_info, \@l_domain_info);
}

sub max {
  my ($x0, $x1) = @_;

  return $x0 if $x0 >= $x1;
  return $x1;
}
