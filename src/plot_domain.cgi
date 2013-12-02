#!/usr/bin/perl -w

# domain_plt.pl - produce SVG plot for aligned domains

use strict;
use Getopt::Long;
use Pod::Usage;

use CGI qw(header param end_html);
#use URI::Escape;

use vars qw($pminx $pmaxx $pminy $pmaxy $lvstr $max_x $max_y
	    $fxscal $fyscal $fxoff $fyoff
	    @block_colors
	    $annot_color %annot_names %color_names);

@block_colors = qw( slategrey lightgreen lightblue pink cyan tan gold plum mediumplum );

$annot_color = 1;
%annot_names = ();

# $max_x, $max_y define the maximum plotting area
# the actual bounding box/view area will be larger if annotation comments are available
($max_x,$max_y)=(540,24);

my @xax_len = (200, 300, 500, 1000, 2000, 3000, 5000, 10000, 20000, 30000, 50000);
my $max_xax = -1;
# tick array - values */
my @tarr = (50, 100,200,500,1000,2000,5000,10000,20000);
my $MAX_INTERVAL=20000;

my $q = new CGI;

my @arg_names = $q->param();

my %args = map { $_ => $q->param($_) } @arg_names;

#unless ($ENV{DOCUMENT_ROOT}) {
#  %args = map { $_ => uri_unescape($args{$_}) } keys %args;
#}

my @regions = split(/\n\s*/,$args{doms});

my @region_info = ();

for my $region ( @regions) {
  $region =~ s/^\s+//;
  next unless ($region =~ m/^Region/);

  my @fields = split(/\s+:\s*/,$region);

  my %data = ();

  @data{qw(descr color)} = @fields[-2,-1];

  my @scores = split(/;\s*/,$fields[1]);

  for my $score (@scores) {
    my ($key, $value) = split(/=/,$score);
    $data{$key} = $value;
  }

  next if ($data{color}==0 && $data{Q} < 30.0);

  @data{qw(beg0 end0 beg1 end1)} = ($fields[0] =~ m/^Region:\s*(\d+)-(\d+):(\d+)-(\d+)$/);

  push @region_info, \%data;
}

if (scalar(@region_info)) {
  openplt($args{n1});
  draw_align(\%args);
  draw_doms(\@region_info, $args{n1});
  closeplt($args{n1});
}

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
  my ($n1) = @_;

  my ($xbound, $ybound) = ($max_x + 24, 40);

  print $q->header() if ($ENV{DOCUMENT_ROOT});
  print("<?xml version=\"1.0\" standalone=\"no\"?>\n");
  print("<!DOCTYPE svg PUBLIC \"-//W3C//DTD SVG 1.1//EN\" \n");
  print("\"http://www.w3.org/Graphics/SVG/1.1/DTD/svg11.dtd\">\n\n");

  print("<svg width=\"$xbound\" height=\"$ybound\" version=\"1.1\"\n");
#  print("<svg version=\"1.1\"\n");
  print("xmlns=\"http://www.w3.org/2000/svg\" xmlns:xlink=\"http://www.w3.org/1999/xlink\">\n\n");

  if ($args{o_pgm} =~ m/lsw/) {
    $max_xax = $n1;
  } else {
    for (my $i=0; $i < scalar(@xax_len); $i++) {
      if ($n1 <= $xax_len[$i]) {
	$max_xax = $xax_len[$i];
	last;
      }
    }
    $max_xax = $xax_len[$#xax_len] if ($max_xax <= 0);
  }


  $fxscal = ($max_x-1)/$max_xax;
  $fyscal = 1;

  $fxscal *= 0.9; $fxoff = 24;
  $fyoff = 0;

  xaxis($n1,0);

  newline("stroke=\"black\" stroke-width=\"1.5\"");
  move(SX(0), SY(12));
  draw(SX($n1),SY(12));

  clsline($n1,0);
}

sub closeplt {
  print "</svg>\n";
}

sub draw_block {
  my ($x, $y, $w, $h, $color, $text) = @_;

  $color = ($color % scalar(@block_colors));
  my $svg_color = $block_colors[$color];
  my $tx = $x + int($w/2);
  my $ty = $y + 9;

  $text = substr($text,0,10);

  print (qq(<rect x="$x" y="$y" width="$w" height="$h" fill="$svg_color" stroke="white" stroke-width="1" />\n));
  print (qq(<text x="$tx" y="$ty" font-size="9" font-family="sans-serif" fill="white" text-anchor="middle">$text</text>\n));
}

sub draw_doms {
  my ($annot_arr_r, $n1) = @_;

  for my $annot ( @$annot_arr_r) {
    draw_block(SX($annot->{beg1}), SY(18), SX($annot->{end1})-SX($annot->{beg1}),
	       12, $annot->{color}, $annot->{descr});
  }
}

sub draw_align {
  my $arg_r = shift;

  my ($x, $y, $w, $h) = (SX($arg_r->{start}), SY(21), SX($arg_r->{stop}) - SX($arg_r->{start}), 18);

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

# void sxy_move(int x, int y)
sub sxy_move
{
  my ($x, $y) = @_;
  move(SX($x), SY($y));
}

# void draw(int x, int y)
sub draw
{
  my ($x, $y) = @_;
  printf(" L %d %d",$x,$y);
}

# void sxy_draw(int x, int y)
sub sxy_draw
{
  my ($x, $y) = @_;
  draw(SX($x),SY($y));
}

# void xaxis(long n, int offset, char *title)
sub xaxis {
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
