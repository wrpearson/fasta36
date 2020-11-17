#!/usr/bin/env perl
#
################################################################
# copyright (c) 2012, 2014 by William R. Pearson and The Rector &
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

use warnings;
use strict;

use vars qw($have_zdb $have_bits);

#define SX(x) (int)((double)(x)*fxscal+fxoff+6)
sub SX {
  my $xx = shift;
  return int($xx*$fxscal+$fxoff+18);
}

#define SY(y) (max_y + 24 - (int)((double)(y)*fyscal+fyoff))
sub SY {
  my $yy = shift;
  return $max_y + 24 - int($yy*$fyscal+$fyoff);
}

# $y_delta used widely to offset for x_domain annotation
my $y_delta = 0;

#void openplt(long n0, long n1, int sq0off, int sq1off, char *xtitle, char *ytitle)
sub openplt
{
  my ($n0, $n1, $sq0off, $sq1off, $xtitle, $ytitle, $x_annot_r, $y_annot_r,$have_zdb, $have_bits) = @_;

  if ($lvstr) {
    @elinval = split(/\s+/,$lvstr);
  }
  elsif ($ENV{LINEVAL}) {
    @elinval = split(/\s+/,$ENV{LINEVAL});
  }

  my ($xbound, $ybound) = ($max_x + 24, $max_y + 48);
  $y_delta = 0;

  if ($x_annot_r) {
    my $x_comments = 0;
    for my $annot (@$x_annot_r) {
      if (length($annot->{sdescr}) > $x_comments) {
	$x_comments = length($annot->{sdescr})
      }
    }
    $ybound += 24 + 6 * $x_comments;
    $y_delta += 24;
  }

  if ($y_annot_r) {
    my $y_comments = 0;
    for my $annot (@$y_annot_r) {
      if (length($annot->{sdescr}) > $y_comments) {
	$y_comments = length($annot->{sdescr})
      }
    }
    $xbound += 6 * $y_comments;
    $max_x -= 6 * $y_comments;
    $max_y -= 6 * $y_comments
  }

  print("<?xml version=\"1.0\" standalone=\"no\"?>\n");
  print("<!DOCTYPE svg PUBLIC \"-//W3C//DTD SVG 1.1//EN\" \n");
  print("\"http://www.w3.org/Graphics/SVG/1.1/DTD/svg11.dtd\">\n\n");

  print("<svg width=\"$xbound\" height=\"$ybound\" version=\"1.1\"\n");
  print("xmlns=\"http://www.w3.org/2000/svg\" xmlns:xlink=\"http://www.w3.org/1999/xlink\">\n\n");

  ($pmaxx, $pmaxy) = ($n0, $n1);

  $fxscal = ($max_x-1)/$n1;
  $fyscal = ($max_y-1)/$n0;

  if ($fxscal > $fyscal) {$fxscal = $fyscal;}
  else {$fyscal = $fxscal;}

  if ($fyscal * $n0 < $max_y/5.0) {
    $fyscal = ($max_y-1)/($n0*5.0);
  }

  $fxscal *= 0.9; $fxoff = ($max_x-1)/11.0;
  $fyscal *= 0.9; $fyoff = ($max_y-1)/11.0;
  if ($x_annot_r) {$fyoff -= (48 + $y_delta);}

  # draw the plot frame
  printf("<rect x=\"%d\" y=\"%d\" width=\"%d\" height=\"%d\"\n",
	 SX(0),SY($n1+1), SX($n0+1)-SX(0), SY(0) - SY($n1+1));
  print("stroke=\"black\" stroke-width=\"2.0\" fill=\"none\" />\n");

  my $n_div = 11;
  xaxis($n0, $sq1off, $xtitle, $n_div);

  $n_div = 21 unless ($n0 == $n1);
  yaxis($n1, $sq0off, $ytitle, $n_div);
  legend($have_zdb, $have_bits, ($x_annot_r));

  if ($x_annot_r) {xgrid($x_annot_r, $n0, $sq0off, $n1, $sq1off);}
  if ($y_annot_r) {ygrid($y_annot_r, $n0, $sq0off, $n1, $sq1off);}
}
	
# void drawdiag(long n0, long n1)
sub drawdiag
{
  my ($n0, $n1) = @_;
  # 	printf("currentlinewidth 1.5 mul setlinewidth\n"); */
  newline("stroke=\"black\" stroke-width=\"1.5\"");
  move(SX(0),SY(0));
  draw(SX($n0+1),SY($n1+1));
  clsline($n0,$n1,10000);
}

# tick array - values */
my @tarr = (10,20,50,100,200,500,1000,2000,5000,10000,20000,50000,100000,200000,500000,1000000);
my $MAX_INTERVAL=1000000;

# void xaxis(long n, int offset, char *title)
sub xaxis {
  my ($n, $offset, $title, $n_div) = @_;

  my ($i, $jm, $tick);
  my ($js, $jo, $jl);
  my $num_len;
  my $numstr;

  $tick = 6;

  # search for the correct increment for the tick array */
  for ($i=0; $i< @tarr; $i++) {
    # seek to divide into 20 or fewer divisions */
    if (($jm = $n/$tarr[$i])< $n_div) {goto found;}
  }
  $i= scalar(@tarr)-1;
  $jm = $n/$tarr[$i];
 found:
  # js is the start of the value - modify to accomodate offset */
  $js = $tarr[$i];

  # jo is the offset */
  $jo = ($offset-1) % $tarr[$i];	# figure out offset in tarr[i] increments */

  # jl is the label */
  $jl = ($offset-1)/$tarr[$i];	# figure out offset in tarr[i] increments */
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

    printf("<text x=\"%d\" y=\"%d\" text-anchor=\"middle\">%s</text>\n",SX($js-$jo),SY(0)+$tick+16,$numstr);

    $numstr = sprintf("%ld",$jm*$js+$jl);
    printf("<text x=\"%d\" y=\"%d\" text-anchor=\"middle\">%s</text>\n",SX($jm*$js-$jo),SY(0)+$tick+16,$numstr);
  }
  else {
    for ($i=1; $i<=$jm; $i++) {
      $numstr = sprintf("%ld",$i*$js+$jl);
      printf("<text x=\"%d\" y=\"%d\" text-anchor=\"middle\">%s</text>\n",SX($i*$js-$jo),SY(0)+$tick+16,$numstr);
    }
  }

  printf("<text x=\"%d\" y=\"%d\" text-anchor=\"middle\">%s</text>\n",SX($n/2),SY(0)-$tick+42, $title);
}

sub xgrid {
  my ($annot_arr_r, $n0, $sq0_off, $n1, $sq1_off) = @_;

  my $sq_off = $sq0_off;

  my $show_block = 1;
  my $text_offset = 8;

  if ($show_block) {$text_offset = 24;}

  for my $annot ( @$annot_arr_r) {
    next unless $annot->{beg} >= $sq_off;
    next if ($annot->{end} > $sq_off + $n0 - 1);
    last if ($annot->{beg} > $sq_off + $n0 - 1);
    newline("stroke=\"black\" stroke-width=\"1.5\" stroke-opacity=\"0.33\" stroke-dasharray=\"3,6\"" );
    move(SX($annot->{beg} - $sq_off),SY(0));
    draw(SX($annot->{beg} - $sq_off),SY($n1));
    clsline();
    newline("stroke=\"black\" stroke-width=\"1.5\" stroke-opacity=\"0.33\" stroke-dasharray=\"6,3\"" );
    move(SX($annot->{end} - $sq_off),SY(0));
    draw(SX($annot->{end} - $sq_off),SY($n1));
    clsline();

    if ($show_block) {
      draw_block(SX($annot->{beg} - $sq_off), SY($n1) - 18,
		 SX($annot->{end} - $sq_off) - SX($annot->{beg} - $sq_off),
		 12, $annot_names{$annot->{sname}});
    }

    # show rotated label
    my $xpos = SX(($annot->{end} - $annot->{beg})/2 + $annot->{beg} - $sq_off) + 4;
    my $ypos = SY($n1) - $text_offset;
    printf("<text x=\"0\" y=\"0\" text-anchor=\"left\" transform=\"translate($xpos, $ypos) rotate(-90,0,0)\">%s</text>\n",$annot->{sdescr});
  }
}

sub ygrid {
  my ($annot_arr_r, $n0, $sq0_off, $n1, $sq1_off) = @_;

  my $sq_off = $sq1_off;

  my $show_block = 1;
  my $text_offset = 8;
  if ($show_block) {$text_offset = 24;}

  for my $annot ( @$annot_arr_r) {
    next unless $annot->{beg} >= $sq_off;
    next if ($annot->{end} > $sq_off + $n1 - 1);
    last if ($annot->{beg} > $sq_off + $n1 - 1);
    newline("stroke=\"black\" stroke-width=\"1.5\" stroke-opacity=\"0.33\" stroke-dasharray=\"3,6\"" );
    move(SX(0), SY($annot->{beg} - $sq_off));
    draw(SX($n0), SY($annot->{beg} - $sq_off));
    clsline();
    newline("stroke=\"black\" stroke-width=\"1.5\" stroke-opacity=\"0.33\" stroke-dasharray=\"6,3\"" );
    move(SX(0), SY($annot->{end} - $sq_off));
    draw(SX($n0), SY($annot->{end} - $sq_off));
    clsline();

    my $xpos = SX($n0) + $text_offset;
    my $ypos = SY(($annot->{end} - $annot->{beg})/2 + $annot->{beg} - $sq_off) + 4;

    if ($show_block) {
      draw_block(SX($n0)+6, SY($annot->{end} - $sq_off), 12,
		 SY($annot->{beg} - $sq_off) - SY($annot->{end} - $sq_off),
		 $annot_names{$annot->{sname}});
    }
    printf("<text x=\"$xpos\" y=\"$ypos\" text-anchor=\"left\">%s</text>\n",$annot->{sdescr});
  }
}

#void yaxis(long n, int offset, char *title)
sub yaxis {
  my ($n, $offset, $title, $n_div) = @_;

  my ($i, $jm, $tick);
  my ($js, $jo, $jl);
  my $num_len;
  my $numstr;

  $tick = 6;

  for ($i=0; $i< @tarr; $i++) {
    if (($jm = $n/$tarr[$i])< $n_div) {goto found;}
  }
  $jm = $n/5000;
  $i= scalar(@tarr)-1;

 found:
  $js = $tarr[$i];

  # jo is the offset */
  $jo = ($offset-1) % $tarr[$i];	# figure out offset in tarr[i] increments */
  # jl is the label */
  $jl = ($offset-1) / $tarr[$i];	# figure out offset in tarr[i] increments */
  $jl *= $tarr[$i];

  newline("stroke=\"black\" stroke-width=\"1.5\"");
  for ($i=1; $i<=$jm; $i++) {
    move(SX(0),SY($i*$js-$jo));
    draw(SX(0)-$tick,SY($i*$js-$jo));
  }
  clsline($n,$n,10000);

  $numstr = sprintf("%d",$js+$jl);
  $num_len = length($numstr);
  if ($num_len > 4) {
    printf("<text x=\"%d\" y=\"%d\" text-anchor=\"end\">%s</text>\n",SX(0)-$tick-4,SY($js-$jo)+4,$numstr);

    $numstr = sprintf("%ld",$jm*$js+$jl);
    printf("<text x=\"%d\" y=\"%d\" text-anchor=\"end\">%s</text>\n",SX(0)-$tick-4,SY($jm*$js-$jo)+4,$numstr);
  }
  else {
    for ($i=1; $i<=$jm; $i++) {
      $numstr = sprintf("%ld",$i*$js+$jl);
      printf("<text x=\"%d\" y=\"%d\" text-anchor=\"end\">%s</text>\n",SX(0)-$tick-4,SY($i*$js-$jo)+4,$numstr);
    }
  }
  # make a path for the label */

  #move(SX(0)-$tick-18,SY($n/2));
  printf(qq(<g transform="rotate(-90, 142, 142)">\n));
  printf(qq(<text x="0" y="12" text-anchor="middle">%s</text>\n),$title);
  print("</g>\n");
}

sub draw_block {
  my ($x, $y, $w, $h, $color) = @_;

  $color = ($color % scalar(@block_colors));
  my $svg_color = $block_colors[$color];

  print (qq(<rect x="$x" y="$y" width="$w" height="$h" fill="$svg_color" stroke="white" stroke-width="1" />));
}

sub legend
{
  my ($have_zdb, $have_bits, $annot_flg) = @_;

  my ($i, $last, $del);
  my ($ixp, $iyp);
  my $numstr;
  my $optstr;
  my @xpos=(144,144,288,288,432);
  my @ypos=(36,18,36,18,27);

  my $y_off = 66;
  if ($annot_flg) {$y_off = 120;}

  if ($have_zdb || $have_bits) {$last = 5;}
  else {$last = 4;}

  if ($have_zdb)  {printf("<text x=\"%d\" y=\"%d\">E(): </text>",54,$max_y + $y_off - 24 + $y_delta);}
  elsif ($have_bits) {printf("<text x=\"%d\" y=\"%d\">bits: </text>",54,$max_y + $y_off - 24 + $y_delta);}

  $del = 10;
  for ($i=0; $i<$last ; $i++) {
    $optstr = sprintf("stroke-width=\"1.5\" stroke=\"%s\"",$line_colors[$i]);
    newline($optstr);
    #    linetype(i);*/
    move($xpos[$i]-48,$max_y + $y_off - $ypos[$i] + $y_delta);
    draw($xpos[$i]+12,$max_y + $y_off - $ypos[$i] + $y_delta);
    clsline(1000,1000,10000);

    if ($have_zdb) {
      if ($i==4) {$numstr = sprintf("&gt;%.1lg",$elinval[3]);}
      else {$numstr = sprintf("&lt;%.1lg",$elinval[$i]);}
    }
    elsif ($have_bits) {
      if ($i==4) {$numstr = sprintf("&lt;%.1lf",$blinval[3]);}
      else {$numstr = sprintf("&gt;=%.1lf",$blinval[$i]);}
    }
    else {
      if ($i==3) {$numstr = sprintf("&lt;%d",$ilinval[3]);}
      else {$numstr = sprintf("&gt;%d",$ilinval[$i]);}
    }

    printf("<text align=\"center\" x=\"%d\" y=\"%d\">%s</text>\n",$xpos[$i] + 18, $max_y + $y_off - $ypos[$i] + $y_delta + 4,$numstr);
  }
}

# void linetype(int type)
sub linetype
{
  my $type = shift;
  printf(" stroke=\"%s\"",$line_colors[$type]);
}

#void closeplt()
sub closeplt
{
  print("</svg>\n");
}

#void opnline(int s, double bits)
sub opnline
{
  my ($s, $bits) = @_;

  my $e_val;

  if ($have_zdb) {
    $e_val = bit_to_E($bits);
    printf("<!-- score: %d; bits: %.1g; E(): %.1g -->\n",$s,$bits,$e_val);
    print("<path ");
    if ($e_val < $elinval[0]) {linetype(0);}
    elsif ($e_val < $elinval[1]) {linetype(1);}
    elsif ($e_val < $elinval[2]) {linetype(2);}
    elsif ($e_val < $elinval[3]) {linetype(3);}
    else {linetype(4);}
  }
  elsif ($have_bits) {
    printf("<!-- score: %d; bits: %.1g -->\n",$s,$bits);
    print("<path ");
    if ($bits >= $blinval[0]) {linetype(0);}
    elsif ($bits >= $blinval[1]) {linetype(1);}
    elsif ($bits >= $blinval[2]) {linetype(2);}
    elsif ($bits >= $blinval[3]) {linetype(3);}
    else {linetype(4);}
  }
  else {
    printf("<!-- score: %d -->\n",$s);
    print("<path ");
    if ($s > $ilinval[0]) {linetype(0);}
    elsif ($s> $ilinval[1]) {linetype(1);}
    elsif ($s> $ilinval[2]) {linetype(2);}
    else {linetype(3);}
  }

  print(" d=\"");
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

#void cal_coord(int n0, int n1, long *a_start0, long *a_stop0, long *a_start1, long *a_stop1 )
sub cal_coord {}

