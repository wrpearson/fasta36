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

#define SX(x) (int)((double)(x)*fxscal+fxoff+24)
sub SX {
  my $xx = shift;
  return int($xx*$fxscal+$fxoff+24);
}

#define SY(y) (int)((double)(y)*fyscal+fyoff+48)
sub SY {
  my $yy = shift;
  return int($yy*$fyscal+$fyoff+84);
}

# alignment lines: black blue cyan green lt_green */
my @rlincol=(0.0,0.0,0.0,0.45,0.0);
my @glincol=(0.0,0.0,0.5,0.30,1.0);
my @blincol=(0.0,0.8,0.5,0.15,0.0);

my @line_colors=qw(black blue cyan green lightgreen);
my @block_colors = qw( slategrey lightgreen lightblue pink cyan tan gold plum mediumplum );

# domain blocks: grey blue cyan green lt_green
my @rblk_col=(0.33, 0.0, 0.0, 0.45, 0.0);
my @gblk_col=(0.33, 0.0, 0.5, 0.30, 1.0);
my @bblk_col=(0.33, 0.8, 0.5, 0.15, 0.0);

# void openplt(long n0, long n1, int sq0off, int sq1off, char *xtitle, char *ytitle)
sub openplt {
  my ($n0, $n1, $sq0off, $sq1off, $xtitle, $ytitle, $x_annot_r, $y_annot_r, $have_zdb, $have_bits) = @_;

  if ($lvstr) {
    @elinval = split(/\s+/,$lvstr);
  }
  elsif ($ENV{LINEVAL}) {
    @elinval = split(/\s+/,$ENV{LINEVAL});
  }
	
## 8.5 x 11 paper is 612 pt x 792 pt -- important to stay on one page, with comments
## max_x, max_y are set to 540 pt -- 7.5 in, 9pt (0.75 in) margins, leaving little extra space
## if comments are provided, max_x, max_y must be reduced (max_x for space, max_y to keep things square)
##

  my ($xbound, $ybound) = ($max_x + 24, $max_y + 72);
  if ($x_annot_r) {$ybound += 64;}
  if ($y_annot_r) {
    $xbound += 100;
    $max_x -= $max_x/10;
    $max_y -= $max_x/10;
  }

  print("%!PS-Adobe-2.0\n");
  print("%%Creator: plalign\n");
  print("%%CreationDate: %s","2012-01-01");
  print("%%DocumentFonts: Courier\n");
  print("%%Pages: 1\n");
  print("%%BoundingBox: 18 18 $xbound $ybound\n");
  print("%%EndComments\n");
  print("%%EndProlog\n");
  print("%%Page: 1 1\n");
  print("/Courier findfont 14 scalefont setfont\n");
  print("/vcprint { gsave 90 rotate dup stringwidth pop 2 div neg 0 rmoveto\n");
  print("show newpath stroke grestore } def\n");
  print("/vprint { gsave 90 rotate\n");
  print("show newpath stroke grestore } def\n");
  print("/hcprint { gsave dup stringwidth pop 2 div neg 0 rmoveto\n");
  print("show newpath stroke grestore } def\n");
  print("/hrprint { gsave dup stringwidth pop neg 0 rmoveto\n");
  print("show newpath stroke grestore } def\n");
  print("/hprint { gsave show newpath stroke grestore } def\n");
  # % x y w h RT - % draw a rectangle size w h at x y
  print("/RT { [ ] 0 setdash newpath 4 -2 roll moveto dup 0 exch rlineto exch 0 rlineto neg 0 exch rlineto closepath }  def \n");

  ($pmaxx, $pmaxy) = ($n0, $n1);

# $max_x, $max_y define the maximum plotting area
# the actual bounding box/view area will be larger if annotation comments are available

  $fxscal = ($max_x-1)/$n1;
  $fyscal = ($max_y-1)/$n0;

  if ($fxscal > $fyscal) {$fxscal = $fyscal;}
  else {$fyscal = $fxscal;}

  if ($fyscal * $n0 < $max_y/5.0) {
    $fyscal = ($max_y-1)/($n0*5.0);
  }

  $fxscal *= 0.9; $fxoff = ($max_x-1)/11.0;
  $fyscal *= 0.9; $fyoff = ($max_y-1)/11.0;

  printf("%% openplt - frame - %ld %ld\n", $n0, $n1);

  # draw the plot frame
  linetype(0);
  print("gsave\n");
  print("currentlinewidth 1.5 mul setlinewidth\n");
  newline();
  move(SX(0),SY(0));
  draw(SX(0),SY($n1+1));
  draw(SX($n0+1),SY($n1+1));
  draw(SX($n0+1),SY(0));
  draw(SX(0),SY(0));
  clsline($n0,$n1,100000);
  print("grestore\n");

  my $n_div = 11;
  xaxis($n0,$sq1off, $xtitle, $n_div);

  $n_div = 21 unless $n0 == $n1;
  yaxis($n1,$sq0off, $ytitle, $n_div);
  legend($have_zdb, $have_bits);

  print("%% openplt done\n");

  if ($x_annot_r) {xgrid($x_annot_r, $n0, $sq0off, $n1, $sq1off);}
  if ($y_annot_r) {ygrid($y_annot_r, $n0, $sq0off, $n1, $sq1off);}
}

#void drawdiag(long n0, long n1)
sub drawdiag  {
  my ($n0, $n1) = @_;
  linetype(0);
  printf("%% drawdiag %ld %ld\n",$n0, $n1);
  print("gsave\n");
  print("currentlinewidth 1.5 mul setlinewidth\n");
  newline();
  move(SX(0),SY(0));
  draw(SX($n0+1),SY($n1+1));
  clsline($n0,$n1,10000);
  print("grestore\n");
  print("%% drawdiag done\n");
}

# tick array - values */
my @tarr = (10,20,50,100,200,500,1000,2000,5000,10000,20000,50000,100000,200000,500000,1000000);

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
    if (($jm = $n/$tarr[$i])<$n_div) {goto found;}
  }
  $i=scalar(@tarr)-1;
  $jm = $n/$tarr[$i];
 found:
  # js is the start of the value - modify to accomodate offset */
  $js = $tarr[$i];

  # jo is the offset */
  $jo = ($offset-1) % $tarr[$i];	# figure out offset in tarr[i] increments */

  # jl is the label */
  $jl = ($offset-1)/$tarr[$i];	# figure out offset in tarr[i] increments */
  $jl *= $tarr[$i];

  newline();
  for ($i=1; $i<=$jm; $i++) {
    move(SX($i*$js - $jo),SY(0));
    draw(SX($i*$js - $jo),SY(0)-$tick);
  }
  clsline($n,$n,10000);

  $numstr = sprintf("%ld",$js + $jl );
  $num_len = length($numstr);

  if ($num_len > 4) {
    move(SX($js-$jo),SY(0)-$tick-16);
    printf("(%s) hcprint\n",$numstr);

    $numstr=sprintf("%ld",$jm*$js+$jl);
    move(SX($jm*$js-$jo),SY(0)-$tick-16);
    printf("(%s) hcprint\n",$numstr);
  }
  else {	# put in all the axis values */
    for ($i=1; $i<=$jm; $i++) {
      $numstr=sprintf("%ld",$i*$js+$jl);
      move(SX($i*$js-$jo),SY(0)-$tick-16);
      printf("(%s) hcprint\n",$numstr);
    }
  }

  print("newpath\n");
  move(SX($n/2),SY(0)-$tick-30);

#  for (bp = strchr(title,'('); (bp!=NULL); bp = strchr(bp+1,'(')) *bp=' ';
#  for (bp = strchr(title,')'); (bp!=NULL); bp = strchr(bp+1,')')) *bp=' ';
  $title =~ s/\(/ /g;
  $title =~ s/\)/ /g;
  printf("(%s) hcprint\n",$title);
}
		
# void yaxis(long n, int offset, char *title)
sub yaxis {
  my ($n, $offset, $title, $n_div) = @_;

  my  ($i, $jm, $tick);
  my ($js, $jo, $jl);
  my ($num_len, $numstr);
	
  $tick = 6;

  for ($i=0; $i<@tarr; $i++) {
    if (($jm = $n/$tarr[$i])<$n_div) {goto found;}
  }
  $jm = $n/5000;
  $i= scalear(@tarr)-1;

 found:
  $js = $tarr[$i];

  # $jo is the offset */
  $jo = ($offset-1) % $tarr[$i];	# figure out offset in tarr[i] increments */
  # jl is the label */
  $jl = ($offset-1)/$tarr[$i];	# figure out offset in tarr[i] increments */
  $jl *= $tarr[$i];

  newline();
  for ($i=1; $i<=$jm; $i++) {
    move(SX(0),SY($i*$js-$jo));
    draw(SX(0)-$tick,SY($i*$js-$jo));
  }
  clsline($n,$n,10000);

  $numstr = sprintf("%ld",$js+$jl);

  $num_len = length($numstr);

  if ($num_len > 4) {
    move(SX(0)-$tick-4,SY($js-$jo)-4);
    printf("(%s) hrprint\n",$numstr);

    $numstr = sprintf("%ld",$jm*$js+$jl);
    move(SX(0)-$tick-4,SY($jm*$js-$jo)-4);
    printf("(%s) hrprint\n",$numstr);
  }
  else {
    for ($i=1; $i<=$jm; $i++) {
      $numstr = sprintf("%ld",$i*$js+$jl);
      move(SX(0)-$tick-4,SY($i*$js-$jo)-4);
      printf("(%s) hrprint\n",$numstr);
    }
  }

  move(SX(0)-$tick-42,SY($n/2));
#  for (bp = strchr(title,'('); (bp!=NULL); bp = strchr(bp+1,'(')) *bp=' ';
#  for (bp = strchr(title,')'); (bp!=NULL); bp = strchr(bp+1,')')) *bp=' ';
  $title =~ s/\(/\\(/g;
  $title =~ s/\)/\\)/g;
  printf("(%s) vcprint\n",$title);

}

sub xgrid {
  my ($annot_arr_r, $n0, $sq0_off, $n1, $sq1_off) = @_;

  my $sq_off = $sq0_off;

  my $show_block = 1;
  my $text_offset = 8;
  if ($show_block) {$text_offset = 24;}
  my $color = 1;

  print("%% xgrid: $n0 $n1\n");
  print("gsave\n");
  print("/Courier findfont 11 scalefont setfont\n");
  print("currentlinewidth 0.5 mul setlinewidth\n");
  for my $annot ( @$annot_arr_r) {
    next unless $annot->{beg} >= $sq_off;
    next if ($annot->{end} > $sq_off + $n0 - 1);
    last if ($annot->{beg} > $sq_off + $n0 - 1);
    newline();
    print("0.33 0.33 0.33 setrgbcolor\n");
    move(SX($annot->{beg}-$sq_off),SY(0));
    print("[3 6] 0 setdash\n");
    draw(SX($annot->{beg}-$sq_off),SY($n1));
    clsline();
    newline();
    print("0.33 0.33 0.33 setrgbcolor\n");
    move(SX($annot->{end}-$sq_off),SY(0));
    print("[6 3] 0 setdash\n");
    draw(SX($annot->{end}-$sq_off),SY($n1));
    clsline();

    # show rotated label
    my $xpos = SX(($annot->{end} - $annot->{beg})/2 + $annot->{beg} - $sq_off) + 4;
    my $ypos = SY($n1) + $text_offset;
    # printf("<text x=\"0\" y=\"0\" text-anchor=\"left\" transform=\"translate($xpos, $ypos) rotate(-90,0,0)\">%s</text>\n",$annot->{sdescr});
    if ($show_block) {
      draw_block(SX($annot->{beg} - $sq_off), SY($n1) + 6, 
		 SX($annot->{end} - $sq_off)-SX($annot->{beg} - $sq_off),12,
		 $annot_names{$annot->{sname}});
    }
    move($xpos, $ypos);
    my $str = $annot->{sdescr};
    $str =~ s/\(/\\(/g;
    $str =~ s/\)/\\)/g;
    print "($str) vprint\n";
  }
  print("grestore\n");
}

sub ygrid {
  my ($annot_arr_r, $n0, $sq0_off, $n1, $sq1_off) = @_;

  my $sq_off = $sq1_off;

  my $show_block = 1;

  my $text_offset = 8;
  if ($show_block) {$text_offset = 24;}
  my $color=4;

  print("gsave\n");
  print("/Courier findfont 11 scalefont setfont\n");
  print("currentlinewidth 0.5 mul setlinewidth\n");
  for my $annot ( @$annot_arr_r) {
    next unless $annot->{beg} >= $sq_off;
    next if ($annot->{end} > $sq_off + $n1 - 1);
    last if ($annot->{beg} > $sq_off + $n1 - 1);
    newline();
    print("0.33 0.33 0.33 setrgbcolor\n");
    move(SX(0), SY($annot->{beg}-$sq_off));
    print("[3 6] 0 setdash\n");
    draw(SX($n0), SY($annot->{beg}-$sq_off));
    clsline();
    newline();
    move(SX(0), SY($annot->{end}-$sq_off));
    print("[6 3] 0 setdash\n");
    draw(SX($n0), SY($annot->{end}-$sq_off));
    clsline();

    my $xpos = SX($n0) + $text_offset;
    my $ypos = SY(($annot->{end} - $annot->{beg})/2 + $annot->{beg} - $sq_off) - 3;

    if ($show_block) {
      draw_block(SX($n0)+6, SY($annot->{beg} - $sq_off), 12,
		 SY($annot->{end} - $sq_off)-SY($annot->{beg} - $sq_off),
		 $annot_names{$annot->{sname}});
    }

#    printf("<text x=\"$xpos\" y=\"$ypos\" text-anchor=\"left\">%s</text>\n",$annot->{sdescr});
    move($xpos, $ypos);
    my $str = $annot->{sdescr};
    $str =~ s/\(/\\(/g;
    $str =~ s/\)/\\)/g;
    print "($str) hprint\n";
  }
  print("grestore\n");
}

sub draw_block {
  my ($x, $y, $w, $h, $color) = @_;

  $color = ($color % scalar(@block_colors));
  my $rgb = $color_names{$block_colors[$color]};

# color is index into @[rgb]blk_color 
  print "gsave\n";

  printf("%.3f %.3f %.3f setrgbcolor\n",
	 $rgb->[0]/255,$rgb->[1]/255,$rgb->[2]/255);

  printf "%d %d %d %d RT\n",$x+1,$y+1,$w-2,$h-2;
  print "fill\n";
  print "stroke\n";

  print "1.0 1.0 1.0 setrgbcolor\n";
  print "$x $y $w $h RT\n";
  print "stroke\n";
  print "grestore\n";
}

# void legend()
sub legend {
  my ($have_zdb, $have_bits) = @_;

  my ($i, $last);
  my ($ixp, $iyp);
  my $numstr;
  my @xpos=(144,144,288,288,432);
  my @ypos=(36,18,36,18,27);

  if ($have_zdb || $have_bits) {$last = 5;}
  else {$last = 4;}

  move(60,27+18);
  if ($have_zdb)  {draw_str("E():  ");}
  elsif ($have_bits) {draw_str("Bits:  ");}

  for ($i=0; $i<$last ; $i++) {
    print("gsave currentlinewidth 1.5 mul setlinewidth\n");
    newline();
    linetype($i);
    move($xpos[$i]-36,$ypos[$i]+18);
    draw($xpos[$i]+24,$ypos[$i]+18);
    clsline(1000,1000,10000);
    print("grestore\n");
    move($xpos[$i]+36,$ypos[$i]-4+18);
    if ($have_zdb) {
      if ($i==4) {$numstr = sprintf(">%.1lg",$elinval[3]);}
      else {$numstr = sprintf("<%.1lg",$elinval[$i]);}
    }
    elsif ($have_bits) {
      if ($i==4) {$numstr = sprintf("<%.1lf",$blinval[3]);}
      else {$numstr = sprintf(">=%.1lf",$blinval[$i]);}
    }
    else {
      if ($i==3) {$numstr = sprintf("<%d",$ilinval[3]);}
      else {$numstr = sprintf(">%d",$ilinval[$i]);}
    }
    printf("(%s) hprint\n",$numstr);
  }
}

#void linetype(type)
sub linetype {
  my $type = shift;

  my $rgb_name = $line_colors[$type];
  my $rgb = $color_names{$rgb_name};

  printf("%5.3f %5.3f %5.3f setrgbcolor\n",
	 $rgb->[0]/256, $rgb->[1]/256, $rgb->[2]/256);
}

#void closeplt()
sub closeplt {
  print("%%Trailer\n");
  print("showpage\n");
  print("%%EOF\n");
}

# void opnline(int s, double bits)
sub opnline {
  my ($s, $bits) = shift;

  my $e_val;

  if ($have_zdb) {
    $e_val = bit_to_E($bits);
    if ($e_val < $elinval[0]) {linetype(0);}
    elsif ($e_val < $elinval[1]) {linetype(1);}
    elsif ($e_val < $elinval[2]) {linetype(2);}
    elsif ($e_val < $elinval[3]) {linetype(3);}
    else {linetype(4);}
  }
  elsif ($have_bits) {
    if ($bits >= $blinval[0]) {linetype(0);}
    elsif ($bits >= $blinval[1]) {linetype(1);}
    elsif ($bits >= $blinval[2]) {linetype(2);}
    elsif ($bits >= $blinval[3]) {linetype(3);}
    else {linetype(4);}
  }
  else {
    if ($s > $ilinval[0]) {linetype(0);}
    elsif ($s > $ilinval[1]) {linetype(1);}
    elsif ($s > $ilinval[2]) {linetype(2);}
    else {linetype(3);}
  }

  print("newpath\n");
}

# void newline()
sub newline {
  print("0 0 0 setrgbcolor\n newpath\n");
}

# void clsline(long x,long y,int s)
sub clsline {
  print("stroke\n");
}

# void move(int x, int y)
sub move {
  my ($xx, $yy) = @_;

  printf("%d %d moveto\n",$xx,$yy);
}

# void sxy_move(int x, int y)
sub sxy_move {
  my ($x, $y) = @_;
  printf("%d %d moveto\n",SX($x),SY($y));
}

# void draw(int x, int y)
sub draw {
  my ($x,$y) = @_;
  printf("%d %d lineto\n",$x,$y);
}

# void sxy_draw(int x, int y)
sub sxy_draw {
  my ($x,$y) = @_;
  printf("%d %d lineto\n",SX($x),SY($y));
}

#void draw_str(char *str)
sub draw_str
{
  my $str = shift;

#  for (bp = strchr(str,'('); (bp!=NULL); bp = strchr(bp+1,'(')) *bp=' ';
#  for (bp = strchr(str,')'); (bp!=NULL); bp = strchr(bp+1,')')) *bp=' ';

  $str =~ s/\(/\\(/g;
  $str =~ s/\)/\\)/g;

  printf("(%s) show\n",$str);
}

#void cal_coord(int n0, int n1, long *a_start0, long *a_stop0, long *a_start1, long *a_stop1 )
sub cal_coord {}
