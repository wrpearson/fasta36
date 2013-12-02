/* functions/variables for postscript plots */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

#include "lav_defs.h"

#define SX(x) (int)((double)(x)*fxscal+fxoff+24)
#define SY(y) (int)((double)(y)*fyscal+fyoff+48)

/* black blue cyan green lt_green */
float rlincol[]={0.0,0.0,0.0,0.45,0.0};

float glincol[]={0.0,0.0,0.5,0.30,1.0};
float blincol[]={0.0,0.8,0.5,0.15,0.0};

void
openplt(long n0, long n1, int sq0off, int sq1off, 
	char *xtitle, char *ytitle)
{
  char *getenv(), *sptr;
  time_t tt;

  tt = time(NULL);

  if (strlen(lvstr)>0) {
    sscanf(lvstr,"%lg %lg %lg",&elinval[0],&elinval[1],&elinval[2]);
  }
  else if ((sptr=getenv("LINEVAL"))!=NULL && strlen(sptr)>0) {
    sscanf(sptr,"%lg %lg %lg",&elinval[0],&elinval[1],&elinval[2]);
  }
	
  printf("%%!PS-Adobe-2.0\n");
  printf("%%%%Creator: plalign\n");
  printf("%%%%CreationDate: %s",ctime(&tt));
  printf("%%%%DocumentFonts: Courier\n");
  printf("%%%%Pages: 1\n");
  printf("%%%%BoundingBox: 18 18 564 588\n");
  printf("%%%%EndComments\n");
  printf("%%%%EndProlog\n");
  printf("%%%%Page: 1 1\n");
  printf("/Courier findfont 14 scalefont setfont\n");
  printf("/vcprint { gsave 90 rotate dup stringwidth pop 2 div neg 0 rmoveto\n");
  printf("show newpath stroke grestore } def\n");
  printf("/hcprint { gsave dup stringwidth pop 2 div neg 0 rmoveto\n");
  printf("show newpath stroke grestore } def\n");
  printf("/hrprint { gsave dup stringwidth pop neg 0 rmoveto\n");
  printf("show newpath stroke grestore } def\n");
  printf("/hprint { gsave show newpath stroke grestore } def\n");

  pmaxx = n0;
  pmaxy = n1;

  fxscal = (double)(max_x-1)/(double)(n1);
  fyscal = (double)(max_y-1)/(double)(n0);

  if (fxscal > fyscal) fxscal = fyscal;
  else fyscal = fxscal;

  if (fyscal * n0 < (double)max_y/5.0) 
    fyscal = (double)(max_y-1)/((double)(n0)*5.0);

  fxscal *= 0.9; fxoff = (double)(max_x-1)/11.0;
  fyscal *= 0.9; fyoff = (double)(max_y-1)/11.0;

  printf("%% openplt - frame - %ld %ld\n", n0, n1);
  linetype(0);
  printf("gsave\n");
  printf("currentlinewidth 1.5 mul setlinewidth\n");
  newline();
  move(SX(0),SY(0));
  draw(SX(0),SY(n1+1));
  draw(SX(n0+1),SY(n1+1));
  draw(SX(n0+1),SY(0));
  draw(SX(0),SY(0));
  clsline(n0,n1,100000);
  printf("grestore\n");
  xaxis(n0,sq1off, xtitle);
  yaxis(n1,sq0off, ytitle);
  legend();
  printf("%% openplt done\n");
}
	
void
drawdiag(long n0, long n1)
{
  
	linetype(0);
	printf("%% drawdiag %ld %ld\n",n0, n1);
	printf("gsave\n");
	printf("currentlinewidth 1.5 mul setlinewidth\n");
	newline();
	move(SX(0),SY(0));
	draw(SX(n0+1),SY(n1+1));
	clsline(n0,n1,10000);
	printf("grestore\n");
	printf("%% drawdiag done\n");
}

/* tick array - values */
static int tarr[] = {10,20,50,100,200,500,1000,2000,5000,10000,20000,50000,100000,200000,500000,1000000};
static int ntarr = sizeof(tarr)/sizeof(int);

void
xaxis(long n, int offset, char *title)
{
  int i, jm, tick;
  long js, jo, jl;
  int num_len;
  char numstr[20],*bp;

  tick = 6;

  /* search for the correct increment for the tick array */
  for (i=0; i<ntarr; i++) {
    /* seek to divide into 20 or fewer divisions */
    if ((jm = n/tarr[i])<21) goto found;
  }
  i=ntarr-1;
  jm = n/tarr[i];
 found:
  /* js is the start of the value - modify to accomodate offset */
  js = tarr[i];

  /* jo is the offset */
  jo = offset%tarr[i];	/* figure out offset in tarr[i] increments */

  /* jl is the label */
  jl = offset/tarr[i];	/* figure out offset in tarr[i] increments */
  jl *= tarr[i];

  newline();
  for (i=1; i<=jm; i++) {
    move(SX((long)i*js - jo),SY(0));
    draw(SX((long)i*js - jo),SY(0)-tick);
  }
  clsline(n,n,10000);


  sprintf(numstr,"%ld",js + jl );
  num_len = strlen(numstr);

  if (num_len > 4) {
    move(SX(js-jo),SY(0)-tick-16);
    printf("(%s) hcprint\n",numstr);

    sprintf(numstr,"%ld",jm*js+jl);
    move(SX((long)jm*js-jo),SY(0)-tick-16);
    printf("(%s) hcprint\n",numstr);
  }
  else {	/* put in all the axis values */
    for (i=1; i<=jm; i++) {
      sprintf(numstr,"%ld",i*js+jl);
      move(SX(i*js-jo),SY(0)-tick-16);
      printf("(%s) hcprint\n",numstr);
    }
  }

  printf("newpath\n");
  move(SX(n/2),SY(0)-tick-30);

  for (bp = strchr(title,'('); (bp!=NULL); bp = strchr(bp+1,'(')) *bp=' ';
  for (bp = strchr(title,')'); (bp!=NULL); bp = strchr(bp+1,')')) *bp=' ';
  printf("(%s) hcprint\n",title);
}
		
void
yaxis(long n, int offset, char *title)
{
  int i, jm, tick;
  long js, jo, jl;
  int num_len;
  char numstr[20],*bp;
	
  tick = 6;

  for (i=0; i<ntarr; i++) {
    if ((jm = n/tarr[i])<21) goto found;
  }
  jm = n/5000l;
  i=ntarr-1;

 found:
  js = (long)tarr[i];

  /* jo is the offset */
  jo = offset%tarr[i];	/* figure out offset in tarr[i] increments */
  /* jl is the label */
  jl = offset/tarr[i];	/* figure out offset in tarr[i] increments */
  jl *= tarr[i];

  newline();
  for (i=1; i<=jm; i++) {
    move(SX(0),SY((long)i*js-jo));
    draw(SX(0)-tick,SY((long)i*js-jo));
  }
  clsline(n,n,10000);


  sprintf(numstr,"%ld",js+jl);

  num_len = strlen(numstr);

  if (num_len > 4) {
    move(SX(0)-tick-4,SY(js-jo)-4);
    printf("(%s) hrprint\n",numstr);

    sprintf(numstr,"%ld",(long)jm*js+jl);
    move(SX(0)-tick-4,SY((long)jm*js-jo)-4);
    printf("(%s) hrprint\n",numstr);
  }
  else {
    for (i=1; i<=jm; i++) {
      sprintf(numstr,"%ld",(long)i*js+jl);
      move(SX(0)-tick-4,SY((long)i*js-jo)-4);
      printf("(%s) hrprint\n",numstr);
    }
  }

  move(SX(0)-tick-42,SY(n/2));
  for (bp = strchr(title,'('); (bp!=NULL); bp = strchr(bp+1,'(')) *bp=' ';
  for (bp = strchr(title,')'); (bp!=NULL); bp = strchr(bp+1,')')) *bp=' ';
  printf("(%s) vcprint\n",title);

}

void
legend()
{
  int i, last, del;
  int ixp, iyp;
  char numstr[20];
  int xpos[]={144,144,288,288,432};
  int ypos[]={36,18,36,18,27};

  if (have_zdb || have_bits) last = 5;
  else last = 4;

  move(72,27);
  if (have_zdb)  draw_sstr("E\\(\\):  ");
  else if (have_bits) draw_str("Bits:  ");

  del = 10;
  for (i=0; i<last ; i++) {
    printf("gsave currentlinewidth 1.5 mul setlinewidth\n");
    newline();
    linetype(i);
    move(xpos[i],ypos[i]);
    draw(xpos[i]+60,ypos[i]);
    clsline(1000,1000,10000);
    printf("grestore\n");
    move(xpos[i]+72,ypos[i]-4);
    if (have_zdb) {
      if (i==4) sprintf(numstr,">%.1lg",elinval[3]);
      else sprintf(numstr,"<%.1lg",elinval[i]);
    }
    else if (have_bits) {
      if (i==4) sprintf(numstr,"<%.1lf",blinval[3]);
      else sprintf(numstr,">=%.1lf",blinval[i]);
    }
    else {
      if (i==3) sprintf(numstr,"<%d",ilinval[3]);
      else sprintf(numstr,">%d",ilinval[i]);
    }
    printf("(%s) hprint\n",numstr);
  }
}

void
linetype(type)
     int type;
{
  printf("%5.3f %5.3f %5.3f setrgbcolor\n",
	 rlincol[type],glincol[type],blincol[type]);
}

void
closeplt()
{
  printf("%%%%Trailer\n");
  printf("showpage\n");
  printf("%%%%EOF\n");
}

void
opnline(int s, double bits)
{
  double e_val;

  if (have_zdb) {
    e_val = bit_to_E(bits);
    if (e_val < elinval[0]) linetype(0);
    else if (e_val < elinval[1]) linetype(1);
    else if (e_val < elinval[2]) linetype(2);
    else if (e_val < elinval[3]) linetype(3);
    else linetype(4);
  }
  else if (have_bits) {
    if (bits >= blinval[0]) linetype(0);
    else if (bits >= blinval[1]) linetype(1);
    else if (bits >= blinval[2]) linetype(2);
    else if (bits >= blinval[3]) linetype(3);
    else linetype(4);
  }
  else {
    if (s > ilinval[0]) linetype(0);
    else if (s> ilinval[1]) linetype(1);
    else if (s> ilinval[2]) linetype(2);
    else linetype(3);
  }

  printf("newpath\n");
}

void
newline()
{
  printf("0 0 0 setrgbcolor\n newpath\n");
}

void
clsline(long x,long y,int s)
{
  printf("stroke\n");
}

void
move(int x, int y)
{
  printf("%d %d moveto\n",x,y);
}

void
sxy_move(int x, int y)
{
  x = SX(x);
  y = SY(y);
  printf("%d %d moveto\n",x,y);
}

void
draw(int x, int y)
{
  printf("%d %d lineto\n",x,y);
}

void
sxy_draw(int x, int y)
{
  x = SX(x);
  y = SY(y);
  printf("%d %d lineto\n",x,y);
}

void
draw_str(char *str)
{
  char *bp;

  for (bp = strchr(str,'('); (bp!=NULL); bp = strchr(bp+1,'(')) *bp=' ';
  for (bp = strchr(str,')'); (bp!=NULL); bp = strchr(bp+1,')')) *bp=' ';

  printf("(%s) show\n",str);
}

void
draw_sstr(char *str)
{
  printf("(%s) show\n",str);
}

void cal_coord(int n0, int n1, 
	       long *a_start0, long *a_stop0, long *a_start1, long *a_stop1 )
{}
