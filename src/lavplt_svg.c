/* functions/variables for postscript plots */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

#include "lav_defs.h"

#define SX(x) (int)((double)(x)*fxscal+fxoff+6)
#define SY(y) (max_y + 24 - (int)((double)(y)*fyscal+fyoff))

/* black blue cyan green lt_green */
char *line_color[]={"black","blue","brown","green","red"};

void
openplt(long n0, long n1, int sq0off, int sq1off, 
	char *xtitle, char *ytitle)
{
  char *sptr;
  time_t tt;
  int tick = 6;

  tt = time(NULL);

  if (strlen(lvstr)>0) {
    sscanf(lvstr,"%lg %lg %lg",&elinval[0],&elinval[1],&elinval[2]);
  }
  else if ((sptr=getenv("LINEVAL"))!=NULL && strlen(sptr)>0) {
    sscanf(sptr,"%lg %lg %lg",&elinval[0],&elinval[1],&elinval[2]);
  }
	
  printf("<?xml version=\"1.0\" standalone=\"no\"?>\n");
  printf("<!DOCTYPE svg PUBLIC \"-//W3C//DTD SVG 1.1//EN\" \n");
  printf("\"http://www.w3.org/Graphics/SVG/1.1/DTD/svg11.dtd\">\n\n");

  printf("<svg width=\"564\" height=\"612\" version=\"1.1\"\n");
  printf("xmlns=\"http://www.w3.org/2000/svg\" xmlns:xlink=\"http://www.w3.org/1999/xlink\">\n\n");

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

  /*   printf("currentlinewidth 1.5 mul setlinewidth\n"); */

  printf("<rect x=\"%d\" y=\"%d\" width=\"%d\" height=\"%d\"\n",
	 SX(0),SY(n1+1), SX(n0+1)-SX(0), SY(0) - SY(n1+1));
  printf("stroke=\"black\" stroke-width=\"2.0\" fill=\"none\" />\n");

  xaxis(n0,sq1off, xtitle);
  yaxis(n1,sq0off, ytitle);
  legend();
}
	
void
drawdiag(n0,n1)
	long n0, n1;
{
  
	/* 	printf("currentlinewidth 1.5 mul setlinewidth\n"); */
	newline("stroke=\"black\" stroke-width=\"1.5\"");
	move(SX(0),SY(0));
	draw(SX(n0+1),SY(n1+1));
	clsline(n0,n1,10000);
}

/* tick array - values */
int tarr[] = {10,20,50,100,200,500,1000,2000,5000,10000,20000,50000,100000,200000,500000,1000000};
#define MAX_INTERVAL 1000000
int ntarr = sizeof(tarr)/sizeof(int);

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

  newline("stroke=\"black\" stroke-width=\"1.5\"");
  for (i=1; i<=jm; i++) {
    move(SX((long)i*js - jo),SY(0));
    draw(SX((long)i*js - jo),SY(0)+tick);
  }
  clsline(n,n,10000);

  sprintf(numstr,"%ld",js + jl );
  num_len = strlen(numstr);
  if (num_len > 4) {

    printf("<text x=\"%d\" y=\"%d\" text-anchor=\"middle\">%s</text>\n",SX((long)js-jo),SY(0)+tick+16,numstr);

    sprintf(numstr,"%ld",jm*js+jl);
    printf("<text x=\"%d\" y=\"%d\" text-anchor=\"middle\">%s</text>\n",SX((long)jm*js-jo),SY(0)+tick+16,numstr);
  }
  else {
    for (i=1; i<=jm; i++) {
      sprintf(numstr,"%ld",i*js+jl);
      printf("<text x=\"%d\" y=\"%d\" text-anchor=\"middle\">%s</text>\n",SX((long)i*js-jo),SY(0)+tick+16,numstr);
    }
  }

  printf("<text x=\"%d\" y=\"%d\" text-anchor=\"middle\">%s</text>\n",SX(n/2),SY(0)-tick+42, title);
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

  newline("stroke=\"black\" stroke-width=\"1.5\"");
  for (i=1; i<=jm; i++) {
    move(SX(0),SY((long)i*js-jo));
    draw(SX(0)-tick,SY((long)i*js-jo));
  }
  clsline(n,n,10000);

  sprintf(numstr,"%ld",js+jl);
  num_len = strlen(numstr);
  if (num_len > 4) {
    printf("<text x=\"%d\" y=\"%d\" text-anchor=\"end\">%s</text>\n",SX(0)-tick-4,SY(js-jo)+4,numstr);

    sprintf(numstr,"%ld",(long)jm*js+jl);
    printf("<text x=\"%d\" y=\"%d\" text-anchor=\"end\">%s</text>\n",SX(0)-tick-4,SY((long)jm*js-jo)+4,numstr);
  }
  else {
    for (i=1; i<=jm; i++) {
      sprintf(numstr,"%ld",(long)i*js+jl);
      printf("<text x=\"%d\" y=\"%d\" text-anchor=\"end\">%s</text>\n",SX(0)-tick-4,SY((long)i*js-jo)+4,numstr);
    }
  }
  /* make a path for the label */

  printf("<g transform=\"rotate(-90, 142, 142)\">\n");
  printf("<text x=\"0\" y=\"16\" text-anchor=\"middle\">%s</text>\n",title);
  printf("</g>\n");
}

void
legend()
{
  int i, last, del;
  int ixp, iyp;
  char numstr[20];
  char optstr[128];
  int xpos[]={144,144,288,288,432};
  int ypos[]={36,18,36,18,27};

  int tick = 6;

  if (have_zdb || have_bits) last = 5;
  else last = 4;

  if (have_zdb)  printf("<text x=\"%d\" y=\"%d\">E(): </text>",54,max_y + 66 - 27);
  else if (have_bits) printf("<text x=\"%d\" y=\"%d\">bits: </text>",54,max_y + 66 - 27);

  del = 10;
  for (i=0; i<last ; i++) {
    
    sprintf(optstr,"stroke-width=\"1.5\" stroke=\"%s\"",line_color[i]);
    newline(optstr);
    /*    linetype(i);*/
    move(xpos[i]-12,max_y + 66 - ypos[i]);
    draw(xpos[i]+48,max_y + 66 - ypos[i]);
    clsline(1000,1000,10000);

    if (have_zdb) {
      if (i==4) sprintf(numstr,"&gt;%.1lg",elinval[3]);
      else sprintf(numstr,"&lt;%.1lg",elinval[i]);
    }
    else if (have_bits) {
      if (i==4) sprintf(numstr,"&lt;%.1lf",blinval[3]);
      else sprintf(numstr,"&gt;=%.1lf",blinval[i]);
    }
    else {
      if (i==3) sprintf(numstr,"&lt;%d",ilinval[3]);
      else sprintf(numstr,"&gt;%d",ilinval[i]);
    }

    printf("<text align=\"center\" x=\"%d\" y=\"%d\">%s</text>\n",xpos[i]+54, max_y + 66 - ypos[i],numstr);
  }
}

void
linetype(type)
     int type;
{
  printf(" stroke=\"%s\"",line_color[type]);
}

void
closeplt()
{
  printf("</svg>\n");
}

void
opnline(int s, double bits)
{
  double e_val;

  if (have_zdb) {
    e_val = bit_to_E(bits);
    printf("<!-- score: %d; bits: %.1g; E(): %.1g -->\n",s,bits,e_val);
    printf("<path ");
    if (e_val < elinval[0]) linetype(0);
    else if (e_val < elinval[1]) linetype(1);
    else if (e_val < elinval[2]) linetype(2);
    else if (e_val < elinval[3]) linetype(3);
    else linetype(4);
  }
  else if (have_bits) {
    printf("<!-- score: %d; bits: %.1g -->\n",s,bits);
    printf("<path ");
    if (bits >= blinval[0]) linetype(0);
    else if (bits >= blinval[1]) linetype(1);
    else if (bits >= blinval[2]) linetype(2);
    else if (bits >= blinval[3]) linetype(3);
    else linetype(4);
  }
  else {
    printf("<!-- score: %d -->\n",s);
    printf("<path ");
    if (s > ilinval[0]) linetype(0);
    else if (s> ilinval[1]) linetype(1);
    else if (s> ilinval[2]) linetype(2);
    else linetype(3);
  }

  printf(" d=\"");
}

void
newline(char *options)
{
  if (options != NULL)  {
      printf("<path %s d=\"",options);
  }
  else printf("<path stroke=\"black\" d=\"");
}

void
clsline(long x, long y, int s)
{
  printf("\" fill=\"none\" />\n");
}

void
move(int x, int y)
{
  printf(" M %d %d",x,y);
}

void
sxy_move(int x, int y)
{
  move(SX(x), SY(y));
}

void
draw(int x, int y)
{
  printf(" L %d %d",x,y);
}

void
sxy_draw(int x, int y)
{
  draw(SX(x),SY(y));
}

void
draw_str(char *str) {
  char *bp;

  for (bp = strchr(str,'('); (bp!=NULL); bp = strchr(bp+1,'(')) *bp=' ';
  for (bp = strchr(str,')'); (bp!=NULL); bp = strchr(bp+1,')')) *bp=' ';

  printf("(%s) show\n",str);
}

void
draw_sstr(char *str) {

  printf("(%s) show\n",str);
}

void cal_coord(int n0, int n1, 
	       long *a_start0, long *a_stop0, long *a_start1, long *a_stop1 )
{}
