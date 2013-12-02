/* lav2ps.c - produce postscript from lav output */

/* $Id: lav2plt.c 625 2011-03-23 17:21:38Z wrp $ */
/* $Revision: 625 $  */

/* 15-May-2009 fix tarr[i] out of range bug */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

#define XTERNAL
#include "lav_defs.h"

int g_n0, g_n1;

void openplt(long, long, int, int, char *, char *);
void closeplt();
void drawdiag(long n0, long n1);
void xaxis(long, int, char *);
void yaxis(long, int, char *);
void legend();
void linetype(int);
void opnline(int s, double bits);
void clsline();
void sxy_move(int, int);
void sxy_draw(int, int);
void draw_str(char *);
void draw_sstr(char *);

int have_bits, have_zdb;
double bit_to_E(double bits);

void get_str(FILE *, char *str, size_t len);
void get_str2(FILE *, char *, size_t, char *, size_t );
void get_seq_info(FILE *file, 
		  char *str0, size_t len0, int *n0_begin, int *n0_end,
		  char *str1, size_t len1, int *n1_begin, int *n1_end);
void del1(char *);
void do_alignment(FILE *, int, int);

/* for getopt() */
#ifdef UNIX
#include <unistd.h>
#else
extern int optind;
extern char *optarg;
#endif

int 
main(int argc, char **argv) {

  char line[MAX_STR];
  char pgm_desc[MAX_STR];
  char s_name0[MAX_STR], s_name1[MAX_STR];
  char s_desc0[MAX_STR], s_desc1[MAX_STR];
  int p0_beg, p1_beg, p0_end, p1_end;
  int open_plt = 0;
  int copt;

  /* check options */
  while ((copt = getopt(argc, argv, "BZ:")) != -1) {
    switch (copt) {
    case 'B':
      have_bits = 1;
      break;
    case 'Z':
      sscanf(optarg, "%ld", &zdb_size);
      have_zdb = 1;
      break;
    case '?':
    default:
      fprintf(stderr," usage -  ps_lav -B -Z db_size\n");
    }
  }

  while (fgets(line,sizeof(line), stdin)!=NULL) {
    if (line[0] == '#' || line[0] == '\n') continue;
    else {
      switch(line[0]) {
      case 'd':
	get_str(stdin, pgm_desc, sizeof(pgm_desc));
	break;
      case 'h':
	get_str2(stdin, 
		 s_desc0, sizeof(s_desc0),
		 s_desc1, sizeof(s_desc1));
	break;
      case 's':
	get_seq_info(stdin, 
		     s_name0, sizeof(s_name0), &p0_beg, &p0_end,
		     s_name1, sizeof(s_name1), &p1_beg, &p1_end);
	g_n0 = p0_end - p0_beg + 1;
	g_n1 = p1_end - p1_beg + 1;
	break;
      case 'a':
	if (!open_plt) {
	  openplt(g_n0, g_n1, p0_beg, p1_beg,  s_desc0, s_desc1);
	  if ((g_n0 == g_n1) && (p0_beg == p1_beg) && (p0_end == p1_end) &&
	      strcmp(s_name0, s_name1) == 0) {
	    drawdiag(p0_end-p0_beg + 1, p1_end - p1_beg + 1);
	  }
	  open_plt = 1;
	}

	do_alignment(stdin, p0_beg, p1_beg);
	break;
      }
    }
  }
  if (!open_plt) {
    openplt(g_n0, g_n1, p0_beg, p1_beg,  s_desc0, s_desc1);
    if ((g_n0 == g_n1) && (p0_beg == p1_beg) && (p0_end == p1_end) &&
	strcmp(s_name0, s_name1) == 0) {
      drawdiag(p0_end-p0_beg + 1, p1_end - p1_beg + 1);
    }
    open_plt = 1;
  }
  closeplt();
  exit(0);
}

void
get_str(FILE *file, char *str, size_t len) {

  char line[MAX_STR], *bp, *bp1;
  int have_quote = 0;

  str[0] = '\0';

  while (fgets(line, sizeof(line), file)!=NULL) {
    if (line[0] == '}') return;
    if ((bp = strchr(line,'\n'))!=NULL) *bp = '\0';
    if (have_quote == 0 && (bp=strchr(line,'\"'))!=NULL) {
      have_quote = 1;
      if ((bp1 = strchr(bp+1, '\"'))!=NULL) {
	*bp1 = '\0';
	have_quote = 2;
      }
      strncat(str, bp+1, len-1);
      len -= strlen(bp+1);
    }
    else if (have_quote == 1) {
      if ((bp = strchr(line, '\"'))!=NULL) *bp = '\0';
      strncat(str, line, len-1);
      len -= strlen(line);
    }
  }
}

void
get_str2(FILE *file, 
	 char *str0, size_t len0,
	 char *str1, size_t len1) {

  char line[MAX_STR], *bp, *bp1;
  int have_quote0 = 0;
  int have_quote1 = 0;

  str0[0] = str1[0] = '\0';

  while (fgets(line, sizeof(line), file)!=NULL) {
    if (line[0] == '}') return;
    if ((bp = strchr(line,'\n'))!=NULL) *bp = '\0';
    if (have_quote0 == 0 && (bp=strchr(line,'\"'))!=NULL) {
      have_quote0 = 1;
      if ((bp1 = strchr(bp+1, '\"'))!=NULL) {
	*bp1 = '\0';
	have_quote0 = 2;
      }
      strncat(str0, bp+1, len0-1);
      len0 -= strlen(bp+1);
    }
    else if (have_quote0 == 1) {
      if ((bp = strchr(line, '\"'))!=NULL) *bp = '\0';
      strncat(str0, line, len0-1);
      len0 -= strlen(line);
    }
    else if (have_quote1 == 0 && (bp=strchr(line,'\"'))!=NULL) {
      have_quote1 = 1;
      if ((bp1 = strchr(bp+1, '\"'))!=NULL) {
	*bp1 = '\0';
	have_quote1 = 2;
      }
      strncat(str1, bp+1, len1-1);
      len1 -= strlen(bp+1);
    }
    else if (have_quote1 == 1) {
      if ((bp = strchr(line, '\"'))!=NULL) *bp = '\0';
      strncat(str1, line, len1-1);
      len1 -= strlen(line);
    }
  }
}

void
get_seq_info(FILE *file, 
	     char *str0, size_t len0, int *n0_begin, int *n0_end,
	     char *str1, size_t len1, int *n1_begin, int *n1_end) {

  char line[MAX_STR], *bp;
  int have_quote0 = 0;
  int have_quote1 = 0;

  str0[0] = str1[0] = '\0';

  fgets(line, sizeof(line), file);
  if (line[0] == '}') return;
  sscanf(line, "%s %d %d", str0, n0_begin, n0_end);
  fgets(line, sizeof(line), file);
  if (line[0] == '}') return;
  sscanf(line, "%s %d %d", str1, n1_begin, n1_end);

  if ((bp = strchr(str0+1,'\"'))!=NULL) *bp = 0;
  if ((bp = strchr(str1+1,'\"'))!=NULL) *bp = 0;
  if (str0[0] == '\"') {del1(str0);}
  if (str1[0] == '\"') {del1(str1);}

  fgets(line, sizeof(line), file);	/* get the last } */
}

void
del1(char *str) {
  char *str_p;

  str_p = str++;

  while (*str++) {*str_p++ = *str++;}

  *str_p = '\0';
}

void
do_alignment(FILE *file, int p0_beg, int p1_beg) {

  char line[MAX_STR], *bp;
  int score, s0_beg, s0_end, s1_beg, s1_end, percent;
  int have_line = 0;
  double bits;

  while (fgets(line, sizeof(line), file)!=NULL) {
    if ( strchr(line,'}') != NULL) {
      clsline();
      return;
    }
    if ((bp=strchr(line,'s')) != NULL)  sscanf(bp+1, "%d %lf", &score, &bits);
    else if ((bp=strchr(line,'b')) != NULL) sscanf(bp+1, "%d %d", &s0_beg, &s1_beg);
    else if ((bp=strchr(line, 'e')) != NULL) sscanf(bp+1, "%d %d", &s0_end, &s1_end);
    else if ((bp=strchr(line, 'l')) != NULL) {
      sscanf(bp+1, "%d %d %d %d %d",
	     &s0_beg, &s1_beg, &s0_end, &s1_end, &percent);
      if (have_line) {
	sxy_draw(s0_beg-p0_beg+1, s1_beg-p1_beg+1);
	sxy_draw(s0_end-p0_beg+1, s1_end-p1_beg+1);
      }
      else {
	opnline(score, bits);
	sxy_move(s0_beg - p0_beg + 1, s1_beg - p1_beg + 1);
	sxy_draw(s0_end - p0_beg + 1, s1_end - p1_beg + 1);
	have_line = 1;
      }
    }
  }
}

#include <math.h>

#ifndef M_LN2
#define M_LN2 0.69314718055994530942
#endif

/* produce e_val from bit score */

double
bit_to_E (double bit)
{
  double a_n0, a_n1, p_val;

  a_n0 = (double)g_n0;
  a_n1 = (double)g_n1;

  p_val = a_n0 * a_n1 / pow(2.0, bit);
  if (p_val > 0.01) p_val = 1.0 - exp(-p_val);

  return (double)zdb_size * p_val;
}

