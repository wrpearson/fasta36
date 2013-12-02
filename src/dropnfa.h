
/* $Id: dropnfa.h 795 2011-07-04 15:36:12Z wrp $ */
/* $Revision: 795 $  */

/* global definitions shared by dropnfa.c and altivec.c */

#ifndef MAXSAV
#define MAXSAV 10
#endif



struct dstruct		/* diagonal structure for saving current run */
{			
   int     score;	/* hash score of current match */
   int     start;	/* start of current match */
   int     stop;	/* end of current match */
   struct savestr *dmax;   /* location in vmax[] where best score data saved */
};

struct savestr
{
   int     score;		/* pam score with segment optimization */
   int     score0;		/* pam score of best single segment */
   int     gscore;		/* score from global match */
   int     dp;			/* diagonal of match */
   int     start;		/* start of match in lib seq */
   int     stop;		/* end of match in lib seq */
};

struct bdstr {
  int CC, DD, CP, DP;
};

struct mtp_str {		/* used to hold previous static values */
  int IP;
  int *MP[3];			/* save crossing points */
  int *FP;			/* forward dividing points */
  int *MT[3];			/* 0: rep, 1: del, 2: ins  -- was char, now int */
  int *FT;			/* was char, now int */
};

struct f_struct {
  struct dstruct *diag;
  int ndo;
  int hmask;			/* hash constants */
  int *pamh1;			/* pam based array */
  int *pamh2;			/* pam based kfact array */
  int *link, *harr;		/* hash arrays */
  int kshft;			/* shift width */
  int c_gap, opt_cut;
#ifdef TFASTA
  unsigned char *aa1x;
  int n10;
#endif
  struct bdstr *bss;
  struct mtp_str mtp;
  int bss_size;
  struct swstr *ss;
  struct swstr *f_ss, *r_ss;
  int *waa_s, *waa_a;
  int **pam2p[2];
  int max_res;
  double aa0_f[MAXSQ];
  double *kar_p;

#ifdef FA_ALTIVEC
  int vec_len;
  vecInt **vec_matrix;
  vector signed ALTIVEC_SIZE *vec_HH;
  vector signed ALTIVEC_SIZE *vec_EE;

  int vec_len2;
  vecInt2 **vec_matrix2;
  vector signed ALTIVEC_SIZE2 *vec_HH2;
  vector signed ALTIVEC_SIZE2 *vec_EE2;
#endif
};

static int
FLOCAL_ALIGN(const unsigned char *A, const unsigned char *B,
	     int M, int N, int low, int up,
	     int **W, int G,int H, int MW,
	     struct f_struct *f_str);
