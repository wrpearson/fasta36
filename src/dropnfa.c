/* $Id: dropnfa.c $ */

/* copyright (c) 1998, 1999, 2014 by William R. Pearson and The Rector &
   Vistors of the University of Virginia */

/* Licensed under the Apache License, Version 2.0 (the "License");
   you may not use this file except in compliance with the License.
   You may obtain a copy of the License at

   http://www.apache.org/licenses/LICENSE-2.0

   Unless required by applicable law or agreed to in writing,
   software distributed under this License is distributed on an "AS
   IS" BASIS, WITHOUT WRRANTIES OR CONDITIONS OF ANY KIND, either
   express or implied.  See the License for the specific language
   governing permissions and limitations under the License. 
*/

/* 17-Sept-2008 - modified for multiple non-overlapping alignments */

/* 18-Sep-2006 - removed global variables for alignment from nw_align
   and bg_align */

/* 18-Oct-2005 - converted to use a_res and aln for alignment coordinates */

/* 14-May-2003 - modified to return alignment start at 0, rather than
   1, for begin:end alignments
*/

/*
  implements the fasta algorithm, see:

  W. R. Pearson, D. J. Lipman (1988) "Improved tools for biological
  sequence comparison" Proc. Natl. Acad. Sci. USA 85:2444-2448

  This version uses Smith-Waterman for final protein alignments

  W. R. Pearson (1996) "Effective protein sequence comparison"
  Methods Enzymol. 266:227-258


  26-April-2001 - -DGAP_OPEN redefines -f, as gap open penalty

  4-Nov-2001 - modify spam() while(1).
*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <math.h>

#include "defs.h"
#include "param.h"

/* this must be consistent with upam.h */
#define MAXHASH 32
#define NMAP MAXHASH+1

/* globals for fasta */
#define MAXWINDOW 64

#ifndef MAXSAV
#define MAXSAV 10
#endif

static char *verstr="3.8 Nov 2011";

extern void w_abort(char *, char *);
int shscore(const unsigned char *aa0, int n0, int **pam2);
extern void init_karlin(const unsigned char *aa0, int n0, struct pstruct *ppst,
			double *aa0_f, double **kp);
extern void init_karlin_a(struct pstruct *, double *, double **);
extern int do_karlin(const unsigned char *, int n1, int **,
		     const struct pstruct *, double *, double *,
		     double *, double *);
extern int ELK_to_s(double E_join, int n0, int n1, double Lambda, double K, double H);

extern void aancpy(char *to, char *from, int count, struct pstruct *ppst);
char *ckalloc(size_t);
unsigned long adler32(unsigned long, const unsigned char *, unsigned int);

#ifdef TFASTA
extern int aatran(const unsigned char *ntseq, unsigned char *aaseq, int maxs, int frame);
#endif

#include "dropnfa.h"

#define DROP_INTERN
#include "drop_func.h"

struct swstr { int H, E;};

int
dmatch (const unsigned char *aa0, int n0,
	const unsigned char *aa1, int n1,
	int hoff, int window, 
	int **pam2, int gdelval, int ggapval,
	struct f_struct *f_str);


extern int sw_walign (int **pam2p, int n0,
		      const unsigned char *aa1, int n1,
		      int q, int r,
		      struct swstr *ss,
		      struct a_res_str *a_res
		      );

int bd_walign (const unsigned char *aa0, int n0,
	       const unsigned char *aa1, int n1,
	       struct pstruct *ppst, 
	       struct f_struct *f_str, int hoff,
	       struct a_res_str *a_res);

struct a_res_str *
merge_ares_chains(struct a_res_str *cur_ares, 
		  struct a_res_str *tmpl_ares,
		  int score_ix, const char *msg);

/* initialize for fasta */

void
init_work (unsigned char *aa0, int n0, 
	   struct pstruct *ppst,
	   struct f_struct **f_arg
	   )
{
  int mhv, phv;
  int hmax;
  int i0, hv;
  int pamfact;
  int btemp;
  struct f_struct *f_str;
  /* these used to be globals, but do not need to be */
  int ktup;		/* word size examined */
  int fact;		/* factor used to scale ktup match value */
  int kt1;		/* ktup-1 */
  int lkt;		/* last ktup - initiall kt1, but can be increased
			   for param_u.fa.hsq >= NMAP */

  int maxn0;		/* used in band alignment */
  int *pwaa_a, *pwaa_s;		/* pam[aa0[]] profile */
  int i, j, e, f;
  struct swstr *ss;
  int **pam2p;
  int *waa;
  int nsq, ip, *hsq;

  if (ppst->ext_sq_set) {
    nsq = ppst->nsqx; ip = 1;
    hsq = ppst->hsqx;
  }
  else {
    nsq = ppst->nsqx; ip = 0;
    hsq = ppst->hsq;
  }

  f_str = (struct f_struct *)calloc(1,sizeof(struct f_struct));

#ifndef TFASTA
  if((ppst->zsflag%10) == 6) {
    f_str->kar_p = NULL;
    init_karlin(aa0, n0, ppst, &f_str->aa0_f[0], &f_str->kar_p);
  }
#endif

  if (!ppst->param_u.fa.use_E_thresholds) {	/* old fashioned thresholds */
    btemp = 2 * ppst->param_u.fa.bestoff / 3 +
      n0 / ppst->param_u.fa.bestscale +
      ppst->param_u.fa.bkfact * 
      (ppst->param_u.fa.bktup - ppst->param_u.fa.ktup);

    if (ppst->nt_align) {
      btemp = (btemp*ppst->pam_h)/5;  /* normalize to standard +5/-4 */
    }
    else {
      /* correct problem produced by allowing bktup=3, which increases
	 btemp by bkfact for proteins 2-May-2011 */
      if (ppst->param_u.fa.ktup < 3) btemp -= ppst->param_u.fa.bkfact;
    }

    btemp = min (btemp, ppst->param_u.fa.bestmax);
    if (btemp > 3 * n0) btemp = 3 * shscore(aa0,n0,ppst->pam2[0]) / 5;

    ppst->param_u.fa.cgap = btemp + ppst->param_u.fa.bestoff / 3;

    if (ppst->param_u.fa.optcut_set != 1) {
#ifndef TFASTA
      ppst->param_u.fa.optcut = btemp;
#else
      ppst->param_u.fa.optcut = (btemp*3)/2;
#endif
    }
  }

#ifndef OLD_FASTA_GAP
  ppst->param_u.fa.pgap = ppst->gdelval + 2*ppst->ggapval;
#else
  ppst->param_u.fa.pgap = ppst->gdelval + ppst->ggapval;
#endif
  pamfact = ppst->param_u.fa.pamfact;
  ktup = ppst->param_u.fa.ktup;
  fact = ppst->param_u.fa.scfact * ktup;

  if (pamfact == -1) pamfact = 0;
  else if (pamfact == -2) pamfact = 1;

  for (i0 = 1, mhv = -1; i0 < ppst->nsq; i0++)
    if (hsq[i0] < NMAP && hsq[i0] > mhv) mhv = hsq[i0];

  if (mhv <= 0) {
    fprintf (stderr, "*** ERROR [%s:%d] maximum hsq <=0 %d\n", __FILE__, __LINE__, mhv);
    exit (1);
  }

  for (f_str->kshft = 0; mhv > 0; mhv /= 2) f_str->kshft++;

  /*      kshft = 2;	*/
  kt1 = ktup - 1;
  hv = 1;
  for (i0 = 0; i0 < ktup; i0++) hv = hv << f_str->kshft;
  hmax = hv;
  f_str->hmask = (hmax >> f_str->kshft) - 1;

  if ((f_str->harr = (int *) calloc (hmax, sizeof (int))) == NULL) {
    fprintf (stderr, "*** ERROR [%s:%d] - cannot allocate hash array: hmax: %d hmask: %d\n",
	     __FILE__,__LINE__,hmax, f_str->hmask);
    exit (1);
  }

  if ((f_str->pamh1 = (int *) calloc (nsq+1, sizeof (int))) == NULL) {
    fprintf (stderr, "*** ERROR [%s:%d] - cannot allocate pamh1 array nsq=%d\n",
	     __FILE__, __LINE__, nsq);
    exit (1);
  }

  if ((f_str->pamh2 = (int *) calloc (hmax, sizeof (int))) == NULL) {
    fprintf (stderr, "*** ERROR [%s:%d] - cannot allocate pamh2 array hmax=%d\n",
	     __FILE__, __LINE__,hmax);
    exit (1);
  }

  if ((f_str->link = (int *) calloc (n0, sizeof (int))) == NULL) {
    fprintf (stderr, "*** ERROR [%s:%d] - cannot allocate hash link array n0=%d",
	     __FILE__, __LINE__, n0);
    exit (1);
  }

  for (i0 = 0; i0 < hmax; i0++) f_str->harr[i0] = -1;
  for (i0 = 0; i0 < n0; i0++) f_str->link[i0] = -1;

  /* encode the aa0 array */
  phv = hv = 0;
  lkt = kt1;
  /* restart hv, phv calculation */
  for (i0 = 0; i0 < min(lkt,n0); i0++) {
    if (hsq[aa0[i0]] >= NMAP) {hv=phv=0; lkt = i0+ ktup; continue;}
    hv = (hv << f_str->kshft) + hsq[aa0[i0]];
    phv += ppst->pam2[ip][aa0[i0]][aa0[i0]]*ktup;
  }

  for (; i0 < n0; i0++) {
    if (hsq[aa0[i0]] >= NMAP) {
      hv=phv=0;
      /* restart hv, phv calculation */
      for (lkt = i0+kt1; (i0 < lkt || hsq[aa0[i0]]>=NMAP) && i0<n0; i0++) {
	if (hsq[aa0[i0]] >= NMAP) {
	  hv=phv=0; 
	  lkt = i0+ktup;
	  continue;
	}
	hv = (hv << f_str->kshft) + hsq[aa0[i0]];
	phv += ppst->pam2[ip][aa0[i0]][aa0[i0]]*ktup;
      }
    }
    if (i0 >= n0) break;
    hv = ((hv & f_str->hmask) << f_str->kshft) + hsq[aa0[i0]];
    f_str->link[i0] = f_str->harr[hv];
    f_str->harr[hv] = i0;
    if (pamfact) {
      f_str->pamh2[hv] = (phv += ppst->pam2[ip][aa0[i0]][aa0[i0]] * ktup);
      /* this check should always be true, but just in case */
      if (hsq[aa0[i0-kt1]]<NMAP)
	phv -= ppst->pam2[ip][aa0[i0 - kt1]][aa0[i0 - kt1]] * ktup;
    }
    else f_str->pamh2[hv] = fact * ktup;
  }

  /* this has been modified from 0..<ppst->nsq to 1..< ppst->nsq because the
     pam2[0][0] is now undefined for consistency with blast
  */

  if (pamfact) {
    for (i0 = 1; i0 < nsq; i0++) {
      f_str->pamh1[i0] = ppst->pam2[ip][i0][i0] * ktup;
    }
  }
  else {
    for (i0 = 1; i0 < nsq; i0++) {
      f_str->pamh1[i0] = fact;
    }
  }
  f_str->ndo = 0;
  if ((f_str->diag = (struct dstruct *) calloc ((size_t)MAXDIAG,
						sizeof (struct dstruct)))==NULL) {
    fprintf (stderr,"*** ERROR [%s:%d] - cannot allocate diagonal arrays: %lu\n",
	     __FILE__, __LINE__, MAXDIAG *sizeof (struct dstruct));
    exit (1);
  };


#ifdef TFASTA
  if ((f_str->aa1x =(unsigned char *)calloc((size_t)ppst->maxlen+2,
					    sizeof(unsigned char)))
      == NULL) {
    fprintf (stderr, "*** ERROR [%s:%d] - cannot allocate aa1x array %d\n",
	     __FILE__, __LINE__, ppst->maxlen+2);
    exit (1);
  }
  f_str->aa1x++;
#endif

  f_str->bss_size = ppst->param_u.fa.optwid*2+4;
  f_str->bss = (struct bdstr *) calloc((size_t)ppst->param_u.fa.optwid*2+4,
				       sizeof(struct bdstr));
  f_str->bss++;

  /* allocate space for the scoring arrays */
  maxn0 = n0 + 4;
  if ((ss = (struct swstr *) calloc (maxn0, sizeof (struct swstr)))
      == NULL) {
    fprintf (stderr, "*** ERROR [%s:%d] - cannot allocate ss array %3d\n",
	     __FILE__, __LINE__, n0);
    exit (1);
  }
  ss++;

  ss[n0].H = -1;	/* this is used as a sentinel - normally H >= 0 */
  ss[n0].E = 1;
  f_str->ss = ss;

  /* initialize variable (-S) pam matrix */
  if ((f_str->waa_s= (int *)calloc((nsq+1)*(n0+1),sizeof(int))) == NULL) {
    fprintf(stderr,"*** ERROR [%s:%d] - cannot allocate waa_s array %3d\n",
	    __FILE__, __LINE__, nsq*n0);
    exit(1);
  }

  /* initialize pam2p[1] pointers */
  if ((f_str->pam2p[1]= (int **)calloc((n0+1),sizeof(int *))) == NULL) {
    fprintf(stderr,"*** ERROR [%s:%d] - cannot allocate pam2p[1] array %3d\n",
	    __FILE__, __LINE__, n0);
    exit(1);
  }

  pam2p = f_str->pam2p[1];
  if ((pam2p[0]=(int *)calloc((nsq+1)*(n0+1),sizeof(int))) == NULL) {
    fprintf(stderr,"*** ERROR [%s:%d] - cannot allocate pam2p[1][] array %3d\n",
	    __FILE__, __LINE__, nsq*n0);
    exit(1);
  }

  for (i=1; i<n0; i++) {
    pam2p[i]= pam2p[0] + (i*(nsq+1));
  }

  /* initialize universal (alignment) matrix */
  if ((f_str->waa_a= (int *)calloc((nsq+1)*(n0+1),sizeof(int))) == NULL) {
    fprintf(stderr,"*** ERROR [%s:%d] - cannot allocate waa_a struct %3d\n",
	    __FILE__, __LINE__, nsq*n0);
    exit(1);
  }
   
  /* initialize pam2p[0] pointers */
  if ((f_str->pam2p[0]= (int **)calloc((n0+1),sizeof(int *))) == NULL) {
    fprintf(stderr,"*** ERROR [%s:%d] - cannot allocate pam2p[1] array %3d\n",
	    __FILE__, __LINE__, n0);
    exit(1);
  }

  pam2p = f_str->pam2p[0];
  if ((pam2p[0]=(int *)calloc((nsq+1)*(n0+1),sizeof(int))) == NULL) {
    fprintf(stderr,"*** ERROR [%s:%d] - cannot allocate pam2p[1][] array %3d\n",
	    __FILE__, __LINE__, nsq*n0);
    exit(1);
  }

  for (i=1; i<n0; i++) {
    pam2p[i]= pam2p[0] + (i*(nsq+1));
  }

  /* 
     pwaa effectively has a sequence profile --
     pwaa[0..n0-1] has pam score for residue 0 (-BIGNUM)
     pwaa[n0..2n0-1] has pam scores for residue 1 (A)
     pwaa[2n0..3n-1] has pam scores for residue 2 (R), ...

     thus: pwaa = f_str->waa_s + (*aa1p++)*n0; sets up pwaa so that
     *pwaa++ rapidly moves though the scores of the aa1p[] position
     without further indexing

     For a real sequence profile, pwaa[0..n0-1] vs ['A'] could have
     a different score in each position.
  */

  if (ppst->pam_pssm) {
    pwaa_s = f_str->waa_s;
    pwaa_a = f_str->waa_a;
    for (e = 0; e < nsq; e++)	{	/* for each residue in the alphabet */
      for (f = 0; f < n0; f++) {	/* for each position in aa0 */
	*pwaa_s++ = f_str->pam2p[ip][f][e] = ppst->pam2p[ip][f][e];
	*pwaa_a++ = f_str->pam2p[0][f][e]  = ppst->pam2p[0][f][e];
      }
    }
  }
  else {	/* initialize scanning matrix */
    pwaa_s = f_str->waa_s;
    pwaa_a = f_str->waa_a;
    for (e = 0; e <nsq; e++)	/* for each residue in the alphabet */
      for (f = 0; f < n0; f++)	{	/* for each position in aa0 */
	*pwaa_s++ = f_str->pam2p[ip][f][e]= ppst->pam2[ip][aa0[f]][e];
	*pwaa_a++ = f_str->pam2p[0][f][e] = ppst->pam2[0][aa0[f]][e];
      }
  }

  f_str->max_res = max(3*n0/2,MIN_RES);

  *f_arg = f_str;
}

/* pstring1 is a message to the manager, currently 512 */
/* pstring2 is the same information, but in a markx==10 format */
void
get_param (const struct pstruct *ppstr, 
	   char **pstring1, char *pstring2, struct score_count_s *s_info)
{
  char options_str1[128];
  char options_str2[128];
#ifndef TFASTA
  char *pg_str="FASTA";
#else
  char *pg_str="TFASTA";
#endif

  if (!ppstr->param_u.fa.use_E_thresholds) {
    sprintf(options_str1,"join: %d (%0.3g), opt: %d (%0.3g)",
	    ppstr->param_u.fa.cgap, (double)s_info->s_cnt[0]/(double)s_info->tot_scores,
	    ppstr->param_u.fa.optcut, (double)s_info->s_cnt[2]/(double)s_info->tot_scores);
    sprintf(options_str2,"pg_join: %d (%.3g)\n; pg_optcut: %d (%.3g)",
	    ppstr->param_u.fa.cgap, (double)s_info->s_cnt[0]/(double)s_info->tot_scores,
	    ppstr->param_u.fa.optcut, (double)s_info->s_cnt[2]/(double)s_info->tot_scores);
  }
  else {
    sprintf(options_str1,"E-join: %.2g (%.3g), E-opt: %.2g (%.3g)",
	    ppstr->param_u.fa.E_join, (double)s_info->s_cnt[0]/(double)s_info->tot_scores,
	    ppstr->param_u.fa.E_band_opt, (double)s_info->s_cnt[2]/(double)s_info->tot_scores);
    sprintf(options_str2,"pg_join_E(): %.2g (%.3g)\n; pg_optcut_E(): %.2g (%.3g)",
	    ppstr->param_u.fa.E_join, (double)s_info->s_cnt[0]/(double)s_info->tot_scores,
	    ppstr->param_u.fa.E_band_opt, (double)s_info->s_cnt[2]/(double)s_info->tot_scores);
  }

  if (!ppstr->param_u.fa.optflag) {
     sprintf (pstring1[0], "%s (%s)", pg_str, verstr);
     if (ppstr->param_u.fa.iniflag) strcat(pstring1[0]," init1");
  }
  else {
     sprintf (pstring1[0], "%s (%s) [optimized]", pg_str, verstr);
  }

  sprintf (pstring1[1], 
#ifdef OLD_FASTA_GAP
	   "%s matrix (%d:%d)%s, gap-pen: %d/%d\n ktup: %d, %s, width: %3d",
#else
	   "%s matrix (%d:%d)%s, open/ext: %d/%d\n ktup: %d, %s, width: %3d",
#endif
	   ppstr->pam_name, ppstr->pam_h,ppstr->pam_l,
	   (ppstr->ext_sq_set) ? "xS":"\0",
	   ppstr->gdelval, ppstr->ggapval,
	   ppstr->param_u.fa.ktup, options_str1,
	   ppstr->param_u.fa.optwid);


   if (pstring2 != NULL) {
     sprintf (pstring2, 
#ifdef OLD_FASTA_GAP
	      "; pg_name_alg: %s\n; pg_ver_rel: %s\n; pg_matrix: %s (%d:%d)\n\
; pg_gap-pen: %d %d\n; pg_ktup: %d\n; %s\n",
#else
     "; pg_name_alg: %s\n; pg_ver_rel: %s\n; pg_matrix: %s (%d:%d)\n\
; pg_open-ext: %d %d\n; pg_ktup: %d\n; %s\n",
#endif
	      pg_str,verstr,ppstr->pam_name, ppstr->pam_h,ppstr->pam_l, ppstr->gdelval,
              ppstr->ggapval,ppstr->param_u.fa.ktup,options_str2);
   }
}

void
close_work (const unsigned char *aa0, int n0,
	    struct pstruct *ppst,
	    struct f_struct **f_arg)
{
  struct f_struct *f_str;

  f_str = *f_arg;

  if (f_str != NULL) {
    if (f_str->kar_p!=NULL) free(f_str->kar_p);
    f_str->ss--;
    free(f_str->ss);

    f_str->bss--;
    free(f_str->bss);
    f_str->bss_size = 0;

    /*    free(f_str->res);  */
    free(f_str->waa_a);
    free(f_str->waa_s);

    free(f_str->pam2p[1][0]);
    free(f_str->pam2p[1]);

    free(f_str->pam2p[0][0]);
    free(f_str->pam2p[0]);


    free(f_str->diag);
    free(f_str->link);
    free(f_str->pamh2); 
    free(f_str->pamh1);
    free(f_str->harr);

    free(f_str);
    *f_arg = NULL;
  }
  else {
    fprintf(stderr, "*** ERROR [%s:%d] - close_work() with NULL f_str ***\n",
	    __FILE__, __LINE__);
  }
}

int savemax (struct dstruct *, int, 
	     struct savestr *vmax, 
	     struct savestr **lowmax);


int spam (const unsigned char *, const unsigned char *, struct savestr *,
	 int **, int, int, int);
int sconn(struct savestr **, int nsave, int cgap, int pgap, int noff);
void kpsort(struct savestr **, int);

extern int
NW_ALIGN(const unsigned char *, const unsigned char *, int, int, 
      int **pam2p, int, int q, int r, int *res, int *nc);

static int
LOCAL_ALIGN(const unsigned char *, const unsigned char *,
	    int, int, int, int,
	    int **, int, int, int *, int *, int *, int *, int,
	    struct f_struct *);

static int
B_ALIGN(const unsigned char *A, const unsigned char *B, int M,
	int N, int low, int up, int **W, int G, int H, int *S,
	int *nS, int MW, int MX, struct bdstr *bss, struct mtp_str *mtp);

static void
do_fasta (const unsigned char *aa0, int n0,
	  const unsigned char *aa1, int n1,
	  const struct pstruct *ppst, struct f_struct *f_str,
	  struct rstruct *rst, int *hoff, int shuff_flg,
	  struct score_count_s *s_info) {
   int     nd;		/* diagonal array size */
   int     lhval;
   int     kfact;
   struct savestr vmax[MAXSAV];	/* best matches saved for one sequence */
   struct savestr *vptr[MAXSAV];
   struct savestr *lowmax;
   int lowscor;
   register struct dstruct *dptr;
   register int tscor;

   register struct dstruct *diagp;
   int noff;
   struct dstruct *dpmax, *dpmin;
   register int lpos;
   int     tpos;
   struct savestr *vmptr;
   int     scor, ib, nsave;
   int xdrop, do_extend;
   int ktup, kt1, lkt, ip, ktup_sq;
   const int *hsq;
   int opt_cut, c_gap;


  if (ppst->ext_sq_set) {
    ip = 1;
    hsq = ppst->hsqx;
  }
  else {
    ip = 0;
    hsq = ppst->hsq;
  }

  xdrop = -ppst->pam_l;
  /* do extended alignment in spam iff protein or short sequences */
  do_extend = !ppst->nt_align || (n0 < 50) || (n1 < 50);

  ktup = ppst->param_u.fa.ktup;
  kt1 = ktup-1;
  if (ktup <= 3) {
    ktup_sq = ktup*ktup;
  }
  else {
    ktup_sq = ktup;
  }
  if (ktup == 1) ktup_sq *= 2;

  if (n1 < ktup) {
     rst->score[0] = rst->score[1] = rst->score[2] = 0;
     return;
   }

   if (n0+n1+1 >= MAXDIAG) {
     fprintf(stderr,"*** ERROR [%s:%d] - n0,n1 too large: %d + %d (%d) > %d \n",
	     __FILE__, __LINE__, n0,n1,n0+n1+1,MAXDIAG);
     rst->score[0] = rst->score[1] = rst->score[2] = -1;
     return;
   }

   /* dynamically set optcut and cgap */
   if (ppst->param_u.fa.use_E_thresholds) {
     c_gap = ELK_to_s(ppst->param_u.fa.E_join*ktup_sq, n0, n1, ppst->pLambda, ppst->pK, ppst->pH);
     opt_cut = ELK_to_s(ppst->param_u.fa.E_band_opt*ktup_sq, n0, n1, ppst->pLambda, ppst->pK, ppst->pH);
     rst->valid_stat = 0;
   }
   else {
     c_gap = ppst->param_u.fa.cgap;
     opt_cut = ppst->param_u.fa.optcut;
     rst->valid_stat = 1;
   }

   s_info->tot_scores++;

   nd = n0 + n1;

   /* here we are cleaning up after the previous run, when the
      structure is initialized (init_work), ndo is set to 0, so the
      initialization keeps going out to nd only when nd increases

      dpmax now remembers the largest diagp with hits, and only cleansa
      up to there
   */

   dpmin = dpmax = &f_str->diag[nd];
   for (dptr = &f_str->diag[f_str->ndo]; dptr < dpmax;){
      dptr->stop = -1;
      dptr->dmax = NULL;
      dptr++->score = 0;
   }
   dpmax = f_str->diag;

   for (vmptr = vmax; vmptr < &vmax[MAXSAV]; vmptr++) {
     vmptr->score = 0;
   }

   lowmax = vmax;
   lowscor = 0;

   /* start hashing */
   lhval = 0;
   lkt = kt1;
   for (lpos = 0; (lpos < lkt || hsq[aa1[lpos]]>=NMAP) && lpos <n1; lpos++) {
     /* restart lhval calculation */
     if (hsq[aa1[lpos]]>=NMAP) {
       lhval = 0; lkt = lpos + ktup;
       continue;
     }
     lhval = ((lhval & f_str->hmask) << f_str->kshft) + hsq[aa1[lpos]];
   }

   noff = n0 - 1;
   diagp = &f_str->diag[noff + lkt];
   for (; lpos < n1; lpos++, diagp++) {
     /* skip over low complexity */
     if (hsq[aa1[lpos]]>=NMAP) {
       lpos++ ; diagp++;
       while (hsq[aa1[lpos]]>=NMAP && lpos < n1 ) {lpos++; diagp++;}
       if (lpos >= n1) break;
       lhval = 0;
     }

     /* lhval is the hash value of the library sequence */
     lhval = ((lhval & f_str->hmask) << f_str->kshft) + hsq[aa1[lpos]];
     /* tpos gives the locations where the library word matches */

     /* get the diagonal of the initial hit */
     tpos = f_str->harr[lhval];
     if (tpos >= 0) {       /* we have a hit */
       /* tpos = link[tpos] means that tpos always gets smaller, do
	  diagp-tpos gets bigger */
       if (diagp-tpos < dpmin) dpmin = diagp - tpos;
       for (; tpos >= 0; tpos = f_str->link[tpos]) {
	 /* here tscor is actually the end of the run */
	 if ((tscor = (dptr = &diagp[-tpos])->stop) >= 0) {
	   /* increased by ktup word length */
	   tscor += ktup;
	   /* now tscor becomes the penalty for the unmatched
	      (non-identical) residues in the diagonal between the
	      previous match and the current match (lpos) */
	   if ((tscor -= lpos) <= 0) {
	     scor = dptr->score;
	     /* the score is getting worse; if it is better than lowscor, save it */
	     if ((tscor += (kfact = f_str->pamh2[lhval])) < 0 && lowscor < scor) {
	       lowscor = savemax(dptr, dptr - f_str->diag, vmax, &lowmax);
	     }
	     /* tscor is now a candidate score of the current run;
		kfact is the score for starting over */
	     if ((tscor += scor) >= kfact) {
	       /* continuing better */
	       dptr->score = tscor;
	       dptr->stop = lpos;
	     }
	     else {
	       /* starting over better */
	       dptr->score = kfact;
	       dptr->start = (dptr->stop = lpos) - kt1;
	     }
	   }
	   else {	 /* continuing a match */
	     dptr->score += f_str->pamh1[aa0[tpos]];
	     dptr->stop = lpos;
	   }
	 }
	 else {       /* no run in this diagonal yet */
	   dptr->score = f_str->pamh2[lhval];
	   dptr->start = (dptr->stop = lpos) - kt1;
	 }
	 /* dptr is biggest at the end of tpos = link[tpos] */
	 if (dptr > dpmax) dpmax = dptr;
       }
     }				/* end tpos */
   }				/* end lpos */

   for (dptr = dpmin; dptr <= dpmax;) {
     if (dptr->score > lowscor) {
      lowscor = savemax (dptr, dptr-f_str->diag, vmax, &lowmax);
     }
     dptr->stop = -1;
     dptr->dmax = NULL;
     dptr++->score = 0;
   }
   f_str->ndo = nd;

/*
        at this point all of the elements of aa1[lpos]
        have been searched for elements of aa0[tpos]
        with the results in diag[dpos]
*/
   for (nsave = 0, vmptr = vmax; vmptr < &vmax[MAXSAV]; vmptr++)    {
     /*
     fprintf(stderr,"0: %4d-%4d  1: %4d-%4d  dp: %d score: %d\n",
	     noff+vmptr->start-vmptr->dp,
	     noff+vmptr->stop-vmptr->dp,
	     vmptr->start,vmptr->stop,
	     vmptr->dp,vmptr->score);

     */
      if (vmptr->score > 0) {
	vmptr->score = spam (aa0, aa1, vmptr, ppst->pam2[ip], xdrop,
			     noff,do_extend);
	if (vmptr->score > 0 ) vptr[nsave++] = vmptr;
      }
   }

   if (nsave <= 0) {
     rst->score[0] = rst->score[1] = rst->score[2] = 0;
     return;
   }
       
   /*
   fprintf(stderr,"n0: %d; n1: %d; noff: %d\n",n0,n1,noff);
   for (ib=0; ib<nsave; ib++) {
     fprintf(stderr,"0: %4d-%4d  1: %4d-%4d  dp: %d score: %d\n",
	     noff+vptr[ib]->start-vptr[ib]->dp,
	     noff+vptr[ib]->stop-vptr[ib]->dp,
	     vptr[ib]->start,vptr[ib]->stop,
	     vptr[ib]->dp,vptr[ib]->score);
   }
   fprintf(stderr,"---\n");
   */

   /* find the best init1 score */
   for (vmptr=vptr[0],ib=1; ib<nsave; ib++) {
     if (vptr[ib]->score > vmptr->score) vmptr=vptr[ib];
   }

   /* sconn does not modify vmptr->score, so only do it if it will help */
   if (vmptr->score >= c_gap ) {
     s_info->s_cnt[0]++;
     scor = sconn (vptr, nsave, c_gap, ppst->param_u.fa.pgap, noff);
     rst->alg_info |= 1;
   }
   else { scor = vmptr->score;}

/*  kssort (vptr, nsave); */

   rst->score[1] = vmptr->score;
   rst->score[0] = max (scor, vmptr->score);
   rst->score[2] = rst->score[0];		/* initn */

   *hoff= noff - vmptr->dp;	/* always need *hoff */
   if (ppst->param_u.fa.optflag) {
     if ( /* shuff_flg || */ rst->score[0] > opt_cut ) {
       rst->score[2] = dmatch (aa0, n0, aa1, n1, *hoff,
			     ppst->param_u.fa.optwid, ppst->pam2[ip],
			     ppst->gdelval,ppst->ggapval,f_str);

       s_info->s_cnt[2]++;
       rst->alg_info |= 2;
       rst->valid_stat = 1;
     }
     /*
     else if (rst->score[0] > c_gap) {
       rst->score[2] = dmatch (aa0, n0, aa1, n1, *hoff, 4, ppst->pam2[ip],
			     ppst->gdelval,ppst->ggapval,f_str);
     }
     */
   }
   else {	/* we never do dmatch, so initn is valid */
     rst->valid_stat = 1;
   }
}

void do_work (const unsigned char *aa0, int n0,
	      const unsigned char *aa1, int n1,
	      int frame,
	      const struct pstruct *ppst, struct f_struct *f_str,
	      int qr_flg, int shuff_flg, struct rstruct *rst, 
	      struct score_count_s *s_info)
{
  int hoff, n10;

  double lambda, H;
  
  rst->score[0] = rst->score[1] = rst->score[2] = 0;
  rst->escore = 1.0;
  rst->segnum = rst->seglen = 1;
  rst->valid_stat = 0;
  rst->alg_info = 0;

  if (n1 < ppst->param_u.fa.ktup) return;

#ifdef TFASTA  
  n10=aatran(aa1,f_str->aa1x,n1,frame);
  do_fasta (aa0, n0, f_str->aa1x, n10, ppst, f_str, rst, &hoff, shuff_flg, s_info);
#else	/* FASTA */
  do_fasta (aa0, n0, aa1, n1, ppst, f_str, rst, &hoff, shuff_flg, s_info);
#endif

#ifndef TFASTA
  if((ppst->zsflag%10) == 6 && 
     (do_karlin(aa1, n1, ppst->pam2[0], ppst,f_str->aa0_f, 
		f_str->kar_p, &lambda, &H)>0)) {
    rst->comp = 1.0/lambda;
    rst->H = H;
  }
  else {rst->comp = rst->H = -1.0;}
#else
  rst->comp = rst->H = -1.0;
#endif
}

void do_opt (const unsigned char *aa0, int n0,
	     const unsigned char *aa1, int n1,
	     int frame,
	     struct pstruct *ppst,
	     struct f_struct *f_str,
	     struct rstruct *rst)
{
  int optflag, tscore, hoff, n10;
  struct score_count_s s_info;

  optflag = ppst->param_u.fa.optflag;
  ppst->param_u.fa.optflag = 1;

#ifdef TFASTA  
  n10=aatran(aa1,f_str->aa1x,n1,frame);
  do_fasta (aa0, n0, f_str->aa1x, n10, ppst, f_str, rst, &hoff, 0, &s_info);
#else	/* FASTA */
  do_fasta(aa0,n0,aa1,n1,ppst,f_str,rst, &hoff, 0, &s_info);
#endif
  ppst->param_u.fa.optflag = optflag;
}

int
savemax (struct dstruct *dptr, int dpos,
	 struct savestr *vmax, struct savestr **lowmax)
{
   struct savestr *vmptr;
   int i;

/* check to see if this is the continuation of a run that is already saved */

   if ((vmptr = dptr->dmax) != NULL && vmptr->dp == dpos &&
	 vmptr->start == dptr->start) {
      vmptr->stop = dptr->stop;
      if ((i = dptr->score) <= vmptr->score) return (*lowmax)->score;
      vmptr->score = i;
      if (vmptr != (*lowmax)) return (*lowmax)->score;
   }
   else {
     i = (*lowmax)->score = dptr->score;
     (*lowmax)->dp = dpos;
     (*lowmax)->start = dptr->start;
     (*lowmax)->stop = dptr->stop;
     dptr->dmax = (*lowmax);
   }

   for (vmptr = vmax; vmptr < vmax+MAXSAV; vmptr++) {
     if (vmptr->score < i) {
       i = vmptr->score;
       *lowmax = vmptr;
     }
   }
   return i;
}

int spam (const unsigned char *aa0, const unsigned char *aa1,
	  struct savestr *dmax, int **pam2, int xdrop,
	  int noff, int do_extend)
{
   register int lpos, tot;
   register const unsigned char *aa0p, *aa1p;

   int drop_thresh;

   struct {
     int     start, stop, score;
   } curv, maxv;

   aa1p = &aa1[lpos= dmax->start];	/* get the start of lib seq */
   aa0p = &aa0[lpos - dmax->dp + noff];	/* start of query */
#ifdef DEBUG
   /* also add check in calling routine */
   if (aa0p < aa0) { return -99; }
#endif
   curv.start = lpos;			/* start index in lib seq */

   tot = curv.score = maxv.score = 0;

   for (; lpos <= dmax->stop; lpos++) {
     tot += pam2[*aa0p++][*aa1p++];
     if (tot > curv.score) {		/* update current score */
       curv.stop = lpos;
       curv.score = tot;
      }
      else if (tot < 0) {
	if (curv.score > maxv.score) {	/* save score, start, stop */
	  maxv.start = curv.start;
	  maxv.stop = curv.stop;
	  maxv.score = curv.score;
	}
	tot = curv.score = 0;		/* reset running score */
	curv.start = lpos+1;		/* reset start */
	if(lpos >= dmax->stop) break;	/* if the zero is beyond stop, quit */
      }
   }

   if (curv.score > maxv.score) {
     maxv.start = curv.start;
     maxv.stop = curv.stop;
     maxv.score = curv.score;
   }

#ifndef NOSPAM_EXT

   /* now check to see if the score gets better by extending */
   if (do_extend && maxv.score > xdrop) {

     if (maxv.stop == dmax->stop) {
       tot = maxv.score;
       drop_thresh = maxv.score - xdrop;
       aa1p = &aa1[lpos= dmax->stop];
       aa0p = &aa0[lpos - dmax->dp + noff];
       while (tot > drop_thresh ) {
	 ++lpos;
	 tot += pam2[*(++aa0p)][*(++aa1p)];
	 if (tot > maxv.score) {
	   maxv.start = lpos;
	   maxv.score = tot;
	   drop_thresh = tot - xdrop;
	 }
       }
     }

     /* scan backwards now */

     if (maxv.start == dmax->start) {
       tot = maxv.score;
       drop_thresh = maxv.score - xdrop;
       aa1p = &aa1[lpos= dmax->start];
       aa0p = &aa0[lpos - dmax->dp + noff];
       while (tot > drop_thresh) {
	 --lpos;
	 tot += pam2[*(--aa0p)][*(--aa1p)];
	 if (tot > maxv.score) {
	   maxv.start = lpos;
	   maxv.score = tot;
	   drop_thresh = tot - xdrop;
	 }
       }
     }
   }
#endif

/*	if (maxv.start != dmax->start || maxv.stop != dmax->stop)
		printf(" new region: %3d %3d %3d %3d\n",maxv.start,
			dmax->start,maxv.stop,dmax->stop);
*/
   dmax->start = maxv.start;
   dmax->stop = maxv.stop;

   return maxv.score;
}

int sconn (struct savestr **v, int n, int cgap, int pgap, int noff)
{
   int     i, si;
   struct slink
   {
      int     score;
      struct savestr *vp;
      struct slink *next;
   }      *start, *sl, *sj, *so, sarr[MAXSAV];
   int     lstart, tstart, plstop, ptstop;

/*	sort the score left to right in lib pos */

   kpsort (v, n);

   start = NULL;

/*	for the remaining runs, see if they fit */

   for (i = 0, si = 0; i < n; i++) {

/*	if the score is less than the gap penalty, it never helps */
     if (v[i]->score < cgap) continue;
     lstart = v[i]->start;
     tstart = lstart - v[i]->dp + noff;

/*	put the run in the group */
     sarr[si].vp = v[i];
     sarr[si].score = v[i]->score;
     sarr[si].next = NULL;

/* 	if it fits, then increase the score */
      for (sl = start; sl != NULL; sl = sl->next)
      {
	 plstop = sl->vp->stop;
	 ptstop = plstop - sl->vp->dp + noff;
	 if (plstop < lstart && ptstop < tstart)
	 {
	    sarr[si].score = sl->score + v[i]->score + pgap;
	    break;
	 }
      }

/*	now recalculate where the score fits */
      if (start == NULL)
	 start = &sarr[si];
      else
	 for (sj = start, so = NULL; sj != NULL; sj = sj->next) {
	    if (sarr[si].score > sj->score) {
	       sarr[si].next = sj;
	       if (so != NULL) so->next = &sarr[si];
	       else  start = &sarr[si];
	       break;
	    }
	    so = sj;
	 }
      si++;
   }

   if (start != NULL)
      return (start->score);
   else
      return (0);
}

void
kssort (v, n)
struct savestr *v[];
int     n;
{
   int     gap, i, j;
   struct savestr *tmp;

   for (gap = n / 2; gap > 0; gap /= 2)
      for (i = gap; i < n; i++)
	 for (j = i - gap; j >= 0; j -= gap)
	 {
	    if (v[j]->score >= v[j + gap]->score)
	       break;
	    tmp = v[j];
	    v[j] = v[j + gap];
	    v[j + gap] = tmp;
	 }
}

void
kpsort (struct savestr **v, int n) {
  int gap, i, j, k;
  int incs[4] = { 21, 7, 3, 1 };
  struct savestr *tmp;
  int v_start;

  for ( k = 0; k < 4; k++) {
    gap = incs[k];
    for (i = gap; i < n; i++) {
      tmp = v[i];
      j = i;
      v_start = v[i]->start;
      while (j >= gap && v[j - gap]->start > v_start) {
	v[j] = v[j - gap];
	j -= gap;
      }
      v[j] = tmp;
    }
  }
}

int dmatch (const unsigned char *aa0, int n0,
		   const unsigned char *aa1, int n1,
		   int hoff, int window, 
		   int **pam2, int gdelval, int ggapval,
		   struct f_struct *f_str)
{
   int low, up;

   window = min (n1, window);

#ifdef DEBUG
   if (window > f_str->bss_size) {
     fprintf(stderr,"*** ERROR [%s:%d] - dmatch window [%d] out of range [%d]\n",
	     __FILE__, __LINE__, window, f_str->bss_size);
     window = f_str->bss_size - 4;
   }
#endif

   /* hoff is the offset found from aa1 to seq 2 by dmatch */

   low = -window/2-hoff;
   up = low+window;

   return FLOCAL_ALIGN(aa0-1,aa1-1,n0,n1, low, up,
		      pam2,
#ifdef OLD_FASTA_GAP
		       -(gdelval-ggapval),
#else
		       -gdelval,
#endif
		       -ggapval,window,f_str);
 }


/* A PACKAGE FOR LOCALLY ALIGNING TWO SEQUENCES WITHIN A BAND:

   To invoke, call LOCAL_ALIGN(A,B,M,N,L,U,W,G,H,MW).
   The parameters are explained as follows:
	A, B : two sequences to be aligned
	M : the length of sequence A
	N : the length of sequence B
	L : lower bound of the band
	U : upper bound of the band
	W : scoring table for matches and mismatches
	G : gap-opening penalty
	H : gap-extension penalty
	MW  : maximum window size
*/

#include <stdio.h>

#define MININT -9999999

int
FLOCAL_ALIGN(const unsigned char *A, const unsigned char *B,
	     int M, int N, int low, int up,
	     int **W, int G,int H, int MW,
	     struct f_struct *f_str)
{
  int band;
  register struct bdstr *bssp;
  int i, j, si, ei;
  int c, d, e, m;
  int leftd, rightd;
  int best_score;
  int *wa, curd;
  int ib;
  
  bssp = f_str->bss;

  m = G+H;
  low = max(-M, low);
  up = min(N, up);
  
  if (N <= 0) return 0;

  if (M <= 0) return 0;

  band = up-low+1;
  if (band < 1) {
    fprintf(stderr,"***  error [%s:%d] - low > up is unacceptable!: M: %d N: %d l/u: %d/%d\n",
	    __FILE__, __LINE__, M, N, low, up);
    return 0;
  }

  if (low > 0) leftd = 1;
  else if (up < 0) leftd = band;
  else leftd = 1-low;
  rightd = band;
  si = max(0,-up);	/* start index -1 */
  ei = min(M,N-low);	/* end index */
  bssp[leftd].CC = 0;
  for (j = leftd+1; j <= rightd; j++) {
    bssp[j].CC = 0;
    bssp[j].DD = -G;
  }

  bssp[rightd+1].CC = MININT;
  bssp[rightd+1].DD = MININT;

  best_score = 0;
  bssp[leftd-1].CC = MININT;
  bssp[leftd].DD = -G;

  for (i = si+1; i <= ei; i++) {
    if (i > N-up) rightd--;
    if (leftd > 1) leftd--;
    wa = W[A[i]];
    if ((c = bssp[leftd+1].CC-m) > (d = bssp[leftd+1].DD-H)) d = c;
    if ((ib = leftd+low-1+i ) > 0) c = bssp[leftd].CC+wa[B[ib]];

    if (d > c) c = d;
    if (c < 0) c = 0;
    e = c-G;
    bssp[leftd].DD = d;
    bssp[leftd].CC = c;
    if (c > best_score) best_score = c;

    for (curd=leftd+1; curd <= rightd; curd++) {
      if ((c = c-m) > (e = e-H)) e = c;
      if ((c = bssp[curd+1].CC-m) > (d = bssp[curd+1].DD-H)) d = c;
      c = bssp[curd].CC + wa[B[curd+low-1+i]];
      if (e > c) c = e;
      if (d > c) c = d;
      if (c < 0) c = 0;
      bssp[curd].CC = c;
      bssp[curd].DD = d;
      if (c > best_score) best_score = c;
    }
  }

  return best_score;
}

/* ckalloc - allocate space; check for success */
char *ckalloc(size_t amount)
{
  char *p;

  if ((p = malloc( (unsigned) amount)) == NULL)
    w_abort("Ran out of memory.","");
  return(p);
}

/* calculate the 100% identical score */
int
shscore(const unsigned char *aa0, int n0, int **pam2)
{
  int i, sum;
  for (i=0,sum=0; i<n0; i++)
    sum += pam2[aa0[i]][aa0[i]];
  return sum;
}

static int
BCHECK_SCORE(const unsigned char *A, const unsigned char *B,
		       int M, int N, int *S, int **w, int g, int h,
		       int *nres)
{ 
  register int   i,  j, op, nc;
  int *Ssave;
  int score;

  score = i = j = op = nc = 0;
  Ssave = S;
  while (i < M || j < N) {
	op = *S++;
	if (op == 0) {
	  score = w[A[++i]][B[++j]] + score;
	  nc++;
/*  	  fprintf(stderr,"op0 %4d %4d %4d %4d\n",i,j,w[A[i]][B[i]],score); */
	}
	else if (op > 0) {
	  score = score - (g+op*h);
/*  	  fprintf(stderr,"op> %4d %4d %4d %4d %4d\n",i,j,op,-(g+op*h),score); */
	  j += op;
	  nc += op;
	} else {
	  score = score - (g-op*h);
/*  	  fprintf(stderr,"op< %4d %4d %4d %4d %4d\n",i,j,op,-(g-op*h),score); */
	  i -= op;
	  nc -= op;
	}
  }
  *nres = nc;
  return score;
}

/* bd_malign is a recursive interface to bd_walign() that is called
   from do_walign(). bd_malign() first does an alignment, then checks
   to see if the score is greater than the threshold. If so, it tries
   doing a left and right alignment.
 */
struct a_res_str *
bd_malign (const unsigned char *aa0, int n0,
	   const unsigned char *aa1, int n1,
	   int score_thresh, int max_res,
	   struct pstruct *ppst,
	   struct f_struct *f_str,
	   struct a_res_str *cur_ares, int first_align)
{
  struct a_res_str *tmpl_ares, *tmpr_ares, *this_ares, *last_ares;
  int hoff, l_min, l_max, window;
  unsigned char *local_aa1;
  int local_n1;
  /*
    unsigned char aa1_min_s, aa1_max_s;
  */
  struct rstruct rst;
  int score_ix;
  char save_res;
  int max_sub_score = -1;
  int min_alen;
  int have_local_aa1;
  struct score_count_s s_info = {0, 0, 0};

  min_alen = min(MIN_LOCAL_LEN,n0);

  score_ix = ppst->score_ix;

  /* now we need alignment storage - get it */
  if ((cur_ares->res = (int *)calloc((size_t)max_res,sizeof(int)))==NULL) {
    fprintf(stderr,"*** ERROR [%s:%d] - cannot allocate alignment results array %d\n",
	    __FILE__, __LINE__, max_res);
    exit(1);
  }

  cur_ares->next = NULL;
  cur_ares->n1 = n1;
    
  /* lots of changes to optcut, optwid, and optflag were made in do_walign */

  do_fasta(aa0, n0, aa1, n1, ppst, f_str, &cur_ares->rst, &hoff,0, &s_info);

  if (first_align || cur_ares->rst.score[score_ix]>score_thresh) {
    if (ppst->sw_flag) {
      /* hoff gives us a projection of the query on the library
	 sequence, which can be used to limit the portions of the
	 library sequence that will be aligned by Smith-Waterman

	 this ensures that the optimal alignment (with pam2p[0], not
	 pam2p[1]), aligns the correct region when the only difference between
	 the regions is the lseg encoding
      */

      window = min(n1, ppst->param_u.fa.optwid);
      /*  this windowing seems inappropriate when Smith-Waterman is used (sw_flag),
	  but it is done to ensure that seg'ed regions are ignored 
      */
      l_min = 0;
      l_max = n1;
      if (ppst->pam_x_set) {
	l_min = max(0, -window-hoff);
	l_max = min(n1, n0-hoff+window);
      }

      have_local_aa1 = 0;
      local_aa1 = (unsigned char *)aa1;
      if (l_min > 0 || l_max < n1 - 1) {
	if (l_max - l_min < 0) {
	  fprintf(stderr,"*** ERROR [%s:%d] - l_min: %d > l_max %d\n",__FILE__, __LINE__, l_min,l_max);
	  exit(1);
	}
	if ((local_aa1 = (unsigned char *)calloc(l_max - l_min +2,sizeof(unsigned char *)))==NULL) {
	  fprintf(stderr,"*** ERROR [%s:%d] - cannot allocate local_aa1\n",__FILE__, __LINE__);
	  exit(1);
	}
	  
	local_aa1++;
	memcpy(local_aa1,aa1+l_min, l_max - l_min);
	have_local_aa1 = 1;
      }
      /*
	if (l_min > 0) {
	aa1_min_s = aa1[l_min-1];
	aa1[l_min-1] = '\0';
	}
	if (l_max < n1 - 1) {
	aa1_max_s = aa1[l_max];
	aa1[l_max] = '\0';
	}
      */
      cur_ares->sw_score = sw_walign(f_str->pam2p[0], n0, local_aa1, l_max - l_min,
				     -ppst->gdelval, -ppst->ggapval,
				     f_str->ss, cur_ares);

      cur_ares->min1 += l_min;
      cur_ares->max1 += l_min;
      cur_ares->n1 += l_min;
      if (have_local_aa1) {free(--local_aa1);}
      /*
	if (l_min > 0) { aa1[l_min-1] = aa1_min_s;}
	if (l_max < n1 - 1) { aa1[l_max] = aa1_max_s;}
      */
    }
    else {
      cur_ares->sw_score = bd_walign(aa0, n0, aa1, n1, ppst, f_str, hoff, cur_ares);
    }
  }
  else {
    cur_ares->nres = 0;
    cur_ares->sw_score=0;
    cur_ares->rst.score[0] = cur_ares->rst.score[1] = cur_ares->rst.score[2] = 0;
  }

  /* check to see if a variant is better */

  if (!ppst->do_rep || cur_ares->rst.score[score_ix] <= score_thresh) {return cur_ares;}

  if (cur_ares->min1 >= min_alen) { /* try the left  */
    /* allocate a_res */	
    tmpl_ares = (struct a_res_str *)calloc(1, sizeof(struct a_res_str));

    local_aa1 = (unsigned char *)calloc(cur_ares->min1+2,sizeof(unsigned char));
    local_aa1++;
    memcpy(local_aa1,aa1,cur_ares->min1);
    /*      
	    save_res = aa1[cur_ares->min1-1];
	    aa1[cur_ares->min1-1] = '\0';
    */
    tmpl_ares = bd_malign(aa0, n0, local_aa1, cur_ares->min1,
			  score_thresh, max_res, 
			  ppst, f_str, tmpl_ares,0);

    free(--local_aa1);
    /*
      aa1[cur_ares->min1-1] = save_res;
    */

    /*      fprintf(stderr," bd_malign ret(%d-%d) left: %d-%d/%d-%d : %d\n",
	    cur_ares->min1, cur_ares->max1,
	    tmpl_ares->min0,tmpl_ares->max0,
	    tmpl_ares->min1,tmpl_ares->max1,
	    tmpl_ares->rst.score[score_ix]);
    */

    if (tmpl_ares->rst.score[score_ix] > score_thresh) {
      max_sub_score = tmpl_ares->rst.score[score_ix];
    }
    else {
      if (tmpl_ares->res) free(tmpl_ares->res);
      free(tmpl_ares);
      tmpl_ares=NULL;
    }
  }
  else {tmpl_ares = NULL;}

  local_n1 = n1 - cur_ares->max1;
  if (local_n1 >= min_alen) { /* try the right  */      
    /* allocate a_res */	
    tmpr_ares = (struct a_res_str *)calloc(1, sizeof(struct a_res_str));

    /* find boundaries */
    local_aa1 = (unsigned char *)calloc(local_n1+2, sizeof(unsigned char));
    local_aa1++;
    memcpy(local_aa1,aa1+cur_ares->max1, local_n1);
    /*
      save_res = aa1[cur_ares->max1-1];
      aa1[cur_ares->max1-1] = '\0';
    */
    tmpr_ares = bd_malign(aa0, n0, local_aa1, local_n1,
			  score_thresh, max_res,
			  ppst, f_str, tmpr_ares,0);
    free(--local_aa1);
    /*
      aa1[cur_ares->max1-1] = save_res;
    */
    /*       fprintf(stderr," bd_malign ret(%d-%d) right: %d-%d/%d-%d : %d\n",
	     cur_ares->min1, cur_ares->max1,
	     tmpr_ares->min0,tmpr_ares->max0,
	     tmpr_ares->min1,tmpr_ares->max1,
	     tmpr_ares->rst.score[score_ix]);
    */

    if (tmpr_ares->rst.score[score_ix] > score_thresh) {
      /* adjust the left boundary */
      for (this_ares = tmpr_ares; this_ares; this_ares = this_ares->next) {
	this_ares->min1 += cur_ares->max1;
	this_ares->max1 += cur_ares->max1;
	last_ares=this_ares;
      }
      if (tmpr_ares->rst.score[score_ix] > max_sub_score) {
	max_sub_score = tmpr_ares->rst.score[score_ix];
      }
    }
    else {
      if (tmpr_ares->res) free(tmpr_ares->res);
      free(tmpr_ares);
      tmpr_ares = NULL;
    }
  }
  else {tmpr_ares = NULL;}

  if (max_sub_score <= score_thresh) return cur_ares;

  /*
    fprintf(stderr, "lr: %d l: %d r:%d\n",
    max_sub_score,
    (tmpl_ares ? tmpl_ares->rst.score[score_ix] : 0),
    (tmpr_ares ? tmpr_ares->rst.score[score_ix] : 0));
  */

  cur_ares = merge_ares_chains(cur_ares, tmpl_ares, score_ix, "left");
  cur_ares = merge_ares_chains(cur_ares, tmpr_ares, score_ix, "right");

  return cur_ares;
}

/* A PACKAGE FOR LOCALLY ALIGNING TWO SEQUENCES WITHIN A BAND:

   To invoke, call LOCAL_ALIGN(A,B,M,N,L,U,W,G,H,S,dflag,&SI,&SJ,&EI,&EJ,MW).
   The parameters are explained as follows:
	A, B : two sequences to be aligned
	M : the length of sequence A
	N : the length of sequence B
	L : lower bound of the band
	U : upper bound of the band
	W : scoring table for matches and mismatches
	G : gap-opening penalty
	H : gap-extension penalty
	dflag : 0 - no display or backward pass
	*SI : starting position of sequence A in the optimal local alignment
	*SJ : starting position of sequence B in the optimal local alignment
	*EI : ending position of sequence A in the optimal local alignment
	*EJ : ending position of sequence B in the optimal local alignment
	MW  : maximum window size
*/

int bd_walign (const unsigned char *aa0, int n0,
	       const unsigned char *aa1, int n1,
	       struct pstruct *ppst, 
	       struct f_struct *f_str, int hoff,
	       struct a_res_str *a_res)
{
   int low, up, score;
   int min0, min1, max0, max1;
   int window;

   window = min (n1, ppst->param_u.fa.optwid);
   if (window > f_str->bss_size) {
     fprintf(stderr,"*** ERROR [%s:%d] - walign window [%d] out of range [%d]\n",
	     __FILE__, __LINE__, window, f_str->bss_size);
     window = f_str->bss_size - 4;
   }

   /* hoff is the offset found from aa1 to seq 2 by dmatch */
   low = -window/2-hoff;
   up = low+window;

   score=LOCAL_ALIGN(aa0-1,aa1-1,n0,n1, low, up,
		    ppst->pam2[0], -ppst->gdelval, -ppst->ggapval,
		    &min0,&min1,&max0,&max1,ppst->param_u.fa.optwid,f_str);
   a_res->n1 = n1;
  
   if (score <=0) {
     fprintf(stderr,"*** [%s:%d] - score <= 0 - n0/n1: %d/%d hoff: %d window: %d\n",
	     __FILE__, __LINE__, n0, n1, hoff, window);
     return 0;
   }
  
/*
  fprintf(stderr," ALIGN: start0: %d start1: %d stop0: %d stop1: %d, bot: %d top: %d, win: %d MX %d\n",
  min0-1,min1-1,max0-min0+1,max1-min1+1,low-(min1-min0),up-(min1-min0),
  ppst->param_u.fa.optwid,n0);
*/

   a_res->min0 = min0-1; a_res->min1 = min1-1;
   a_res->max0 = max0; a_res->max1 = max1; 

   B_ALIGN(aa0-1+min0-1,aa1-1+min1-1,max0-min0+1,max1-min1+1,
	   low-(min1-min0),up-(min1-min0),
	   ppst->pam2[0], -ppst->gdelval, -ppst->ggapval,
	   a_res->res,&a_res->nres,ppst->param_u.fa.optwid,n0,f_str->bss, &f_str->mtp);

   return score;
}

static int
LOCAL_ALIGN(const unsigned char *A, const unsigned char *B,
	    int M, int N,
	    int low, int up, int **W, int G,int H,
	    int *psi, int *psj, int *pei, int *pej, int MW,
	    struct f_struct *f_str)
{ 
  int band;
  register struct bdstr *bssp;
  int i, j, si, ei;
  int c, d, e, t, m;
  int leftd, rightd;
  int best_score, starti, startj, endi, endj;
  int *wa, curd;
  int ib;
  char flag;
  
  bssp = f_str->bss;

  m = G+H;
  low = max(-M, low);
  up = min(N, up);
  
  if (N <= 0) { 
    *psi = *psj = *pei = *pej;
    return 0;
  }
  if (M <= 0) {
    *psi = *psj = *pei = *pej;
    return 0;
  }
  band = up-low+1;
  if (band < 1) {
    fprintf(stderr,"*** ERROR [%s:%d] low > up is unacceptable!: M: %d N: %d l/u: %d/%d\n",
	    __FILE__, __LINE__,  M, N, low, up);
    return -1;
  }

  /* already done by init_work(); 
  j = (MW + 2 + 2) * sizeof(struct bdstr);
  if (f_str->bss==NULL) f_str->bss = (struct bdstr *) ckalloc(j);
  */

  if (low > 0) leftd = 1;
  else if (up < 0) leftd = band;
  else leftd = 1-low;
  rightd = band;
  si = max(0,-up);
  ei = min(M,N-low);
  bssp[leftd].CC = 0;
  for (j = leftd+1; j <= rightd; j++) {
    bssp[j].CC = 0;
    bssp[j].DD = -G;
  }
  bssp[rightd+1].CC = MININT;
  bssp[rightd+1].DD = MININT;
  best_score = 0;
  endi = si;
  endj = si+low;
  bssp[leftd-1].CC = MININT;
  bssp[leftd].DD = -G;
  for (i = si+1; i <= ei; i++) {
    if (i > N-up) rightd--;
    if (leftd > 1) leftd--;
    wa = W[A[i]];
    if ((c = bssp[leftd+1].CC-m) > (d = bssp[leftd+1].DD-H)) d = c;
    if ((ib = leftd+low-1+i ) > 0) c = bssp[leftd].CC+wa[B[ib]];
/*
    if (ib > N) fprintf(stderr,"B[%d] out of range %d\n",ib,N);
*/
    if (d > c) c = d;
    if (c < 0) c = 0;
    e = c-G;
    bssp[leftd].DD = d;
    bssp[leftd].CC = c;
    if (c > best_score) {
      best_score = c;
      endi = i;
      endj = ib;
    }
    for (curd=leftd+1; curd <= rightd; curd++) {
      if ((c = c-m) > (e = e-H)) e = c;
      if ((c = bssp[curd+1].CC-m) > (d = bssp[curd+1].DD-H)) d = c;
/*
      if ((ib=curd+low-1+i) <= 0 || ib > N)
	fprintf(stderr,"B[%d]:%d\n",ib,B[ib]);
*/
      c = bssp[curd].CC + wa[B[curd+low-1+i]];
      if (e > c) c = e;
      if (d > c) c = d;
      if (c < 0) c = 0;
      bssp[curd].CC = c;
      bssp[curd].DD = d;
      if (c > best_score) {
	best_score = c;
	endi = i;
	endj = curd+low-1+i;
      }
    }
  }
  
  leftd = max(1,-endi-low+1);
  rightd = band-(up-(endj-endi));
  bssp[rightd].CC = 0;
  t = -G;
  for (j = rightd-1; j >= leftd; j--) {
    bssp[j].CC = t = t-H;
    bssp[j].DD = t-G;
  }
  for (j = rightd+1; j <= band; ++j) bssp[j].CC = MININT;
  bssp[leftd-1].CC = bssp[leftd-1].DD = MININT;
  bssp[rightd].DD = -G;
  flag = 0;
  for (i = endi; i >= 1; i--) {
    if (i+low <= 0) leftd++;
    if (rightd < band) rightd++;
    wa = W[A[i]];
    if ((c = bssp[rightd-1].CC-m) > (d = bssp[rightd-1].DD-H)) d = c;
    if ((ib = rightd+low-1+i) <= N) c = bssp[rightd].CC+wa[B[ib]];

/*
    if (ib <= 0) fprintf(stderr,"rB[%d] <1\n",ib);
*/
    if (d > c) c = d;
    e = c-G;
    bssp[rightd].DD = d;
    bssp[rightd].CC = c;
    if (c == best_score) {
      starti = i;
      startj = ib;
      flag = 1;
      break;
    }
    for (curd=rightd-1; curd >= leftd; curd--) {
      if ((c = c-m) > (e = e-H)) e = c;
      if ((c = bssp[curd-1].CC-m) > (d = bssp[curd-1].DD-H)) d = c;

/*
      if ((ib=curd+low-1+i) <= 0 || ib > N)
	fprintf(stderr,"i: %d, B[%d]:%d\n",i,ib,B[ib]);
*/
      c = bssp[curd].CC + wa[B[curd+low-1+i]];
      if (e > c) c = e;
      if (d > c) c = d;
      bssp[curd].CC = c;
      bssp[curd].DD = d;
      if (c == best_score) {
	starti = i;
	startj = curd+low-1+i;
	flag = 1;
	break;
      }
    }
    if (flag == 1) break;
  }
  
  if (starti < 0 || starti > M || startj < 0 || startj > N) {
    printf("starti=%d, startj=%d\n",starti,startj);
    *psi = *psj = *pei = *pej;
    exit(1);
  }
  *psi = starti;
  *psj = startj;
  *pei = endi;
  *pej = endj;
  return best_score;
}

/* A PACKAGE FOR GLOBALLY ALIGNING TWO SEQUENCES WITHIN A BAND:

   To invoke, call B_ALIGN(A,B,M,N,L,U,W,G,H,S,MW,MX).
   The parameters are explained as follows:
	A, B : two sequences to be aligned
	M : the length of sequence A
	N : the length of sequence B
	L : lower bound of the band
	U : upper bound of the band
	W : scoring table for matches and mismatches
	G : gap-opening penalty
	H : gap-extension penalty
	S : script for DISPLAY routine
	MW : maximum window size
	MX : maximum length sequence M to be aligned
*/

/* static int IP; */
/* static int *MP[3]; */		/* save crossing points */
/* static int *FP; */			/* forward dividing points */
/* static char *MT[3]; */		/* 0: rep, 1: del, 2: ins */
/* static char *FT; */

#define gap(k)  ((k) <= 0 ? 0 : g+h*(k))	/* k-symbol indel cost */

/* Append "Delete k" op */
#define DEL(k)				\
{ if (*last < 0)			\
    *last = (*sapp)[-1] -= (k);		\
  else {				\
    *last = (*sapp)[0] = -(k);		\
    (*sapp)++;				\
  }					\
}

/* Append "Insert k" op */
#define INS(k)				\
{ if (*last > 0)			\
    *last = (*sapp)[-1] += (k);		\
  else {				\
    *last = (*sapp)[0] = (k);		\
    (*sapp)++;				\
  }					\
}

#define REP { *last = (*sapp)[0] = 0; (*sapp)++;} /* Append "Replace" op */

/* bg_align(A,B,M,N,up,low,tb,te) returns the cost of an optimum conversion between
  A[1..M] and B[1..N] and appends such a conversion to the current script.
  tb(te)= 1  no gap-open penalty if the conversion begins(ends) with a delete.
  tb(te)= 2  no gap-open penalty if the conversion begins(ends) with an insert.
*/
static int
bg_align(const unsigned char *A, const unsigned char *B, 
	 int M, int N,
	 int low, int up, int tb, int te,
	 int **w, int g, int h,
	 struct bdstr *bss, struct mtp_str *mtp,
	 int **sapp, int *last)
{
  int rmid, k, l, r, v, kt;
  int t1, t2, t3;

  {
  int band, midd;
  int leftd, rightd;	/* for CC, DD, CP and DP */
  register int curd;	/* current index for CC, DD CP and DP */
  register int i, j;
  register int c, d, e;
  int t, fr, *wa, ib, m;

  /* Boundary cases: M <= 0 , N <= 0, or up-low <= 0 */
  if (N <= 0) { 
    if (M > 0) { DEL(M) }
    return 0;
  }
  if (M <= 0) {
    INS(N)
    return 0;
  }
  if ((band = up-low+1) <= 1) {
    for (i = 1; i <= M; i++) { REP }
    return 0;
  }

  /* Divide: Find all crossing points */

  /* Initialization */
  m = g + h;

  midd = band/2 + 1;
  rmid = low + midd - 1;
  leftd = 1-low;
  rightd = up-low+1;
  if (leftd < midd) {
    fr = -1;
    for (j = 0; j < midd; j++) 
      bss[j].CP = bss[j].DP = -1;
    for (j = midd; j <= rightd; j++) {
      bss[j].CP = bss[j].DP = 0;
    }
    mtp->MP[0][0] = -1;
    mtp->MP[1][0] = -1;
    mtp->MP[2][0] = -1;
    mtp->MT[0][0] = mtp->MT[1][0] = mtp->MT[2][0] = 0;
  } else if (leftd > midd) {
    fr = leftd-midd;
    for (j = 0; j <= midd; j++) {
      bss[j].CP = bss[j].DP = fr;
    }
    for (j = midd+1; j <= rightd; j++) 
      bss[j].CP = bss[j].DP = -1;
    mtp->MP[0][fr] = -1;
    mtp->MP[1][fr] = -1;
    mtp->MP[2][fr] = -1;
    mtp->MT[0][fr] = mtp->MT[1][fr] = mtp->MT[2][fr] = 0;
  } else {
    fr = 0;
    for (j = 0; j < midd; j++) {
      bss[j].CP = bss[j].DP = 0;
    }
    for (j = midd; j <= rightd; j++) {
      bss[j].CP = bss[j].DP = 0;
    }
    mtp->MP[0][0] = -1;
    mtp->MP[1][0] = -1;
    mtp->MP[2][0] = -1;
    mtp->MT[0][0] = mtp->MT[1][0] = mtp->MT[2][0] = 0;
  }

  bss[leftd].CC = 0;
  if (tb == 2) t = 0;
  else t = -g;
  for (j = leftd+1; j <= rightd; j++) {
    bss[j].CC = t = t-h;
    bss[j].DD = t-g;
  }
  bss[rightd+1].CC = MININT;
  bss[rightd+1].DD = MININT;
  if (tb == 1) bss[leftd].DD = 0;
  else bss[leftd].DD = -g;
  bss[leftd-1].CC = MININT;
  for (i = 1; i <= M; i++) {
    if (i > N-up) rightd--;
    if (leftd > 1) leftd--;
    wa = w[A[i]];
    if ((c = bss[leftd+1].CC-m) > (d = bss[leftd+1].DD-h)) {
      d = c;
      bss[leftd].DP = bss[leftd+1].CP;
    } else bss[leftd].DP = bss[leftd+1].DP;
    if ((ib = leftd+low-1+i) > 0) c = bss[leftd].CC+wa[B[ib]];
    if (d > c || ib <= 0) {
      c = d;
      bss[leftd].CP = bss[leftd].DP;
    }
    e = c-g;
    bss[leftd].DD = d;
    bss[leftd].CC = c;
    mtp->IP = bss[leftd].CP;
    if (leftd == midd) bss[leftd].CP = bss[leftd].DP = mtp->IP = i;
    for (curd=leftd+1; curd <= rightd; curd++) {
      if (curd != midd) {
	if ((c = c-m) > (e = e-h)) {
	  e = c;
	  mtp->IP = bss[curd-1].CP;
	}  /* otherwise, mtp->IP is unchanged */
	if ((c = bss[curd+1].CC-m) > (d = bss[curd+1].DD-h)) {
	  d = c;
	  bss[curd].DP = bss[curd+1].CP;
	} else {
	  bss[curd].DP = bss[curd+1].DP;
	}
	c = bss[curd].CC + wa[B[curd+low-1+i]];
	if (c < d || c < e) {
	  if (e > d) {
	    c = e;
	    bss[curd].CP = mtp->IP;
	  } else {
	    c = d;
	    bss[curd].CP = bss[curd].DP;
	  }
	} /* otherwise, CP is unchanged */
	bss[curd].CC = c;
	bss[curd].DD = d;
      } else {
	if ((c = c-m) > (e = e-h)) {
	  e = c;
	  mtp->MP[1][i] = bss[curd-1].CP;
	  mtp->MT[1][i] = 2;
	} else {
	  mtp->MP[1][i] = mtp->IP;
	  mtp->MT[1][i] = 2;
	}
	if ((c = bss[curd+1].CC-m) > (d = bss[curd+1].DD-h)) {
	  d = c;
	  mtp->MP[2][i] = bss[curd+1].CP;
	  mtp->MT[2][i] = 1;
	} else {
	  mtp->MP[2][i] = bss[curd+1].DP;
	  mtp->MT[2][i] = 1;
	}
	c = bss[curd].CC + wa[B[curd+low-1+i]];
	if (c < d || c < e) {
	  if (e > d) {
	    c = e;
	    mtp->MP[0][i] = mtp->MP[1][i];
	    mtp->MT[0][i] = 2;
	  } else {
	    c = d;
	    mtp->MP[0][i] = mtp->MP[2][i];
	    mtp->MT[0][i] = 1;
	  }
	} else {
	  mtp->MP[0][i] = i-1;
	  mtp->MT[0][i] = 0;
	}
	if (c-g > e) {
	  mtp->MP[1][i] = mtp->MP[0][i];
	  mtp->MT[1][i] = mtp->MT[0][i];
	}
	if (c-g > d) {
	  mtp->MP[2][i] = mtp->MP[0][i];
	  mtp->MT[2][i] = mtp->MT[0][i];
	}
	bss[curd].CP = bss[curd].DP = mtp->IP = i;
	bss[curd].CC = c;
	bss[curd].DD = d;
      }
    }
  }

  /* decide which path to be traced back */
  if (te == 1 && d+g > c) {
    k = bss[rightd].DP;
    l = 2;
  } else if (te == 2 && e+g > c) {
    k = mtp->IP;
    l = 1;
  } else {
    k = bss[rightd].CP;
    l = 0;
  }
  if (rmid > N-M) l = 2;
  else if (rmid < N-M) l = 1;
  v = c;
  }
  /* Conquer: Solve subproblems recursively */

  /* trace back */
  r = -1;	
  for (; k > -1; r=k, k=mtp->MP[l][r], l=mtp->MT[l][r]){
    mtp->FP[k] = r;
    mtp->FT[k] = l;	/* l=0,1,2 */
  }
  /* forward dividing */
  if (r == -1) { /* optimal alignment did not cross the middle diagonal */
    if (rmid < 0) {
      bg_align(A,B,M,N,rmid+1,up,tb,te,w,g,h,bss, mtp, sapp, last);
    }
    else {
      bg_align(A,B,M,N,low,rmid-1,tb,te,w,g,h,bss, mtp, sapp, last);
    }
  } else {
    k = r;
    l = mtp->FP[k];
    kt = mtp->FT[k];

    /* first block */
    if (rmid < 0) {
      bg_align(A,B,r-1,r+rmid,rmid+1,min(up,r+rmid),tb,1,w,g,h,bss, mtp, sapp,last);
      DEL(1)
    } else if (rmid > 0) {
      bg_align(A,B,r,r+rmid-1,max(-r,low),rmid-1,tb,2,w,g,h,bss, mtp, sapp,last);
      INS(1)
    }

    /* intermediate blocks */
    t2 = up-rmid-1;
    t3 = low-rmid+1;
    for (; l > -1; k = l, l = mtp->FP[k], kt = mtp->FT[k]) {
      if (kt == 0) { REP }
      else if (kt == 1) { /* right-hand side triangle */
	INS(1)
	t1 = l-k-1;
	bg_align(A+k,B+k+rmid+1,t1,t1,0,min(t1,t2),2,1,w,g,h,bss, mtp, sapp,last);
	DEL(1)
      }
      else { /* kt == 2, left-hand side triangle */
	DEL(1)
	t1 = l-k-1;
	bg_align(A+k+1,B+k+rmid,t1,t1,max(-t1,t3),0,1,2,w,g,h,bss, mtp, sapp,last);
	INS(1)
      }
    }

    /* last block */
    if (N-M > rmid) {
      INS(1)
      t1 = k+rmid+1;
      bg_align(A+k,B+t1,M-k,N-t1,0,min(N-t1,t2),2,te,w,g,h,bss, mtp, sapp,last);
    } else if (N-M < rmid) {
      DEL(1)
      t1 = M-(k+1);
      bg_align(A+k+1,B+k+rmid,t1,N-(k+rmid),max(-t1,t3),0,1,te,w,g,h,
	       bss,mtp,sapp,last);
    }
  }
  return(v);
}

int B_ALIGN(const unsigned char *A, const unsigned char *B,
	    int M, int N,
	    int low, int up, int **W, int G, int H, int *S, int *nS,
	    int MW, int MX, struct bdstr *bss, struct mtp_str *mtp)
{ 
  int c, i, j;
  int g, h;
  size_t mj;
  int check_score;
  int **sapp, *sapp_v, *last, last_v;

  g = G;
  h = H;
  sapp_v = S;
  sapp = &sapp_v;

  last_v = 0;
  last = &last_v;

  low = min(max(-M, low),min(N-M,0));
  up = max(min(N, up),max(N-M,0));

  if (N <= 0) { 
    if (M > 0) { DEL(M); }
    return -gap(M);
  }
  if (M <= 0) {
    INS(N);
    return -gap(N);
  }
  if (up-low+1 <= 1) {
    c = 0;
    for (i = 1; i <= M; i++) {
      REP;
      c += W[A[i]][B[i]];
    }
    return c;
  }

  if (mtp->MT[0]==NULL) {
    mj = (MX+1) * sizeof(int);
    mtp->MT[0] = (int *) ckalloc(mj);
    mtp->MT[1] = (int *) ckalloc(mj);
    mtp->MT[2] = (int *) ckalloc(mj);
    mtp->FT = (int *) ckalloc(mj);

    /* mj *= sizeof(int); (already done) */
    mtp->MP[0] = (int *) ckalloc(mj);
    mtp->MP[1] = (int *) ckalloc(mj);
    mtp->MP[2] = (int *) ckalloc(mj);
    mtp->FP = (int *) ckalloc(mj);
  }

  c = bg_align(A,B,M,N,low,up,0,0,W,G,H,bss, mtp, sapp, last);

  check_score = BCHECK_SCORE(A,B,M,N,S,W,G,H,nS);

  free(mtp->FP); free(mtp->MP[2]); free(mtp->MP[1]); free(mtp->MP[0]);
  free(mtp->FT); free(mtp->MT[2]); free(mtp->MT[1]); free(mtp->MT[0]);
  mtp->MT[0]=NULL;

  if (check_score != c)
    printf("\nBCheck_score=%d != %d\n", check_score,c);
  return c;
}

struct a_res_str *
do_walign (const unsigned char *aa0, int n0,
	   const unsigned char *aa1, int n1,
	   int frame, int repeat_thresh,
	   struct pstruct *ppst, 
	   struct f_struct *f_str, 
	   int *have_ares)
{
  int hoff, use_E_thresholds_s, optflag_s, optcut_s, optwid_s, n10, score;
  const unsigned char *aa1p;
  struct rstruct rst;
  struct a_res_str *a_res, *tmp_a_res;
  int a_res_index;
#ifdef DEBUG
  unsigned long adler32_crc;
#endif

#ifdef TFASTA  
  f_str->n10 = n10=aatran(aa1,f_str->aa1x,n1,frame);
  aa1p = f_str->aa1x;
#else
  n10 = n1;
  aa1p = aa1;
#endif

#ifdef DEBUG
  adler32_crc = adler32(1L,aa1,n1);
#endif

  *have_ares = 0x3;	/* set 0x2 bit to indicate local copy */

  if ((a_res = (struct a_res_str *)calloc(1, sizeof(struct a_res_str)))==NULL) {
    fprintf(stderr,"*** ERROR [%s:%d] - cannot allocate a_res", __FILE__, __LINE__);
    return NULL;
  }

  use_E_thresholds_s = ppst->param_u.fa.use_E_thresholds;
  optflag_s = ppst->param_u.fa.optflag;
  optcut_s = ppst->param_u.fa.optcut;
  optwid_s = ppst->param_u.fa.optwid;
  ppst->param_u.fa.use_E_thresholds = 0;
  ppst->param_u.fa.optflag = 1;
  ppst->param_u.fa.optcut = 0;
  if (!ppst->param_u.fa.optwid_set) {
    ppst->param_u.fa.optwid *= 2;
  }

  a_res = bd_malign(aa0, n0, aa1p, n10,
		    repeat_thresh, f_str->max_res,
		    ppst, f_str, a_res,1);

  ppst->param_u.fa.use_E_thresholds = use_E_thresholds_s;
  ppst->param_u.fa.optflag = optflag_s;
  ppst->param_u.fa.optcut = optcut_s;
  ppst->param_u.fa.optwid = optwid_s;

#ifdef DEBUG
  if (adler32(1L,aa1,n1) != adler32_crc) {
    fprintf(stderr,"*** ERROR [%s:%d] - adler32_crc mismatch n1: %d\n",__FILE__, __LINE__, n1);
  }
#endif

  a_res_index = 0;
  for (tmp_a_res=a_res; tmp_a_res; tmp_a_res = tmp_a_res->next) {
    tmp_a_res->index = a_res_index++;
  }

  return a_res;
}

void
pre_cons(const unsigned char *aa1, int n1, int frame, struct f_struct *f_str) {

#ifdef TFASTA
  f_str->n10 = aatran(aa1,f_str->aa1x,n1,frame);
#endif
}

/* aln_func_vals - set up aln.qlfact, qlrev, llfact, llmult, frame, llrev */
/* call from calcons, calc_id, calc_code */
void 
aln_func_vals(int frame, struct a_struct *aln) {

#ifdef TFASTA
  aln->qlfact = 1;
  aln->llfact = 3;
  aln->llmult = 3;
  aln->qlrev = 0;
  aln->frame = frame;
  if (frame > 2) {
    aln->llrev = 1;
    aln->frame = 3 - frame;
  }
  else aln->llrev = 0;
#else	/* FASTA */
  aln->llfact = aln->qlfact = aln->llmult = 1;
  aln->llrev = 0;
  if (frame > 0) aln->qlrev = 1;
  else aln->qlrev = 0;
  aln->frame = 0;
#endif
}
