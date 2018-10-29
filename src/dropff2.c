/* $Id: dropff2.c 989 2012-07-24 19:37:38Z wrp $ */

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

/* this code implements the "fastf" algorithm, which is designed to
   deconvolve mixtures of protein sequences derived from mixed-peptide
   Edman sequencing.  The expected input is:

   >test | 40001 90043 | mgstm1
   MGCEN,
   MIDYP,
   MLLAY,
   MLLGY

   Where the ','s indicate the length/end of the sequencing cycle
   data.  Thus, in this example, the sequence is from a mixture of 4
   peptides, M was found in the first position, G,I, and L(2) at the second,
   C,D, L(2) at the third, etc.

   Because the sequences are derived from mixtures, there need not be
   any partial sequence "MGCEN", the actual deconvolved sequence might be
   "MLDGN".
*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include <ctype.h>

#include "defs.h"
#include "param.h"
#include "structs.h"
#include "tatstats.h"

#define EOSEQ 0 
#define ESS 59
#define MAXHASH 32
#define NMAP MAXHASH+1
#define NMAP_X 23	/* re-code NMAP for 'X' */
#define NMAP_Z 24	/* re-code NMAP for '*' */

#ifndef MAXSAV
#define MAXSAV 10
#endif

#define DROP_INTERN
#include "drop_func.h"

static char *verstr="4.21 May 2006 (ajm/wrp)";

int shscore(unsigned char *aa0, const int n0, int **pam2, int nsq);

#ifdef TFAST
extern int aatran(const unsigned char *ntseq, unsigned char *aaseq, 
		  const int maxs, const int frame);
#endif

struct hlstr { int next, pos;};

void savemax(struct dstruct *, struct f_struct *);

static int m0_spam(unsigned char *, const unsigned char *, int, struct savestr *,
	    int **, struct f_struct *);
static int m1_spam(unsigned char *, int,
		   const unsigned char *, int,
		   struct savestr *, int **, int, struct f_struct *);

int sconn(struct savestr **v, int nsave, int cgap,
	  struct f_struct *, struct rstruct *, const struct pstruct *,
	  const unsigned char *aa0, int n0,
	  const unsigned char *aa1, int n1,
	  int opt_prob);

void kpsort(struct savestr **, int);
void kssort(struct savestr **, int);
void kpsort(struct savestr **, int);

int
sconn_a(unsigned char *, int, int, struct f_struct *,
	struct a_res_str *);

/* initialize for fasta */

void
init_work (unsigned char *aa0, int n0, 
	   struct pstruct *ppst,
	   struct f_struct **f_arg)
{
   int mhv, phv;
   int hmax;
   int i0, ii0, hv;
   struct f_struct *f_str;

   int maxn0;
   int i, j, q;
   struct savestr *vmptr;
   int *res;
   int nsq;

   nsq = ppst->nsqx;

   f_str = (struct f_struct *) calloc(1, sizeof(struct f_struct));
   if(f_str == NULL) {
     fprintf(stderr, "*** error [%s:%d] - cannot calloc f_str [%lu]\n",
	     __FILE__, __LINE__, sizeof(struct f_struct));
     exit(1);
   }

   ppst->sw_flag = 0;

   /* fastf3 cannot work with lowercase symbols as low complexity;
      thus, NMAP must be disabled; this depends on aascii['X']  */
   if (ppst->hsq[NMAP_X] == NMAP ) {ppst->hsq[NMAP_X]=1;}
   if (ppst->hsq[NMAP_Z] == NMAP ) {ppst->hsq[NMAP_Z]=1;}

   /*   this does not work for share ppst structs, as in threads */
   /*else {fprintf(stderr," cannot find 'X'==NMAP\n");} */

   for (i0 = 1, mhv = -1; i0 <= ppst->nsq; i0++)
      if (ppst->hsq[i0] < NMAP && ppst->hsq[i0] > mhv) mhv = ppst->hsq[i0];

   if (mhv <= 0) {
      fprintf (stderr, "*** error [%s:%d] -  maximum hsq <=0 %d\n",
	       __FILE__, __LINE__, mhv);
      exit (1);
   }

   for (f_str->kshft = 0; mhv > 0; mhv /= 2)
      f_str->kshft++;

/*      kshft = 2;	*/
   hmax = hv = (1 << f_str->kshft);
   f_str->hmask = (hmax >> f_str->kshft) - 1;

   if ((f_str->aa0 = (unsigned char *) calloc(n0+1, sizeof(char))) == NULL) {
     fprintf (stderr, "*** error [%s:%d] - cannot allocate f_str->aa0 array; %d\n",
	      __FILE__, __LINE__, n0+1);
     exit (1);
   }
   for (i=0; i<n0; i++) f_str->aa0[i] = aa0[i];
   aa0 = f_str->aa0;

   if ((f_str->aa0t = (unsigned char *) calloc(n0+1, sizeof(char))) == NULL) {
     fprintf (stderr, "*** error [%s:%d] - cannot allocate f_str0->aa0t array; %d\n",
	      __FILE__, __LINE__, n0+1);
     exit (1);
   }
   f_str->aa0ix = 0;

   if ((f_str->harr = (struct hlstr *) calloc (hmax, sizeof (struct hlstr))) == NULL) {
     fprintf (stderr, "*** error [%s:%d] -  cannot allocate hash array; hmax: %d hmask: %d\n",
	      __FILE__, __LINE__, hmax,f_str->hmask);
     exit (1);
   }
   if ((f_str->pamh1 = (int *) calloc (nsq+1, sizeof (int))) == NULL) {
     fprintf (stderr, "*** error [%s:%d] -  cannot allocate pamh1 array [%d]\n",
	      __FILE__, __LINE__, nsq+1);
     exit (1);
   }
   if ((f_str->pamh2 = (int *) calloc (hmax, sizeof (int))) == NULL) {
     fprintf (stderr, "*** error [%s:%d] - cannot allocate pamh2 array [%d]\n",
	      __FILE__, __LINE__, hmax);
     exit (1);
   }
   if ((f_str->link = (struct hlstr *) calloc (n0, sizeof (struct hlstr))) == NULL) {
     fprintf (stderr, "*** error [%s:%d] - cannot allocate hash link array [%d]",
	      __FILE__, __LINE__, n0);
     exit (1);
   }

   for (i0 = 0; i0 < hmax; i0++) {
      f_str->harr[i0].next = -1;
      f_str->harr[i0].pos = -1;
   }

   for (i0 = 0; i0 < n0; i0++) {
      f_str->link[i0].next = -1;
      f_str->link[i0].pos = -1;
   }

   /* encode the aa0 array */
   /*
     this code has been modified to allow for mixed peptide sequences
      aa0[] = 5 8 9 3 4 NULL 5 12 3 7 2 NULL
      the 'NULL' character resets the hash position counter, to indicate that
      any of several residues can be in the same position.
      We also need to keep track of the number of times this has happened, so that
      we can redivide the sequence later

      i0 counts through the sequence
      ii0 counts through the hashed sequence

      */

   f_str->nm0 = 1;
   f_str->nmoff = -1;
   phv = hv = 0;
   for (i0= ii0 = 0; i0 < n0; i0++, ii0++) {
     /* reset the counter and start hashing again */
     if (aa0[i0] == ESS || aa0[i0] == 0) {
       aa0[i0] = 0;	/* set ESS to 0 */
       /*       fprintf(stderr," converted ',' to 0\n");*/
       i0++;	/* skip over the blank */
       f_str->nm0++;
       if (f_str->nmoff < 0) f_str->nmoff = i0;
       phv = hv = 0;
       ii0 = 0;
     }
     hv = ppst->hsq[aa0[i0]];
     f_str->link[i0].next = f_str->harr[hv].next;
     f_str->link[i0].pos = f_str->harr[hv].pos;
     f_str->harr[hv].next = i0;
     f_str->harr[hv].pos = ii0;
     f_str->pamh2[hv] = ppst->pam2[0][aa0[i0]][aa0[i0]];
   }
   if (f_str-> nmoff < 0) f_str->nmoff = n0;


#ifdef DEBUG
   /*
   fprintf(stderr," nmoff: %d/%d nm0: %d\n", f_str->nmoff, n0,f_str->nm0);
   */
#endif

/*
#ifdef DEBUG
   fprintf(stderr," hmax: %d\n",hmax);
   for ( hv=0; hv<hmax; hv++)
       fprintf(stderr,"%2d %c %3d %3d\n",hv,
	       (hv > 0 && hv < ppst->nsq ) ? ppst->sq[ppst->hsq[hv]] : ' ',
	       f_str->harr[hv].pos,f_str->harr[hv].next);
   fprintf(stderr,"----\n");
   for ( hv=0; hv<n0; hv++)
       fprintf(stderr,"%2d: %3d %3d\n",hv,
	       f_str->link[hv].pos,f_str->link[hv].next);
#endif
*/

   f_str->maxsav = MAXSAV;
   if ((f_str->vmax = (struct savestr *)
	calloc(MAXSAV,sizeof(struct savestr)))==NULL) {
     fprintf(stderr, "*** error [%s:%d] - cannot allocate vmax[%d].\n",
	     __FILE__, __LINE__, f_str->maxsav);
     exit(1);
   }

   if ((f_str->vptr = (struct savestr **)
	calloc(MAXSAV,sizeof(struct savestr *)))==NULL) {
     fprintf(stderr, "*** error [%s:%d] - cannot allocate vptr[%d].\n",
	     __FILE__, __LINE__, f_str->maxsav);
     exit(1);
   }

   for (vmptr = f_str->vmax; vmptr < &f_str->vmax[MAXSAV]; vmptr++) {
     vmptr->used = (int *) calloc(n0, sizeof(int));
     if(vmptr->used == NULL) {
       fprintf(stderr, "*** error [%s:%d] - cannot alloc vmptr->used [%d]\n",
	       __FILE__, __LINE__, n0);
       exit(1);
     }
   }

/* this has been modified from 0..<ppst->nsq to 1..<=ppst->nsq because the
   pam2[0][0] is now undefined for consistency with blast
*/

   for (i0 = 1; i0 <= ppst->nsq; i0++)
     f_str->pamh1[i0] = ppst->pam2[0][i0][i0];

   ppst->param_u.fa.cgap = shscore(aa0,f_str->nmoff-1,ppst->pam2[0],ppst->nsq)/3;
   if (ppst->param_u.fa.cgap > ppst->param_u.fa.bestmax/4)
     ppst->param_u.fa.cgap = ppst->param_u.fa.bestmax/4;

   f_str->ndo = 0;
   f_str->noff = n0-1;
   if (f_str->diag==NULL) 
     f_str->diag = (struct dstruct *) calloc ((size_t)MAXDIAG,
					      sizeof (struct dstruct));

   if (f_str->diag == NULL)
   {
      fprintf (stderr, "*** error [%s:%d] - cannot allocate diagonal arrays: %ld\n",
	       __FILE__, __LINE__, (long) MAXDIAG * (long) (sizeof (struct dstruct)));
      exit (1);
   }

#ifdef TFAST
   if ((f_str->aa1x =(unsigned char *)calloc((size_t)ppst->maxlen+2,
					     sizeof(unsigned char)))
       == NULL) {
     fprintf (stderr, "*** error [%s:%d] - cannot allocate aa1x array %d\n",
	      __FILE__, __LINE__, ppst->maxlen+2);
     exit (1);
   }
   f_str->aa1x++;
#endif

   /* allocate space for the scoring arrays */
   maxn0 = n0 + 4;

   maxn0 = max(3*n0/2,MIN_RES);
   if ((res = (int *)calloc((size_t)maxn0,sizeof(int)))==NULL) {
     fprintf(stderr,"*** error [%s:%d] - cannot allocate alignment results array %d\n",
	     __FILE__, __LINE__, maxn0);
     exit(1);
   }
   f_str->res = res;
   f_str->max_res = maxn0;

   /* Tatusov Statistics Setup */

   /* initialize priors array. */
   if((f_str->priors = (double *)calloc(ppst->nsq+1, sizeof(double))) == NULL) {
     fprintf(stderr, "*** error [%s:%d] - cannot allocate priors array [%d]\n",
	     __FILE__, __LINE__, ppst->nsq+1);
     exit(1);
   }
   calc_priors(f_str->priors, ppst, f_str, NULL, 0, ppst->pseudocts);

   f_str->dotat = 0;
   f_str->shuff_cnt = ppst->shuff_node;

   /* End of Tatusov Statistics Setup */

   *f_arg = f_str;
}


/* pstring1 is a message to the manager, currently 512 */
/* pstring2 is the same information, but in a markx==10 format */
void
get_param (const struct pstruct *ppstr,
	   char **pstring1, char *pstring2,
	   struct score_count_s *s_info)
{
#ifndef TFAST
  char *pg_str="FASTF";
#else
  char *pg_str="TFASTF";
#endif

  sprintf (pstring1[0], "%s (%s)",pg_str,verstr);
  sprintf (pstring1[1], "%s matrix (%d:%d), join: %d",
	   ppstr->pam_name, ppstr->pam_h,ppstr->pam_l,ppstr->param_u.fa.cgap);

  if (ppstr->param_u.fa.iniflag) strcat(pstring1[0]," init1");

  if (pstring2 != NULL) {
    sprintf (pstring2, "; pg_name_alg: %s\n; pg_ver_rel: %s\n; pg_matrix: %s (%d:%d)\n\
; pg_join: %d\n",
	     pg_str,verstr, ppstr->pam_name, ppstr->pam_h,ppstr->pam_l,
	     ppstr->param_u.fa.cgap);
   }
}

void
close_work (const unsigned char *aa0, const int n0,
	    struct pstruct *ppst,
	    struct f_struct **f_arg)
{
  struct f_struct *f_str;
  struct savestr *vmptr;

  f_str = *f_arg;

  if (f_str != NULL) {

    for (vmptr = f_str->vmax; vmptr < &f_str->vmax[MAXSAV]; vmptr++)
      free(vmptr->used);

    free(f_str->res);
#ifdef TFAST
    free(f_str->aa1x - 1); /* allocated, then aa1x++'ed */
#endif
    free(f_str->diag);
    free(f_str->link);
    free(f_str->pamh2); 
    free(f_str->pamh1);
    free(f_str->harr);
    free(f_str->aa0t);
    free(f_str->aa0);
    free(f_str->priors);
    free(f_str->vmax);
    free(f_str->vptr);
    free(f_str);
    *f_arg = NULL;
  }
}

int do_fastf (unsigned char *aa0, int n0,
	      const unsigned char *aa1, int n1,
	      const struct pstruct *ppst, struct f_struct *f_str,
	      struct rstruct *rst, int *hoff, int opt_prob)
{
   int     nd;		/* diagonal array size */
   int     lhval;
   int     kfact;
   register struct dstruct *dptr;
   register int tscor;
   register struct dstruct *diagp;
   struct dstruct *dpmax;
   register int lpos;
   int     tpos, npos;
   struct savestr *vmptr;
   int     scor, tmp;
   int     im, ib, nsave;
   int     cmps ();		/* comparison routine for ksort */
   const int *hsq;

   hsq = ppst->hsq;

   if (n1 < 1) {
     rst->score[0] = rst->score[1] = rst->score[2] = 0;
     rst->escore = 1.0;
     rst->segnum = 0;
     rst->seglen = 0;
     return 1;
   }

   if (n0+n1+1 >= MAXDIAG) {
     fprintf(stderr,"*** error [%s:%d] - n0,n1 too large  %d +  %d > %d\n",
	     __FILE__, __LINE__, n0,n1, MAXDIAG);
     rst->score[0] = rst->score[1] = rst->score[2] = -1;
     rst->escore = 2.0;
     rst->segnum = 0;
     rst->seglen = 0;
     return -1;
   }

   nd = n0 + n1;

   dpmax = &f_str->diag[nd];
   for (dptr = &f_str->diag[f_str->ndo]; dptr < dpmax;) {
      dptr->stop = -1;
      dptr->dmax = NULL;
      dptr++->score = 0;
   }

   /* initialize the saved segment structures */
   for (vmptr = f_str->vmax; vmptr < &f_str->vmax[MAXSAV]; vmptr++) {
      vmptr->score = 0;
      memset(vmptr->used, 0, n0 * sizeof(int));
   }

   f_str->lowmax = f_str->vmax;
   f_str->lowscor = 0;

   /* start hashing */

   diagp = &f_str->diag[f_str->noff];
   for (lhval = lpos = 0; lpos < n1; lpos++, diagp++) {
     if (hsq[aa1[lpos]]>=NMAP) {
       lpos++ ; diagp++;
       while (lpos < n1 && hsq[aa1[lpos]]>=NMAP) {lpos++; diagp++;}
       if (lpos >= n1) break;
       lhval = 0;
     }
     lhval = hsq[aa1[lpos]];
     for (tpos = f_str->harr[lhval].pos, npos = f_str->harr[lhval].next;
	  tpos >= 0; tpos = f_str->link[npos].pos, npos = f_str->link[npos].next) {
       /* tscor gets position of end of current lpos diag run */
       if ((tscor = (dptr = &diagp[-tpos])->stop) >= 0) {
	 tscor++;		 /* move forward one */
	 if ((tscor -= lpos) <= 0) { /* check for size of gap to this hit - */
	                             /* includes implicit -1 mismatch penalty */
	   scor = dptr->score;	     /* current score of this run */
	   if ((tscor += (kfact = f_str->pamh2[lhval])) < 0 &&
	       f_str->lowscor < scor)	/* if updating tscor makes run worse, */
	     savemax (dptr, f_str);     /* save it */

	   if ((tscor += scor) >= kfact) {  /* add to current run if continuing */
	     				    /* is better than restart (kfact) */
	       dptr->score = tscor;
	       dptr->stop = lpos;
	     }
	     else {
	       dptr->score = kfact;	/* starting over is better */
	       dptr->start = (dptr->stop = lpos);
	     }
	 }
	 else {				/* continue current run */
	   dptr->score += f_str->pamh1[aa0[tpos]]; 
	   dptr->stop = lpos;
	 }
       }
       else {				/* no diagonal run yet */
	 dptr->score = f_str->pamh2[lhval];
	 dptr->start = (dptr->stop = lpos);
       }
     }				/* end tpos */
   }				/* end lpos */

   for (dptr = f_str->diag; dptr < dpmax;) {
     if (dptr->score > f_str->lowscor) savemax (dptr, f_str);
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

   /* set up pointers for sorting */

   for (nsave = 0, vmptr = f_str->vmax; vmptr < &f_str->vmax[MAXSAV]; vmptr++) {
     if (vmptr->score > 0) {
       vmptr->score = m0_spam (aa0, aa1, n1, vmptr, ppst->pam2[0], f_str);
       f_str->vptr[nsave++] = vmptr;
     }
   }

   /* sort them */
   kssort (f_str->vptr, nsave);

   
#ifdef DEBUG
   /*
   for (ib=0; ib<nsave; ib++) {
     fprintf(stderr,"0: %4d-%4d  1: %4d-%4d  dp: %d score: %d\n",
  	     f_str->noff+f_str->vptr[ib]->start-f_str->vptr[ib]->dp,
  	     f_str->noff+f_str->vptr[ib]->stop-f_str->vptr[ib]->dp,
  	     f_str->vptr[ib]->start,f_str->vptr[ib]->stop,
  	     f_str->vptr[ib]->dp,f_str->vptr[ib]->score);
       for (im=f_str->vptr[ib]->start; im<=f_str->vptr[ib]->stop; im++)
         fprintf(stderr," %c:%c",ppst->sq[aa0[f_str->noff+im-f_str->vptr[ib]->dp]],
		 ppst->sq[aa1[im]]);
       fputc('\n',stderr);
   }
   fprintf(stderr,"---\n");
   */
   /* now use m_spam to re-evaluate */
   /*
   for (tpos = 0; tpos < n0; tpos++) {
     fprintf(stderr,"%c:%2d ",ppst->sq[aa0[tpos]],aa0[tpos]);
     if (tpos %10 == 9) fputc('\n',stderr);
   }
   fputc('\n',stderr);
   */
#endif   
   
   f_str->aa0ix = 0;
   for (ib=0; ib < nsave; ib++) {
     if ((vmptr=f_str->vptr[ib])->score > 0) {
       vmptr->score = m1_spam (aa0, n0, aa1, n1, vmptr,
			       ppst->pam2[0], ppst->pam_l, f_str);
     }
   }
   /* reset aa0 - modified by m1_spam */
   for (tpos = 0; tpos < n0; tpos++) {
     if (aa0[tpos] >= 32) aa0[tpos] -= 32;
   }

   kssort(f_str->vptr,nsave);

   for ( ; nsave > 0; nsave--) 
     if (f_str->vptr[nsave-1]->score >0) break;

   if (nsave <= 0) {
     f_str->nsave = 0;
     rst->score[0] = rst->score[1] = rst->score[2] = 0;
     rst->escore = 1.0;

     return 1;
   }
   else f_str->nsave = nsave;

   
#ifdef DEBUG
   /*
   fprintf(stderr,"n0: %d; n1: %d; noff: %d\n",n0,n1,f_str->noff);
   for (ib=0; ib<nsave; ib++) {
     fprintf(stderr,"0: %4d-%4d  1: %4d-%4d  dp: %d score: %d\n",
  	     f_str->noff+f_str->vptr[ib]->start-f_str->vptr[ib]->dp,
  	     f_str->noff+f_str->vptr[ib]->stop-f_str->vptr[ib]->dp,
  	     f_str->vptr[ib]->start,f_str->vptr[ib]->stop,
  	     f_str->vptr[ib]->dp,f_str->vptr[ib]->score);
     for (im=f_str->vptr[ib]->start; im<=f_str->vptr[ib]->stop; im++)
       fprintf(stderr," %c:%c",ppst->sq[aa0[f_str->noff+im-f_str->vptr[ib]->dp]],
  	       ppst->sq[aa1[im]]);
     fputc('\n',stderr);
   }

   fprintf(stderr,"---\n");
   */
#endif   

   scor = sconn (f_str->vptr, nsave, ppst->param_u.fa.cgap, f_str,
		 rst, ppst, aa0, n0, aa1, n1, opt_prob);

   for (vmptr=f_str->vptr[0],ib=1; ib<nsave; ib++)
     if (f_str->vptr[ib]->score > vmptr->score) vmptr=f_str->vptr[ib];

   rst->score[1] = vmptr->score;
   rst->score[0] = rst->score[2] = max (scor, vmptr->score);

   return 1;
}

void do_work (const unsigned char *aa0, int n0,
	      const unsigned char *aa1, int n1,
	      int frame,
	      const struct pstruct *ppst, struct f_struct *f_str,
	      int qr_flg, int shuff_flg, struct rstruct *rst,
	      struct score_count_s *s_info)
{
  int opt_prob;
  int hoff, n10, i;

  if (qr_flg==1 && f_str->shuff_cnt <= 0) {
    rst->escore = 2.0;
    rst->score[0]=rst->score[1]=rst->score[2]= -1;
    rst->valid_stat = 0;
    return;
  }

  s_info->s_cnt[ppst->score_ix]++;
  s_info->tot_scores++;

  rst->valid_stat = 1;
  if (f_str->dotat || ppst->zsflag == 4 || ppst->zsflag == 14 ) opt_prob=1;
  else opt_prob = 0;
  if (ppst->zsflag == 2 || ppst->zsflag == 12) opt_prob = 0;
  if (qr_flg) {
    opt_prob=1;
    /*    if (frame==1) */
      f_str->shuff_cnt--;
  }

  if (n1 < 1) {
    rst->score[0] = rst->score[1] = rst->score[2] = -1;
    rst->escore = 2.0;
    return;
  }

#ifdef TFAST 
  n10=aatran(aa1,f_str->aa1x,n1,frame);
  if (ppst->debug_lib)
    for (i=0; i<n10; i++)
      if (f_str->aa1x[i]>ppst->nsq) {
	fprintf(stderr, "*** error [%s:%d] - residue[%d/%d] %d range (%d)\n",
		__FILE__, __LINE__, i,n1, f_str->aa1x[i],ppst->nsq);
	f_str->aa1x[i]=0;
	n10=i-1;
      }

  do_fastf (f_str->aa0, n0, f_str->aa1x, n10, ppst, f_str, rst, &hoff, opt_prob);
#else	/* FASTF */
  do_fastf (f_str->aa0, n0, aa1, n1, ppst, f_str, rst, &hoff, opt_prob);
#endif

  rst->comp = rst->H = -1.0;

}

void do_opt (const unsigned char *aa0, int n0,
	     const unsigned char *aa1, int n1,
	     int frame,
	     struct pstruct *ppst,
	     struct f_struct *f_str,
	     struct rstruct *rst)
{
  int optflag, tscore, hoff, n10;

  optflag = ppst->param_u.fa.optflag;
  ppst->param_u.fa.optflag = 1;

#ifdef TFAST  
  n10=aatran(aa1,f_str->aa1x,n1,frame);
  do_fastf (f_str->aa0, n0, f_str->aa1x, n10, ppst, f_str, rst, &hoff, 1);
#else	/* FASTA */
  do_fastf(f_str->aa0, n0, aa1, n1, ppst, f_str, rst, &hoff, 1);
#endif
  ppst->param_u.fa.optflag = optflag;
}

void
savemax (dptr, f_str)
  register struct dstruct *dptr;
  struct f_struct *f_str;
{
   register int dpos;
   register struct savestr *vmptr;
   register int i;

   dpos = (int) (dptr - f_str->diag);

/* check to see if this is the continuation of a run that is already saved */

   if ((vmptr = dptr->dmax) != NULL && vmptr->dp == dpos &&
	 vmptr->start == dptr->start)
   {
      vmptr->stop = dptr->stop;
      if ((i = dptr->score) <= vmptr->score)
	 return;
      vmptr->score = i;
      if (vmptr != f_str->lowmax)
	 return;
   }
   else
   {
      i = f_str->lowmax->score = dptr->score;
      f_str->lowmax->dp = dpos;
      f_str->lowmax->start = dptr->start;
      f_str->lowmax->stop = dptr->stop;
      dptr->dmax = f_str->lowmax;
   }

   for (vmptr = f_str->vmax; vmptr < &f_str->vmax[MAXSAV]; vmptr++)
      if (vmptr->score < i)
      {
	 i = vmptr->score;
	 f_str->lowmax = vmptr;
      }
   f_str->lowscor = i;
}

/* this version of spam() is designed to work with a collection of
   subfragments, selecting the best amino acid at each position so
   that, from each subfragment, each position is only used once.

   As a result, m_spam needs to know the number of fragments.

   In addition, it now requires a global alignment to the fragment
   and resets the start and stop positions

   */

static int
m1_spam (unsigned char *aa0, int n0,
	 const unsigned char *aa1, int n1,
	 struct savestr *dmax, int **pam2, int pam_l,
	 struct f_struct *f_str)
{
  int     tpos, lpos, im, ii, nm, ci;
  int     tot, ctot, pv;

  struct {
    int     start, stop, score;
  } curv, maxv;
  unsigned char *aa0p;
  const unsigned char *aa1p;

  lpos = dmax->start;                   /* position in library sequence */
   tpos = lpos - dmax->dp + f_str->noff; /* position in query sequence */
   /* force global alignment, reset start*/
   if (tpos < lpos) {
     lpos = dmax->start -= tpos;
     tpos = 0;
   }
   else {
     tpos -= lpos;
     lpos = dmax->start = 0;
   }

   dmax->stop = dmax->start + (f_str->nmoff -2 - tpos);
   if (dmax->stop > n1) dmax->stop = n1;

   /*
   if (dmax->start < 0) {
     tpos = -dmax->start;
     lpos = dmax->start=0;
   }
   else tpos = 0;
   */

   aa1p = &aa1[lpos];
   aa0p = &aa0[tpos];

   nm = f_str->nm0;

   tot = curv.score = maxv.score = 0;
   for (; lpos <= dmax->stop; lpos++,aa0p++,aa1p++) {
     ctot = pam_l;
     ci = -1;
     for (im = 0, ii=0; im < nm; im++,ii+=f_str->nmoff) {
       if (aa0p[ii] < 32 && (pv = pam2[aa0p[ii]][*aa1p]) > ctot) {
	 ctot = pv;
	 ci = ii;
/*  	 fprintf(stderr, "lpos: %d im: %d ii: %d ci: %d ctot: %d pi: %d pv: %d\n", lpos, im, ii, ci, ctot, aa0p[ii], pam2[aa0p[ii]][*aa1p]); */
       }
     }
     tot += ctot;
     if (ci >= 0 && aa0p[ci] < 32) {
#ifdef DEBUG
/*         fprintf(stderr, "used: lpos: %d ci: %d : %c\n", lpos, ci, sq[aa0p[ci]]); */
#endif
       aa0p[ci] +=  32;
       dmax->used[&aa0p[ci] - aa0] = 1;
     }
   }
   return tot;
}

int ma_spam (unsigned char *aa0, int n0, const unsigned char *aa1,
	     struct savestr *dmax, struct pstruct *ppst,
	     struct f_struct *f_str)
{
  int **pam2;
  int     tpos, lpos, im, ii, nm, ci, lp0;
  int     tot, ctot, pv;
  struct {
    int     start, stop, score;
  } curv, maxv;
   const unsigned char *aa1p;
   unsigned char *aa0p, *aa0pt;
   int aa0t_flg;

   pam2 = ppst->pam2[0];
   aa0t_flg = 0;

   lpos = dmax->start;			/* position in library sequence */
   tpos = lpos - dmax->dp + f_str->noff; /* position in query sequence */
   lp0 = lpos = dmax->start;
   aa1p = &aa1[lpos];
   aa0p = &aa0[tpos];			/* real aa0 sequence */

			/* the destination aa0 sequence (without nulls) */
   aa0pt = &f_str->aa0t[f_str->aa0ix];

   curv.start = lpos;
   nm = f_str->nm0;

   /* sometimes, tpos may be > 0, with lpos = 0 - fill with 'X' */
   if (lpos == 0 && tpos > 0)
     for (ii = 0; ii < tpos; ii++) *aa0pt++ = 31;  /* filler character */

   tot = curv.score = maxv.score = 0;
   for (; lpos <= dmax->stop; lpos++) {
     ctot = ppst->pam_l;
     ci = -1;
     for (im = 0, ii=0; im < nm; im++,ii+=f_str->nmoff) {
       if (aa0p[ii] < 32 && (pv = pam2[aa0p[ii]][*aa1p]) > ctot) {
	 ctot = pv;
	 ci = ii;
       }
     }
     tot += ctot;
     if (ci >= 0) {
       if (ci >= n0) {fprintf(stderr,"*** warning [%s:%d] - ci off end %d/%d\n",
			      __FILE__, __LINE__, ci,n0);}
       else {
	 *aa0pt++ = aa0p[ci];
	 aa0p[ci] +=  32;
	 aa0t_flg=1;
       }
     }
     aa0p++; aa1p++;
   }

   if (aa0t_flg) {
     dmax->dp -= f_str->aa0ix;		/* shift ->dp for aa0t */
     if ((ci=(int)(aa0pt-f_str->aa0t)) > n0) {
       fprintf(stderr,"*** warning [%s:%d] - aapt off %d/%d end\n",
	       __FILE__, __LINE__, ci,n0);
     }
     else 
       *aa0pt++ = 0;			/* skip over NULL */

     aa0pt = &f_str->aa0t[f_str->aa0ix];
     aa1p = &aa1[lp0];

     /*
     for (im = 0; im < f_str->nmoff; im++)
       fprintf(stderr,"%c:%c,",ppst->sq[aa0pt[im]],ppst->sq[aa1p[im]]);
     fprintf(stderr,"- %3d (%3d:%3d)\n",dmax->score,f_str->aa0ix,lp0);
     */

     f_str->aa0ix += f_str->nmoff;	/* update offset into aa0t */
   }
   /*
      fprintf(stderr," ma_spam returning: %d\n",tot);
   */
   return tot;
}

static int
m0_spam (unsigned char *aa0, const unsigned char *aa1, int n1,
	 struct savestr *dmax, int **pam2,
	 struct f_struct *f_str)
{
   int tpos, lpos, lend, im, ii, nm;
   int     tot, ctot, pv;
   struct {
     int     start, stop, score;
   } curv, maxv;
   const unsigned char *aa0p, *aa1p;

   lpos = dmax->start;			/* position in library sequence */
   tpos = lpos - dmax->dp + f_str->noff; /* position in query sequence */
   if (tpos > 0) {
     if (lpos-tpos >= 0) {
       lpos = dmax->start -= tpos;    /* force global alignment, reset start*/
       tpos = 0;
     }
     else {
       tpos -= lpos;
       lpos = dmax->start = 0;
     }
   }

   nm = f_str->nm0;
   lend = dmax->stop;
   if (n1 - (lpos + f_str->nmoff-2) < 0 ) {
     lend = dmax->stop = (lpos - tpos) + f_str->nmoff-2;
     if (lend >= n1) lend = n1-1;
   }

   aa1p = &aa1[lpos];
   aa0p = &aa0[tpos];

   curv.start = lpos;

   tot = curv.score = maxv.score = 0;
   for (; lpos <= lend; lpos++) {
     ctot = -10000;
     for (im = 0, ii=0; im < nm; im++,ii+=f_str->nmoff) {
       if ((pv = pam2[aa0p[ii]][*aa1p]) > ctot) {
	 ctot = pv;
       }
     }
     tot += ctot;
     aa0p++; aa1p++;
   }

   /* reset dmax if necessary */

   return tot;
}

/* sconn links up non-overlapping alignments and calculates the score */

int sconn (struct savestr **v, int n, int cgap, struct f_struct *f_str,
	   struct rstruct *rst, const struct pstruct *ppst,
	   const unsigned char *aa0, int n0,
	   const unsigned char *aa1, int n1,
	   int opt_prob)
{
   int     i, si, cmpp ();
   struct slink *start, *sl, *sj, *so, sarr[MAXSAV];
   int     lstart, plstop;
   double tatprob;

   /* sarr[] saves each alignment score/position, and provides a link
      back to the previous alignment that maximizes the score */

   /*	sort the score left to right in lib pos */
   kpsort (v, n);

   start = NULL;

   /* for the remaining runs, see if they fit */
   for (i = 0, si = 0; i < n; i++) {

     /* if the score is less than the gap penalty, it never helps */
     if (!opt_prob && (v[i]->score < cgap) ){ continue; }

     lstart = v[i]->start;

     /* put the run in the group */
     sarr[si].vp = v[i];
     sarr[si].score = v[i]->score;
     sarr[si].next = NULL;
     sarr[si].prev = NULL;
     sarr[si].tat = NULL;
     
     if(opt_prob) {
       sarr[si].tatprob = 
	 calc_tatusov(NULL, &sarr[si], aa0, n0, aa1, n1,
		      ppst->pam2[0],ppst->nsq, f_str,
		      ppst->pseudocts, opt_prob,ppst->zsflag);
       sarr[si].tat = sarr[si].newtat;
     }

     /* if it fits, then increase the score */
     for (sl = start; sl != NULL; sl = sl->next) {
       plstop = sl->vp->stop;
       /* if end < start or start > end, add score */
       if (plstop < lstart ) {
	 if(!opt_prob) {
	   sarr[si].score = sl->score + v[i]->score;
	   sarr[si].prev = sl;
	   /*
	     fprintf(stderr,"sconn %d added %d/%d getting %d; si: %d, tat: %g\n",
	     i,v[i]->start, v[i]->score,sarr[si].score,si, 2.0);
	   */
	   break;
	 } else {
	   tatprob = 
	     calc_tatusov(sl, &sarr[si], aa0, n0, aa1, n1,
			  ppst->pam2[0], ppst->nsq, f_str,
			  ppst->pseudocts, opt_prob, ppst->zsflag);
	   /* if our tatprob gets worse when we add this, forget it */
	   if(tatprob > sarr[si].tatprob) {
	     free(sarr[si].newtat->probs); /* get rid of new tat struct */
	     free(sarr[si].newtat);
	     continue;
	   } else {
	     sarr[si].tatprob = tatprob;
	     free(sarr[si].tat->probs); /* get rid of old tat struct */
	     free(sarr[si].tat);
	     sarr[si].tat = sarr[si].newtat;
	     sarr[si].prev = sl;
	     sarr[si].score = sl->score + v[i]->score;
	     /*
	       fprintf(stderr,"sconn TAT %d added %d/%d getting %d; si: %d, tat: %g\n",
	       i,v[i]->start, v[i]->score,sarr[si].score,si, tatprob);
	     */
	     break;
	   }
	 }
       }
     }

     /*	now recalculate where the score fits - resort the scores */
     if (start == NULL) {
       start = &sarr[si];
     } else {
       if(!opt_prob) { /* sort by scores */
	 for (sj = start, so = NULL; sj != NULL; sj = sj->next) {
	   if (sarr[si].score > sj->score) { /* if new score > best score */
	     sarr[si].next = sj;	     /* previous best linked to best */
	     if (so != NULL)		
	       so->next = &sarr[si];	     /* old best points to new best */
	     else
	       start = &sarr[si];
	     break;
	   }
	   so = sj;			     /* old-best saved in so */
	 }
       } else { /* sort by tatprobs */
	 for (sj = start, so = NULL; sj != NULL; sj = sj->next) {
	   if ( sarr[si].tatprob < sj->tatprob ||
		((sarr[si].tatprob == sj->tatprob) && sarr[si].score > sj->score) ) {
	     sarr[si].next = sj;
	     if (so != NULL)
	       so->next = &sarr[si];
	     else
	       start = &sarr[si];
	     break;
	   }
	   so = sj;
	 }
       }
     }
     si++;
   }
   
   if(opt_prob) {
     for (i = 0 ; i < si ; i++) {
       free(sarr[i].tat->probs);
       free(sarr[i].tat);
     }
   }

   if (start != NULL) {

     if(opt_prob)
       rst->escore = start->tatprob;
     else
       rst->escore = 2.0;

     rst->segnum = rst->seglen = 0;
     for(sj = start ; sj != NULL; sj = sj->prev) {
       rst->segnum++;
       rst->seglen += sj->vp->stop - sj->vp->start + 1;
     }
     return (start->score);
   } else {

     if(opt_prob)
       rst->escore = 1.0;
     else
       rst->escore = 2.0;

     rst->segnum = rst->seglen = 0;
     return (0);
   }
}

void
kssort (struct savestr **v, int n)
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
kpsort (v, n)
struct savestr *v[];
int     n;
{
   int     gap, i, j;
   struct savestr *tmp;

   for (gap = n / 2; gap > 0; gap /= 2)
      for (i = gap; i < n; i++)
	 for (j = i - gap; j >= 0; j -= gap)
	 {
	    if (v[j]->start <= v[j + gap]->start)
	       break;
	    tmp = v[j];
	    v[j] = v[j + gap];
	    v[j + gap] = tmp;
	 }
}

/* sorts alignments from right to left (back to front) based on stop */

void
krsort (v, n)
struct savestr *v[];
int     n;
{
   int     gap, i, j;
   struct savestr *tmp;

   for (gap = n / 2; gap > 0; gap /= 2)
      for (i = gap; i < n; i++)
	 for (j = i - gap; j >= 0; j -= gap)
	 {
	    if (v[j]->stop > v[j + gap]->stop)
	       break;
	    tmp = v[j];
	    v[j] = v[j + gap];
	    v[j + gap] = tmp;
	 }
}

struct a_res_str *
do_walign (const unsigned char *aa0, int n0,
	   const unsigned char *aa1, int n1,
	   int frame, int repeat_thresh,
	   struct pstruct *ppst, 
	   struct f_struct *f_str, 
	   int *have_ares)
{
  struct a_res_str *a_res;
  int hoff, n10;
  int ib;
  unsigned char *aa0t;
  const unsigned char *aa1p;

  *have_ares = 0x2;	/* set 0x2 bit to indicate local copy */

  if ((a_res = (struct a_res_str *)calloc(1, sizeof(struct a_res_str)))==NULL) {
    fprintf(stderr,"*** error [%s:%d] - cannot allocate a_res [%lu]",
	    __FILE__, __LINE__, sizeof(struct a_res_str));
    return NULL;
  }

#ifdef TFAST
  f_str->n10 = n10 = aatran(aa1,f_str->aa1x,n1,frame);
  aa1p = f_str->aa1x;
#else
  n10 = n1;
  aa1p = aa1;
#endif

  do_fastf(f_str->aa0, n0, aa1p, n10, ppst, f_str, &a_res->rst, &hoff, 1);

  /* the alignment portion takes advantage of the information left
     over in f_str after do_fastf is done.  in particular, it is
     easy to run a modified sconn() to produce the alignments.

     unfortunately, the alignment display routine wants to have
     things encoded as with bd_align and sw_align, so we need to do that.
     */

  if ((aa0t = (unsigned char *)calloc(n0+1,sizeof(unsigned char)))==NULL) {
    fprintf(stderr,"*** error [%s:%d] - cannot allocate aa0t %d\n",
	    __FILE__, __LINE__, n0+1);
    exit(1);
  }

   kssort (f_str->vptr, f_str->nsave);
   f_str->aa0ix = 0;
   if (f_str->nsave > f_str->nm0) f_str->nsave = f_str->nm0;
   for (ib=0; ib < f_str->nm0; ib++) {
     if (f_str->vptr[ib]->score > 0) {
       f_str->vptr[ib]->score = 
	 ma_spam (f_str->aa0, n0, aa1p, f_str->vptr[ib], ppst, f_str);
     }
   }

   /* after ma_spam is over, we need to reset aa0 */
   for (ib = 0; ib < n0; ib++) {
     if (f_str->aa0[ib] >= 32) f_str->aa0[ib] -= 32;
   }

   kssort(f_str->vptr,f_str->nsave);

   for ( ; f_str->nsave > 0; f_str->nsave--) 
     if (f_str->vptr[f_str->nsave-1]->score >0) break;

  a_res->nres = sconn_a (aa0t,n0, ppst->param_u.fa.cgap, f_str,a_res);
  free(aa0t);

  a_res->res = f_str->res;
  a_res->sw_score = a_res->rst.score[0];
  return a_res;
}

/* this version of sconn is modified to provide alignment information */

int sconn_a (unsigned char *aa0, int n0, int cgap, 
	     struct f_struct *f_str,
	     struct a_res_str *a_res)
{
   int     i, si, cmpp (), n;
   unsigned char *aa0p;
   int sx, dx, doff;

   struct savestr **v;
   struct slink {
     int     score;
     struct savestr *vp;
     struct slink *snext;
     struct slink *aprev;
   } *start, *sl, *sj, *so, sarr[MAXSAV];
   int     lstop, plstart;
   int *res, nres, tres;

/*	sort the score left to right in lib pos */

   v = f_str->vptr;
   n = f_str->nsave;

   krsort (v, n);	/* sort from left to right in library */

   start = NULL;

/*	for each alignment, see if it fits */

   for (i = 0, si = 0; i < n; i++) {

/*	if the score is less than the join threshold, skip it */
     if (v[i]->score < cgap) continue;

     lstop = v[i]->stop;		/* have right-most lstart */

/*	put the alignment in the group */

     sarr[si].vp = v[i];
     sarr[si].score = v[i]->score;
     sarr[si].snext = NULL;
     sarr[si].aprev = NULL;

/* 	if it fits, then increase the score */
/* start points to a sorted (by total score) list of candidate
   overlaps */

     for (sl = start; sl != NULL; sl = sl->snext) { 
       plstart = sl->vp->start;
       if (plstart > lstop ) {
	 sarr[si].score = sl->score + v[i]->score;
	 sarr[si].aprev = sl;
	 break;		/* quit as soon as the alignment has been added */
       }
     }

/* now recalculate the list of best scores */
     if (start == NULL)
       start = &sarr[si];	/* put the first one in the list */
     else
       for (sj = start, so = NULL; sj != NULL; sj = sj->snext) {
	 if (sarr[si].score > sj->score) { /* new score better than old */
	   sarr[si].snext = sj;		/* snext best after new score */
	   if (so != NULL)
	     so->snext = &sarr[si];	/* prev_best->snext points to best */
	   else  start = &sarr[si];	/* start points to best */
	   break;			/* stop looking */
	 }
	 so = sj;		/* previous candidate best */
       }
     si++;				/* increment to snext alignment */
   }

   /* we have the best set of alignments, write them to *res */
   if (start != NULL) {
     res = f_str->res;	/* set a destination for the alignment ops */
     tres = nres = 0;	/* alignment op length = 0 */
     aa0p = aa0;	/* point into query (needed for calcons later) */
     a_res->min1 = start->vp->start;	/* start in library */
     a_res->min0 = 0;			/* start in query */
     for (sj = start; sj != NULL; sj = sj->aprev ) {
       doff = (int)(aa0p-aa0) - (sj->vp->start-sj->vp->dp+f_str->noff);
       /*
       fprintf(stderr,"doff: %3d\n",doff);
       */
       for (dx=sj->vp->start,sx=sj->vp->start-sj->vp->dp+f_str->noff;
	    dx <= sj->vp->stop; dx++) {
	 *aa0p++ = f_str->aa0t[sx++];	/* copy residue into aa0 */
	 tres++;			/* bump alignment counter */
	 res[nres++] = 0;		/* put 0-op in res */
       }
       sj->vp->dp -= doff;
       if (sj->aprev != NULL) {
	 if (sj->aprev->vp->start - sj->vp->stop - 1 > 0 )
	 /* put an insert op into res to get to next aligned block */
	   tres += res[nres++] = (sj->aprev->vp->start - sj->vp->stop - 1);
       }
       /*
       fprintf(stderr,"t0: %3d, tx: %3d, l0: %3d, lx: %3d, dp: %3d noff: %3d, score: %3d\n",
	       sj->vp->start - sj->vp->dp + f_str->noff,
	       sj->vp->stop - sj->vp->dp + f_str->noff,
	       sj->vp->start,sj->vp->stop,sj->vp->dp,
	       f_str->noff,sj->vp->score);
       fprintf(stderr,"%3d - %3d: %3d\n",
	       sj->vp->start,sj->vp->stop,sj->vp->score);
       */
       a_res->max1 = sj->vp->stop;
       a_res->max0 = a_res->max1 - sj->vp->dp + f_str->noff;
     }

     /*
     fprintf(stderr,"(%3d - %3d):(%3d - %3d)\n",
     a_res->min0,a_res->max0,a_res->min1,a_res->max1);
     */

     /* now replace f_str->aa0t with aa0 */
     for (i=0; i<n0; i++) f_str->aa0t[i] = aa0[i];

     return tres;
   }
   else return (0);
}

/* calculate the 100% identical score */
int
shscore(unsigned char *aa0, int n0, int **pam2, int nsq)
{
  int i, sum;
  for (i=0,sum=0; i<n0; i++)
    if (aa0[i]!=0 && aa0[i]<=nsq) sum += pam2[aa0[i]][aa0[i]];
  return sum;
}

void
pre_cons(const unsigned char *aa1, int n1, int frame, struct f_struct *f_str) {

#ifdef TFAST
  f_str->n10=aatran(aa1,f_str->aa1x,n1,frame);
#endif
}

/* aln_func_vals - set up aln.qlfact, qlrev, llfact, llmult, frame, llrev */
/* call from calcons, calc_id, calc_code */
void 
aln_func_vals(int frame, struct a_struct *aln) {

#ifdef TFAST
  aln->qlrev = 0;
  aln->qlfact = 1;
  aln->llfact = aln->llmult = 3;
  aln->frame = 0;
  if (frame > 3) aln->llrev = 1;
#else	/* FASTF */
  aln->llfact = aln->qlfact = aln->llmult = 1;
  aln->llrev = aln->qlrev = 0;
  aln->frame = 0;
#endif
}

void aa0shuffle(unsigned char *aa0, int n0, struct f_struct *f_str) {

  int i, j, k;
  unsigned char tmp;

  for (i = f_str->nmoff-1 ; --i ; ) {

    /* j = nrand(i); if (i == j) continue;*/       /* shuffle columns */ 
    j = (f_str->nmoff - 2) - i; if (i <= j) break; /* reverse columns */

    /* swap all i'th column residues for all j'th column residues */
    for(k = 0 ; k < f_str->nm0 ; k++) {
      tmp = aa0[(k * (f_str->nmoff)) + i];
      aa0[(k * (f_str->nmoff)) + i] = aa0[(k * (f_str->nmoff)) + j];
      aa0[(k * (f_str->nmoff)) + j] = tmp;
    }
  }
}
