/*   $Id: dropfs2.c 1254 2014-01-29 16:03:40Z wrp $ */
/* $Revision: 1254 $  */

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

/* changed to return 2.0, rather than -1.0, for failure */

/* Feb 4, 2005 - modifications to allow searches with ktup=2 for very
   long queries.  This is a temporary solution to savemax(), spam()
   which do not preserve exact matches

   do_fasts() has been modified to allow higher maxsav for do_walign
   than for do_work (2*nsegs, 6*nsegs)
 */

/* this code implements the "fasts" algorithm, which compares a set of
   protein fragments to a protein sequence.  Comma's are used to separate
   the sequence fragments, which need not be the same length.

   The expected input is:

   >mgstm1
   MGDAPDFD,
   MILGYW,
   MLLEYTDS

   The fragments do not need to be in the correct order (which is
   presumably unknown from the peptide sequencing.
*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <math.h>

#include "defs.h"
#include "param.h"
#include "tatstats.h"

#define EOSEQ 0
#define ESS 59
#define NMAP_X 21 /* for 'X' - changed for NCBI */
#define NMAP_Z 25 /* for '*' - changed for NCBI */
#define MAXHASH 32
#define NMAP MAXHASH+1

static char *verstr="4.32 Feb 2007";

#define DROP_INTERN
#include "drop_func.h"

int shscore(const unsigned char *aa0, const int n0, int **pam2, int nsq);
extern void aancpy(char *to, char *from, int count, struct pstruct *ppst);

#ifdef TFAST
extern int aatran(const unsigned char *ntseq, unsigned char *aaseq, const int maxs, const int frame);
#endif

void savemax(struct dstruct *, struct f_struct *, int maxsav, int exact,int t_end);

int spam(const unsigned char *, const unsigned char *, int, struct savestr *, int **, struct f_struct *);
int sconn(struct savestr **v,
	  int nsave,
	  struct f_struct *,
	  struct rstruct *,
	  const struct pstruct *,
	  const unsigned char *aa0, int n0,
	  const unsigned char *aa1, int n1,
	  int opt_prob);

void kpsort(struct savestr **, int);
void kssort(struct savestr **, int);	/* sort by score */
int sconn_a(unsigned char *, int, 
	    const unsigned char *, int,
	    struct f_struct *, 
	    struct a_res_str *,
	    struct pstruct *);
void kpsort(struct savestr **, int);

/* initialize for fasta */

void
init_work (unsigned char *aa0, const int n0, 
	   struct pstruct *ppst,
	   struct f_struct **f_arg
	   )
{
   int mhv, phv;
   int hmax, nsegs;
   int i0, ib, hv, old_hv;
   int pamfact;
   struct f_struct *f_str;
   /* these used to be globals, but do not need to be */
   int ktup, fact, kt1;

   int maxn0;
   int stmp;	/* temporary score */
   int tmp_zsflag;
   int i, j, q;
   int tat_size;
   int *res;

   unsigned char *query;
   int k, l, m, n, N, length, index;

   double *tatprobptr;

   f_str = (struct f_struct *)calloc(1,sizeof(struct f_struct));

   ppst->param_u.fa.pgap = ppst->gdelval + ppst->ggapval;
   ktup = ppst->param_u.fa.ktup;
   if ( ktup  > ppst->param_u.fa.bktup ) {
     ktup = ppst->param_u.fa.ktup = ppst->param_u.fa.bktup;
   }
   fact = ppst->param_u.fa.scfact;

   /* fasts3 cannot work with lowercase symbols as low complexity;
      thus, NMAP must be disabled; this depends on aascii['X']  */
   if (ppst->hsq[NMAP_X] == NMAP ) {ppst->hsq[NMAP_X]=1;}
   if (ppst->hsq[NMAP_Z] == NMAP ) {ppst->hsq[NMAP_Z]=1;}
   /* this does not work in a threaded environment */
   /*    else {fprintf(stderr," cannot find 'X'==NMAP\n");} */

   for (i0 = 1, mhv = -1; i0 < ppst->nsq; i0++)
      if (ppst->hsq[i0] < NMAP && ppst->hsq[i0] > mhv)  mhv = ppst->hsq[i0];

   if (mhv <= 0) {
      fprintf (stderr, " maximum hsq <=0 %d\n", mhv);
      exit (1);
   }

   for (f_str->kshft = 0; mhv > 0; mhv /= 2) f_str->kshft++;

/*      kshft = 2;	*/
   kt1 = ktup-1;
   hv = 1;
   for (i0 = 0; i0 < ktup; i0++) hv = hv << f_str->kshft;
   hmax = hv;
   f_str->hmask = (hmax >> f_str->kshft) - 1;

   if ((f_str->aa0t = (unsigned char *) calloc(n0+1, sizeof(char))) == NULL) {
     fprintf (stderr, " cannot allocate f_str0->aa0t array; %d\n",n0+1);
     exit (1);
   }

   if ((f_str->aa0ti = (int *) calloc(n0+1, sizeof(int))) == NULL) {
     fprintf (stderr, " cannot allocate f_str0->aa0ti array; %d\n",n0+1);
     exit (1);
   }

   if ((f_str->aa0b = (int *) calloc(n0+1, sizeof(int))) == NULL) {
     fprintf (stderr, " cannot allocate f_str0->aa0b array; %d\n",n0+1);
     exit (1);
   }

   if ((f_str->aa0e = (int *) calloc(n0+1, sizeof(int))) == NULL) {
     fprintf (stderr, " cannot allocate f_str0->aa0e array; %d\n",n0+1);
     exit (1);
   }

   if ((f_str->aa0i = (int *) calloc(n0+1, sizeof(int))) == NULL) {
     fprintf (stderr, " cannot allocate f_str0->aa0i array; %d\n",n0+1);
     exit (1);
   }

   if ((f_str->aa0s = (int *) calloc(n0+1, sizeof(int))) == NULL) {
     fprintf (stderr, " cannot allocate f_str0->aa0s array; %d\n",n0+1);
     exit (1);
   }

   if ((f_str->aa0l = (int *) calloc(n0+1, sizeof(int))) == NULL) {
     fprintf (stderr, " cannot allocate f_str0->aa0l array; %d\n",n0+1);
     exit (1);
   }

   if ((f_str->harr = (int *) calloc (hmax, sizeof (int))) == NULL) {
     fprintf (stderr, " cannot allocate hash array: hmax: %d hmask: %d\n",
	      hmax, f_str->hmask);
     exit (1);
   }
   if ((f_str->pamh1 = (int *) calloc (ppst->nsq+1, sizeof (int))) == NULL) {
     fprintf (stderr, " cannot allocate pamh1 array\n");
     exit (1);
   }
   if ((f_str->pamh2 = (int *) calloc (hmax, sizeof (int))) == NULL) {
     fprintf (stderr, " cannot allocate pamh2 array\n");
     exit (1);
   }

   if ((f_str->link = (int *) calloc (n0, sizeof (int))) == NULL) {
     fprintf (stderr, " cannot allocate hash link array");
     exit (1);
   }

   /* for FASTS/FASTM, we want to know when we get to the end of a peptide,
      so we can ensure that we set the end and restart */

   if ((f_str->l_end = (int *) calloc (n0, sizeof (int))) == NULL) {
     fprintf (stderr, " cannot allocate link end array");
     exit (1);
   }

   for (i0 = 0; i0 < hmax; i0++) f_str->harr[i0] = -1;
   for (i0 = 0; i0 < n0; i0++) f_str->link[i0] = -1;
   for (i0 = 0; i0 < n0; i0++) f_str->l_end[i0] = 0;

   /* count the number of peptides */
   nsegs = 1;
   for (i0 = 0; i0 < n0; i0++) {
     if (aa0[i0] == ESS || aa0[i0] == 0) nsegs++;
   }

   /* allocate space for peptides offsets, nm_u */
   if ((f_str->nmoff = (int *)calloc(nsegs+1, sizeof(int)))==NULL) {
     fprintf(stderr, " cannot allocat nmoff array: %d\n", nsegs);
     exit(1);
   }

   if ((f_str->nm_u = (int *)calloc(nsegs+1, sizeof(int)))==NULL) {
     fprintf(stderr, " cannot allocat nm_u array: %d\n", nsegs);
     exit(1);
   }

   phv = hv = 0;
   f_str->nmoff[0] = 0;
   f_str->nm0 = 1;

   /* encode the aa0 array */
   if (kt1 > 0) {
     hv = ppst->hsq[aa0[0]];
     phv = ppst->pam2[0][aa0[0]][aa0[0]];
   }

   for (i0=kt1 ; i0 < n0; i0++) {
     if (aa0[i0] == ESS || aa0[i0] == 0) {
       /*       fprintf(stderr," converted %d to 0\n",aa0[i0]); */
       aa0[i0] = EOSEQ;	/* set ESS to 0 */
       f_str->nmoff[f_str->nm0++] = i0+1; 
       f_str->l_end[i0-1] = 1;
       phv = hv = 0;
       if (kt1 > 0) {
	 i0++;
	 hv = ppst->hsq[aa0[i0]];
	 phv = ppst->pam2[0][aa0[i0]][aa0[i0]];
       }
       continue;
     }

     hv = ((hv & f_str->hmask) << f_str->kshft) + ppst->hsq[aa0[i0]];
     f_str->link[i0] = f_str->harr[hv];
     f_str->harr[hv] = i0;
     f_str->pamh2[hv] = (phv += ppst->pam2[0][aa0[i0]][aa0[i0]]);
     phv -= ppst->pam2[0][aa0[i0 - kt1]][aa0[i0 - kt1]];
   }
   f_str->l_end[n0-1] = 1;

   f_str->nmoff[f_str->nm0] = n0+1;

   /*
#ifdef DEBUG
   fprintf(stderr, ">>%s\n",qtitle);
   for (j=0; j<f_str->nm0; j++) {
     for (i=f_str->nmoff[j]; i < f_str->nmoff[j+1]-1; i++) {
       fprintf(stderr,"%c",ppst->sq[aa0[i]]);
     }
     fprintf(stderr," %d\n",aa0[i]);
   }

   for (j=1; j<=ppst->nsq; j++) {
     fprintf(stderr, "%c %d\n", ppst->sq[j], f_str->harr[j]);
   }

   for (j=0; j<=n0; j++) {
     fprintf(stderr, "%c %d\n", ppst->sq[aa0[j]], f_str->link[j]);
   }

#endif
   */

   /* build an integer array of the max score that can be achieved
      from that position - use in savemax to mark some segments as
      fixed */

   /* setup aa0b[], aa0e[], which specify the begining and end of each
      segment */

   stmp = 0;
   q = -1;
   for (ib = i0 = 0; i0 < n0; i0++) {
     f_str->aa0l[i0] = i0 - q;
     if (aa0[i0]==EOSEQ) {
       f_str->aa0b[i0] = -1;
       f_str->aa0e[i0] = -1;
       f_str->aa0i[i0] = -1;
       f_str->aa0l[i0] = -1;
       q = i0;
       if (i0 > 0)f_str->aa0s[i0-1] = stmp;
       stmp = 0;
       ib++;
     }
     else {
       stmp += ppst->pam2[0][aa0[i0]][aa0[i0]];
     }

     f_str->aa0b[i0] =  f_str->nmoff[ib];
     f_str->aa0e[i0] =  f_str->nmoff[ib+1]-2;
     f_str->aa0i[i0] =  ib;

     /*
     fprintf(stderr,"%2d %c: %2d %2d %2d\n",i0,ppst->sq[aa0[i0]],
	     f_str->aa0b[i0],f_str->aa0e[i0],f_str->aa0i[i0]);
     */
   }
   f_str->aa0s[n0-1]=stmp;	/* save last best possible score */

   /* maxsav - maximum number of peptide alignments saved in search */
   /* maxsav_w - maximum number of peptide alignments saved in
      alignment */

   f_str->maxsav = max(MAXSAV,2*f_str->nm0);
   f_str->maxsav_w = max(MAXSAV,6*f_str->nm0);

   if ((f_str->vmax = (struct savestr *)
	calloc(f_str->maxsav_w,sizeof(struct savestr)))==NULL) {
     fprintf(stderr, "Couldn't allocate vmax[%d].\n",f_str->maxsav_w);
     exit(1);
   }

   if ((f_str->vptr = (struct savestr **)
	calloc(f_str->maxsav_w,sizeof(struct savestr *)))==NULL) {
     fprintf(stderr, "Couldn't allocate vptr[%d].\n",f_str->maxsav_w);
     exit(1);
   }

   if ((f_str->sarr = (struct slink *)
	calloc(f_str->maxsav_w,sizeof(struct slink)))==NULL) {
     fprintf(stderr, "Couldn't allocate sarr[%d].\n",f_str->maxsav_w);
     exit(1);
   }

   /* Tatusov Statistics Setup */

   /* initialize priors array. */
   if((f_str->priors = (double *)calloc(ppst->nsq+1, sizeof(double))) == NULL) {
     fprintf(stderr, "Couldn't allocate priors array.\n");
     exit(1);
   }

   calc_priors(f_str->priors, ppst, f_str, NULL, 0, ppst->pseudocts);

   /* pre-calculate the Tatusov probability array for each full segment */

   tmp_zsflag = ppst->zsflag;
   if (tmp_zsflag > 20) tmp_zsflag -= 20;
   if (tmp_zsflag > 10) tmp_zsflag -= 10;
   if (tmp_zsflag >= 1 && tmp_zsflag <= 3 && f_str->nm0 <= 10) {

     tat_size = (1<<f_str->nm0) -1;
     f_str->dotat = 1;
     f_str->tatprobs = (struct tat_str **) malloc((size_t)tat_size*sizeof(struct tat_str *));
     if (f_str->tatprobs == NULL) {
       fprintf (stderr, " cannot allocate tatprobs array: %ld\n",
		tat_size * sizeof(struct tat_str *));
       exit (1);
     }

     f_str->intprobs = (double **) malloc((size_t)tat_size * sizeof(double *));
     if(f_str->intprobs == NULL) {
       fprintf(stderr, "Couldn't allocate intprobs array.\n");
       exit(1);
     }

     for(k = 0, l = f_str->nm0 ; k < l ; k++) {
       query = &(aa0[f_str->nmoff[k]]);
       length = f_str->nmoff[k+1] - f_str->nmoff[k] - 1;

       /* this segment alone */
       index = (1 << k) - 1;
       generate_tatprobs(query, 0, length - 1, f_str->priors, ppst->pam2[0], ppst->nsq, &(f_str->tatprobs[index]), NULL);

       /* integrate the probabilities */
       N = f_str->tatprobs[index]->highscore - f_str->tatprobs[index]->lowscore;
       tatprobptr = (double *) calloc(N+1, sizeof(double));
       if(tatprobptr == NULL) {
	 fprintf(stderr, "Couldn't calloc tatprobptr.\n");
	 exit(1);
       }
       f_str->intprobs[index] = tatprobptr;

       for (i = 0; i <= N ; i++ ) {
	 tatprobptr[i] = f_str->tatprobs[index]->probs[i];
	 for (j = i + 1 ; j <= N ; j++ ) {
	   tatprobptr[i] += f_str->tatprobs[index]->probs[j];
	 }
       }

       /* this segment built on top of all other subcombinations */
       for(i = 0, j = (1 << k) - 1 ; i < j ; i++) {
	 index = (1 << k) + i;
	 generate_tatprobs(query, 0, length - 1, f_str->priors, ppst->pam2[0], ppst->nsq, &(f_str->tatprobs[index]), f_str->tatprobs[i]);

	 /* integrate the probabilities */
	 N = f_str->tatprobs[index]->highscore - f_str->tatprobs[index]->lowscore;
	 tatprobptr = (double *) calloc(N+1, sizeof(double));
	 if(tatprobptr == NULL) {
	   fprintf(stderr, "Couldn't calloc tatprobptr.\n");
	   exit(1);
	 }
	 f_str->intprobs[index] = tatprobptr;
       
	 for (m = 0; m <= N ; m++ ) {
	   tatprobptr[m] = f_str->tatprobs[index]->probs[m];
	   for (n = m + 1 ; n <= N ; n++ ) {
	     tatprobptr[m] += f_str->tatprobs[index]->probs[n];
	   }
	 }
       }
     }
   } else {
     f_str->dotat = 0;
     f_str->shuff_cnt = ppst->shuff_node;
   }

   /* End of Tatusov Statistics Setup */

   /*
   for (i0=1; i0<=ppst->nsq; i0++) {
     fprintf(stderr," %c: %2d ",ppst->sq[i0],f_str->harr[i0]);
     hv = f_str->harr[i0];
     while (hv >= 0) {
       fprintf(stderr," %2d",f_str->link[hv]);
       hv = f_str->link[hv];
     }
     fprintf(stderr,"\n");
   }
   */

/* this has been modified from 0..<ppst->nsq to 1..<=ppst->nsq because the
   pam2[0][0] is now undefined for consistency with blast
*/
   for (i0 = 1; i0 <= ppst->nsq; i0++)
     f_str->pamh1[i0] = ppst->pam2[0][i0][i0];

   f_str->ndo = 0;
   f_str->noff = n0-1;
   if (f_str->diag==NULL) 
     f_str->diag = (struct dstruct *) calloc ((size_t)MAXDIAG,
					      sizeof (struct dstruct));
   if (f_str->diag == NULL) {
      fprintf (stderr, " cannot allocate diagonal arrays: %ld\n",
	      (long) MAXDIAG * (long) (sizeof (struct dstruct)));
      exit (1);
   }

#ifdef TFAST
   if ((f_str->aa1x =(unsigned char *)calloc((size_t)ppst->maxlen+2,
					     sizeof(unsigned char)))
       == NULL) {
     fprintf (stderr, "cannot allocate aa1x array %d\n", ppst->maxlen+2);
     exit (1);
   }
   f_str->aa1x++;
#endif

   maxn0 = max(3*n0/2,MIN_RES);
   if ((res = (int *)calloc((size_t)maxn0,sizeof(int)))==NULL) {
     fprintf(stderr,"cannot allocate alignment results array %d\n",maxn0);
     exit(1);
   }
   f_str->res = res;
   f_str->max_res = maxn0;

   *f_arg = f_str;
}


/* pstring1 is a message to the manager, currently 512 */
/* pstring2 is the same information, but in a markx==10 format */
void
get_param (const struct pstruct *ppstr,
	   char **pstring1, char *pstring2,
	   struct score_count_s *s_info)
{
#ifdef FASTS
#ifndef TFAST
  char *pg_str="FASTS";
#else
  char *pg_str="TFASTS";
#endif
#endif

#ifdef FASTM
#ifndef TFAST
  char *pg_str="FASTM";
#else
  char *pg_str="TFASTM";
#endif
#endif

  sprintf (pstring1[0], "%s (%s)",pg_str,verstr);
  sprintf (pstring1[1], "%s matrix (%d:%d), ktup=%d",
	   ppstr->pam_name, ppstr->pam_h,ppstr->pam_l, ppstr->param_u.fa.ktup);

  if (ppstr->param_u.fa.iniflag) strcat(pstring1[0]," init1");

  if (pstring2 != NULL) {
    sprintf (pstring2, "; pg_name_alg: %s\n; pg_ver_rel: %s\n; pg_matrix: %s (%d:%d)\n\
; pg_gap-pen: %d %d\n; pg_ktup: %d\n",
	     pg_str,verstr,ppstr->pam_name, ppstr->pam_h,ppstr->pam_l, ppstr->gdelval,
	     ppstr->ggapval,ppstr->param_u.fa.ktup);
   }
}

void
close_work (const unsigned char *aa0, const int n0,
	    struct pstruct *ppst,
	    struct f_struct **f_arg)
{
  struct f_struct *f_str;
  int i, j;

  f_str = *f_arg;

  if (f_str != NULL) {

    free(f_str->res);
#ifdef TFAST
    free(f_str->aa1x - 1); /* because f_str->aa1x got ++'ed when allocated! */
#endif
    free(f_str->diag);
    free(f_str->l_end);
    free(f_str->link);
    free(f_str->pamh2); 
    free(f_str->pamh1);
    free(f_str->harr);
    free(f_str->vmax);
    free(f_str->vptr);
    free(f_str->sarr);
    free(f_str->aa0l);
    free(f_str->aa0s);
    free(f_str->aa0i);
    free(f_str->aa0e);
    free(f_str->aa0b);
    free(f_str->aa0ti);
    free(f_str->aa0t);
    free(f_str->nmoff);
    free(f_str->nm_u);

    if(f_str->dotat) {
      for(i = 0, j = (1 << f_str->nm0) - 1 ; i < j ; i++) {
	free(f_str->tatprobs[i]->probs);
	free(f_str->tatprobs[i]);
	free(f_str->intprobs[i]);
      }
      free(f_str->tatprobs);
      free(f_str->intprobs);
    }

    free(f_str->priors);
    free(f_str);
    *f_arg = NULL;
  }
}

void do_fasts (const unsigned char *aa0, const int n0,
	       const unsigned char *aa1, const int n1,
	       const struct pstruct *ppst, struct f_struct *f_str,
	       struct rstruct *rst, int *hoff, int opt_prob,
	       int maxsav)
{
   int     nd;		/* diagonal array size */
   int     lhval;
   int     kfact;
   register struct dstruct *dptr;
   register int tscor;
   register struct dstruct *diagp;
   struct dstruct *dpmax;
   register int lpos;
   int     tpos;
   struct savestr *vmptr, *vmaxmax;
   int     scor, tmp;
   int     im, ib, nsave, i;
   int     cmps ();		/* comparison routine for ksort */
   int ktup;
   int doffset;


   vmaxmax = &f_str->vmax[maxsav];

   ktup = ppst->param_u.fa.ktup;

   if (n1 < ktup) {
     rst->score[0] = rst->score[1] = rst->score[2] = 0;
     rst->escore = 1.0;
     rst->segnum = 0;
     rst->seglen = 0;
     return;
   }

   if (n0+n1+1 >= MAXDIAG) {
     fprintf(stderr,"n0,n1 too large: %d, %d\n",n0,n1);
     rst->score[0] = rst->score[1] = rst->score[2] = -1;
     rst->escore = 2.0;
     rst->segnum = 0;
     rst->seglen = 0;
     return;
   }

   nd = n0 + n1;

   dpmax = &f_str->diag[nd];
   for (dptr = &f_str->diag[f_str->ndo]; dptr < dpmax;)
   {
      dptr->stop = -1;
      dptr->dmax = NULL;
      dptr++->score = 0;
   }

   for (vmptr = f_str->vmax; vmptr < vmaxmax; vmptr++) {
      vmptr->score = 0;
      vmptr->exact = 0;
   }
   f_str->lowmax = f_str->vmax;
   f_str->lowscor = 0;

   /* start hashing */
   diagp = &f_str->diag[f_str->noff];
   for (lhval=lpos=0; lpos < n1; lpos++, diagp++) {
     if (ppst->hsq[aa1[lpos]]>=NMAP) {	/* skip residue */
       lpos++ ; diagp++;
       while (lpos < n1 && ppst->hsq[aa1[lpos]]>=NMAP) {lpos++; diagp++;}
       if (lpos >= n1) break;
       lhval = 0;
     }

     lhval = ((lhval & f_str->hmask) << f_str->kshft) + ppst->hsq[aa1[lpos]];

     for (tpos = f_str->harr[lhval]; tpos >= 0; tpos = f_str->link[tpos]) {

       dptr = &diagp[-tpos];

       if (f_str->l_end[tpos]) {
	 if (dptr->score + f_str->pamh1[aa0[tpos]] == f_str->aa0s[tpos]) {
	   dptr->stop = lpos;
	   dptr->score = f_str->aa0s[tpos];
	   savemax(dptr, f_str, maxsav, 1, tpos);
	   dptr->dmax = NULL;
	 }

	 else if (dptr->score + f_str->pamh1[aa0[tpos]] > f_str->aa0s[tpos]) {
	   /*
	   fprintf(stderr,"exact match score too high: %d:%d %d < %d + %d - %d:%d - %d > %d\n",
		   tpos, lpos, f_str->aa0s[tpos],dptr->score, f_str->pamh1[aa0[tpos]],
		   dptr->start, dptr->stop,
		   dptr->stop - dptr->start, f_str->aa0l[tpos]);
	   */
	   dptr->stop = lpos;
	   dptr->start = lpos - f_str->aa0l[tpos];
	   dptr->score = f_str->aa0s[tpos];
	   savemax(dptr, f_str, maxsav, 1, tpos);
	   dptr->dmax = NULL;
	 }
       }
       else if ((tscor = dptr->stop) >= 0) {
	 tscor++;	/* tscor is stop of current, increment it */
	 if ((tscor -= lpos) <= 0) {  /* tscor, the end of the current
					 match, is before lpos, so there
					 is a mismatch - this is also the
					 mismatch cost */
	   tscor *= 2;
	   scor = dptr->score;	/* save the run score on the diag */
	   if ((tscor += (kfact = f_str->pamh2[lhval])) < 0 
	       && f_str->lowscor < scor) {
	     /* if what we will get (tscor + kfact) is < 0 and the
		score is better than the worst savemax() score, save
		it */
	     savemax (dptr, f_str, maxsav,0,-1);
	   }

	   /* if extending is better than starting over, extend */
	   if ((tscor += scor) >= kfact) {
	     dptr->score = tscor;
	     dptr->stop = lpos;
	     if (f_str->l_end[tpos]) {
	       if (dptr->score == f_str->aa0s[tpos]) {
		 savemax(dptr, f_str, maxsav,1,tpos);
		 dptr->dmax = NULL;
	       }
	       else if (dptr->score > f_str->lowscor)
		 savemax(dptr, f_str, maxsav,0,tpos);
	     }
	   }
	   else {     /* otherwise, start new */
	     dptr->score = kfact;
	     dptr->start = dptr->stop = lpos;
	   }
	 } 
	 else { /* tscor is after lpos, so extend one residue */
	   dptr->score += f_str->pamh1[aa0[tpos]];
	   dptr->stop = lpos;
	   if (f_str->l_end[tpos]) {
	     if (dptr->score == f_str->aa0s[tpos]) {
	       savemax(dptr, f_str, maxsav,1,tpos);
	       dptr->dmax = NULL;
	     }
	     else if (dptr->score > f_str->lowscor)
	       savemax(dptr, f_str, maxsav,0,tpos);
	   }
	 }
       }
       else {	/* start new */
	 dptr->score = f_str->pamh2[lhval];
	 dptr->start = dptr->stop = lpos;
       }
     }				/* end tpos */
   }				/* end lpos */

   for (dptr = f_str->diag; dptr < dpmax;) {
     if (dptr->score > f_str->lowscor) savemax (dptr, f_str, maxsav,0,-1);
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

   for (nsave=0, vmptr=f_str->vmax; vmptr< vmaxmax; vmptr++) {
      if (vmptr->score > 0) {
	/*

	fprintf(stderr,"%c 0: %4d-%4d  1: %4d-%4d  dp: %d score: %d",
		(vmptr->exact ? 'x' : ' '),
		f_str->noff+vmptr->start-vmptr->dp,
		f_str->noff+vmptr->stop-vmptr->dp,
		vmptr->start,vmptr->stop,
		vmptr->dp,vmptr->score);
	*/
	vmptr->score = spam (aa0, aa1, n1, vmptr, ppst->pam2[0], f_str);
	/*
	fprintf(stderr,"  sscore: %d %d-%d\n",vmptr->score,vmptr->start,vmptr->stop);
	*/
	if (vmptr->score > 0) f_str->vptr[nsave++] = vmptr;
      }
   }

   if (nsave <= 0) {
     rst->score[0] = rst->score[1] = rst->score[2] = 0;
     rst->escore = 1.0;
     rst->segnum = 0;
     rst->seglen = 0;
     f_str->nsave = 0;
     return;
   }

   /*
   fprintf(stderr,"n0: %d; n1: %d; noff: %d\n",n0,n1,f_str->noff);
   for (ib=0; ib<nsave; ib++) {
     fprintf(stderr,"%c 0: %4d-%4d  1: %4d-%4d  dp: %d score: %d\n",
	     f_str->vptr[ib]->exact ? 'x' : ' ',
	     f_str->noff+f_str->vptr[ib]->start-f_str->vptr[ib]->dp,
	     f_str->noff+f_str->vptr[ib]->stop-f_str->vptr[ib]->dp,
	     f_str->vptr[ib]->start,f_str->vptr[ib]->stop,
	     f_str->vptr[ib]->dp,f_str->vptr[ib]->score);
   }

   fprintf(stderr,"---\n");
   */
   kssort(f_str->vptr,nsave);

   /* make certain each seg is used only once */

   for (ib=0; ib<f_str->nm0; ib++) f_str->nm_u[ib]=0;
   for (ib=0; ib < nsave; ib++) {
     doffset = f_str->vptr[ib]->dp - f_str->noff;
     tpos=f_str->aa0i[f_str->vptr[ib]->start - doffset];
     if (f_str->nm_u[tpos] == 0) {
       f_str->nm_u[tpos]=1;
     } else {
       f_str->vptr[ib]->score = -1;
     }
   }

   kssort(f_str->vptr,nsave);
   for (ib = nsave-1; ib >= 0; ib--) {
     if (f_str->vptr[ib]->score > -1) break;
   }
   nsave = ib+1;

#ifdef DEBUG
   /*
   for (ib = 0; ib < nsave; ib++) {
     if (f_str->vptr[ib]->score > 1000) {
       fprintf(stderr," score[%d] too high: %d\n",ib, f_str->vptr[ib]->score);
       for (i=0; i< 10; i++) {
	 fprintf(stderr, "%c:%d ",ppst->sq[aa1[i]],aa1[i]);
       }
       fprintf(stderr,"\n");

       f_str->vptr[ib]->score = 0;
     }
   }
   */
#endif

   scor = sconn (f_str->vptr, nsave, 
		 f_str, rst, ppst, aa0, n0, aa1, n1,
		 opt_prob);

   if (rst->escore < 0.0) rst->escore = 2.0;
   kssort(f_str->vptr,nsave);

   /* here we should use an nsave that is consistent with sconn and nm0 */

   f_str->nsave = nsave;
   if (nsave > f_str->nm0) f_str->nsave = f_str->nm0;

   rst->score[1] = f_str->vptr[0]->score;
   rst->score[0] = rst->score[2] = max(scor, f_str->vptr[0]->score);

}

void do_work (const unsigned char *aa0, const int n0,
	      const unsigned char *aa1, const int n1,
	      int frame,
	      const struct pstruct *ppst, struct f_struct *f_str,
	      int qr_flg, int shuff_flg, struct rstruct *rst,
	      struct score_count_s *s_info)
{
  int opt_prob;
  int hoff, n10, i;

  if (qr_flg==1 && f_str->shuff_cnt <= 0) {
    rst->valid_stat = 0;
    rst->escore = 2.0;
    rst->score[0]=rst->score[1]=rst->score[2]= -1;
    return;
  }
  rst->valid_stat = 1;

  s_info->s_cnt[ppst->score_ix]++;
  s_info->tot_scores++;

  if (f_str->dotat || ppst->zsflag == 4 || ppst->zsflag == 14 ) opt_prob=1;
  else opt_prob = 0;
  if (ppst->zsflag == 2 || ppst->zsflag == 12) opt_prob = 0;
  if (qr_flg) {
    opt_prob=1;
    /*    if (frame==1) */
      f_str->shuff_cnt--;
  }

  if (n1 < ppst->param_u.fa.ktup) {
    rst->score[0] = rst->score[1] = rst->score[2] = -1;
    rst->escore = 2.0;
    return;
  }
#ifdef TFAST
  n10=aatran(aa1,f_str->aa1x,n1,frame);
  if (ppst->debug_lib)
    for (i=0; i<n10; i++)
      if (f_str->aa1x[i]>ppst->nsq) {
	fprintf(stderr,
		"residue[%d/%d] %d range (%d)\n",i,n1,
		f_str->aa1x[i],ppst->nsq);
	f_str->aa1x[i]=0;
	n10=i-1;
      }

  do_fasts (aa0, n0, f_str->aa1x, n10, ppst, f_str, rst, &hoff, opt_prob, f_str->maxsav);
#else	/* FASTA */
  do_fasts (aa0, n0, aa1, n1, ppst, f_str, rst, &hoff, opt_prob, f_str->maxsav);
#endif

  rst->comp = rst->H = -1.0;
}

void do_opt (const unsigned char *aa0, const int n0,
	     const unsigned char *aa1, const int n1,
	     int frame,
	     struct pstruct *ppst, struct f_struct *f_str,
	     struct rstruct *rst)
{
  int lag, tscore, hoff, n10;

#ifdef TFAST
  n10=aatran(aa1,f_str->aa1x,n1,frame);
  do_fasts (aa0, n0, f_str->aa1x, n10, ppst, f_str, rst, &hoff, 1, f_str->maxsav);
#else	/* FASTA */
  do_fasts(aa0,n0,aa1,n1,ppst,f_str,rst, &hoff, 1, f_str->maxsav);
#endif
}


/* modify savemax() so that full length 100% matches are marked
   so that they cannot be removed - if we have a 100% match, mark "exact"

   modify savemax() to split alignments that include a comma
*/

/* savemax(dptr, f_str, maxsav) takes a current diagonal run (saved in dptr),
   and places it in the set of runs to be saved (in  f_str->vmax[])
*/

void 
savemax (struct dstruct *dptr, struct f_struct *f_str, int maxsav,
	 int exact, int tpos)
{
  register int dpos;	/* position along the diagonal, -n0 .. n1 */
  int i, j, lowj;
  register struct savestr *vmptr;
  struct savestr *vmaxmax;

  vmaxmax = &f_str->vmax[maxsav];

  dpos = (int) (dptr - f_str->diag);	/* current diagonal */

/* check to see if this is the continuation of a run that is already saved */
/* if we are at the end of the query, save it regardless */

/*  if (t_end > 0 && t_end < dptr->stop - dptr->start) {return;} */

  if ((vmptr = dptr->dmax) != NULL	/* have an active run */
      && vmptr->dp == dpos &&		/* on the correct diagonal */
      vmptr->start == dptr->start) {	/* and it starts at the same place */
    vmptr->stop = dptr->stop;	/* update the end of the match in vmax[] */

    if (exact == 1) {
    /*
      fprintf(stderr,"have cont exact match: %d - %d:%d %d:%d = %d\n",
	      dptr->score, dptr->start, dptr->stop,
	      vmptr->start, vmptr->stop, dptr->stop - dptr->start+1);
    */
      exact = 1;
    }


/* if the score is worse, don't update, return - if the score gets bad
   enough, it will restart in the diagonal scan */
    if ((i = dptr->score) <= vmptr->score) { return;} 

/* score is better, update */
    vmptr->score = i;

    vmptr->exact = exact;
/* if the score is not the worst, return */
    if (vmptr != f_str->lowmax) { return;}
  }
  else {	/* not a continuation */
    /* save in the lowest place */
    /*
    fprintf(stderr," Replacing: %d - %d:%d => %d - %d:%d",
	    f_str->lowmax->score, f_str->lowmax->start, f_str->lowmax->stop,
	    dptr->score, dptr->start, dptr->stop);
    */

    vmptr = f_str->lowmax;

    /*
    if (exact == 1) {
      fprintf(stderr,"have new exact match: %d - %d:%d = %d\n",
	      dptr->score, dptr->start, dptr->stop, dptr->stop - dptr->start+1);
    }
    */
    vmptr->exact = exact;

    i = vmptr->score = dptr->score;   /* 'i' is used as a bound */
    vmptr->dp = dpos;
    vmptr->start = dptr->start;
    vmptr->stop = dptr->stop;
    dptr->dmax = vmptr;
  }

  /* rescan the list for the worst score */
  for (vmptr = f_str->vmax;  vmptr < &f_str->vmax[maxsav] ; vmptr++) {
    if (vmptr->score < i && !vmptr->exact) {
      i = vmptr->score;
      f_str->lowmax = vmptr;
    }
  }

  f_str->lowscor = i;
}

/* this version of spam scans the diagonal to find the best local score,
   then resets the boundaries for a global alignment and re-scans */

/* NOOVERHANG allows one to score any overhanging alignment as zero.
   Useful for SAGE alignments.  Normally, one allows overhangs because
   of the possibility of partial sequences.
*/

#undef NOOVERHANG

/* 
   May, 2005 - spam() has an intesting bug that occurs when two
   peptides match in order, separated by one position (the comma).  In
   this case, spam() splits the match, and only returns the better of
   the two matches.  So, if spam splits an alignment at a comma, it
   needs the ability to insert the missing match.

*/

int spam (const unsigned char *aa0, const unsigned char *aa1,int n1,
	  struct savestr *dmax, int **pam2,
	  struct f_struct *f_str)
{
   int     lpos, doffset;
   int     tot, mtot;
   struct {
     int  start, stop, score;
   } curv, maxv;
   register const unsigned char *aa0p, *aa1p;

   curv.start = dmax->start;
   aa1p = &aa1[dmax->start];
   doffset = dmax->dp - f_str->noff;
   aa0p = &aa0[dmax->start - doffset];

   tot = curv.score = maxv.score = 0;
   for (lpos = dmax->start; lpos <= dmax->stop; lpos++) {
     tot += pam2[*aa0p++][*aa1p++];
     if (tot > curv.score) {
       curv.stop = lpos;	/* here, curv.stop is actually curv.max */
       curv.score = tot;
      }
      else if (tot < 0) {
	if (curv.score > maxv.score) {
	  maxv.start = curv.start;
	  maxv.stop = curv.stop;
	  maxv.score = curv.score;
	}
	tot = curv.score = 0;
	curv.start = lpos+1;
      }
   }

   if (curv.score > maxv.score) {
     maxv.start = curv.start;
     maxv.stop = curv.stop;
     maxv.score = curv.score;
   }

   if (maxv.score <= 0) return 0;

   /* now, reset the boundaries of the alignment using aa0b[]
      and aa0e[], which specify the residues that start and end
      the segment */
      
   maxv.start = f_str->aa0b[maxv.stop-doffset] + doffset;
   if (maxv.start < 0) {
     maxv.start = 0;
#ifdef NOOVERHANG
     return 0;
#endif
   }

   maxv.stop = f_str->aa0e[maxv.stop-doffset] + doffset;
   if (maxv.stop > n1) {
     maxv.stop = n1-1;
#ifdef NOOVERHANG
     return 0;
#endif
   }
   aa1p = &aa1[lpos = maxv.start];
   aa0p = &aa0[lpos - doffset];

   for (tot=0; lpos <= maxv.stop; lpos++) {
     tot += pam2[*aa0p++][*aa1p++];
   }

   maxv.score = tot;

/*	if (maxv.start != dmax->start || maxv.stop != dmax->stop)
		printf(" new region: %3d %3d %3d %3d\n",maxv.start,
			dmax->start,maxv.stop,dmax->stop);
*/
   dmax->start = maxv.start;
   dmax->stop = maxv.stop;

   return maxv.score;
}

int sconn (struct savestr **v, int n, 
	   struct f_struct *f_str,
	   struct rstruct *rst, const struct pstruct *ppst,
	   const unsigned char *aa0, int n0,
	   const unsigned char *aa1, int n1, int opt_prob)
{
   int     i, si, cmpp ();
   struct slink *start, *sl, *sj, *so, *sarr;
   int     lstart, ltmp, tstart, plstop, ptstop, ptstart, tstop;
   double  tatprob;
   int     dotat;

   sarr = f_str->sarr;

   /*  sort the score left to right in lib pos */
   kpsort (v, n);

   start = NULL;
   rst->score[0] = 0;
   rst->escore = 2.0;

/*  for the remaining runs, see if they fit */
/*  lstart/lstop -> start/stop in library sequence
    tstart/tstop -> start/stop in query sequence
    plstart/plstop ->
*/

   for (i = 0, si = 0; i < n; i++) {

     /* the segment is worth adding; find out where? */
     lstart = v[i]->start;
     ltmp = v[i]->stop;
     tstart = lstart - v[i]->dp + f_str->noff;
     tstop = ltmp - v[i]->dp + f_str->noff;

     /*	put the run in the group */
     sarr[si].vp = v[i];
     sarr[si].score = v[i]->score;
     sarr[si].next = NULL;
     sarr[si].prev = NULL;
     sarr[si].tat = NULL;

/*
  opt_prob for FASTS only has to do with using aa1 for priors,
  i.e. we always calculate tatprobs for segments in FASTS (unlike
  FASTF)
*/
     if(opt_prob) {
       sarr[si].tatprob = 
	 calc_tatusov(NULL, &sarr[si], aa0, n0, aa1, n1, 
		      ppst->pam2[0], ppst->nsq, f_str, 
		      ppst->pseudocts, opt_prob, ppst->zsflag);
       if (sarr[si].tatprob < 0.0) {
	 fprintf(stderr," negative tatprob: %lg\n",sarr[si].tatprob);
	 sarr[si].tatprob = 1.0;
       }
       sarr[si].tat = sarr[si].newtat;
     }

/*  if it fits, then increase the score

    start points to the highest scoring run
    -> next is the second highest, etc.
    put the segment into the highest scoring run that it fits into
*/
     for (sl = start; sl != NULL; sl = sl->next) {
       ltmp = sl->vp->start;
 /* plstop -> previous lstop */
       plstop = sl->vp->stop;
 /* ptstart -> previous t(query) start */
       ptstart = ltmp - sl->vp->dp + f_str->noff;
 /* ptstop -> previous t(query) stop */
       ptstop = plstop - sl->vp->dp + f_str->noff;
#ifndef FASTM
 /* if the previous library stop is before the current library start */
       if (plstop < lstart && ( ptstop < tstart || ptstart > tstop))
#else
 /* if the previous library stop is before the current library start */
       if (plstop < lstart && ptstop < tstart)
#endif
       {
	 if(!opt_prob) {
	    sarr[si].score = sl->score + v[i]->score;
	    sarr[si].prev = sl;
	    break;
	  } else {
	    tatprob = calc_tatusov(sl, &sarr[si], aa0, n0, aa1, n1, 
				   ppst->pam2[0], ppst->nsq, f_str, 
				   ppst->pseudocts, opt_prob, ppst->zsflag);
	    /* if our tatprob gets worse when we add this, forget it */
	    if(tatprob > sarr[si].tatprob) {
	      free(sarr[si].newtat->probs); /* get rid of new tat struct */
	      free(sarr[si].newtat);
	      continue; /* reuse this sarr[si] */
	    } else {
	      sarr[si].tatprob = tatprob;
	      free(sarr[si].tat->probs); /* get rid of old tat struct */
	      free(sarr[si].tat);
	      sarr[si].tat = sarr[si].newtat;
	      sarr[si].prev = sl;
	      sarr[si].score = sl->score + v[i]->score;
	      /*
		fprintf(stderr,"sconn %d added %d:%d getting %d; si: %d, tat: %g\n",
		i,v[i]->start, v[i]->score,sarr[si].score,si, tatprob);
	      */
	      break;
	    }
	  }
	}
      }
      
      /* now recalculate where the score fits */
      if (start == NULL) start = &sarr[si];
      else {
	if(!opt_prob) {
	  for (sj = start, so = NULL; sj != NULL; sj = sj->next) {
	    if (sarr[si].score > sj->score) {
	      sarr[si].next = sj;
	      if (so != NULL)
		so->next = &sarr[si];
	      else
		start = &sarr[si];
	      break;
	    }
	    so = sj;
	  }
	} else {
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
     if(opt_prob) {
       rst->escore = start->tatprob;
     } else {
       rst->escore = 2.0;
     }

     rst->segnum = rst->seglen = 0;
     for(sj = start ; sj != NULL; sj = sj->prev) {
       rst->segnum++;
       rst->seglen += sj->vp->stop - sj->vp->start + 1;
     }
     return (start->score);
   } else {
     rst->escore = 1.0;
   }

   rst->segnum = rst->seglen = 0;
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

/* calculate the 100% identical score */
int
shscore(const unsigned char *aa0, const int n0, int **pam2, int nsq)
{
  int i, sum;
  for (i=0,sum=0; i<n0; i++)
    if (aa0[i] != EOSEQ && aa0[i]<=nsq) sum += pam2[aa0[i]][aa0[i]];
  return sum;
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
  struct rstruct rst;
  int ib, i;
  unsigned char *aa0t;
  const unsigned char *aa1p;
  struct savestr *vmptr;

  /* the a_res for this function must always be re-calculated */

  *have_ares = 0x02;	/* set 0x2 bit to indicate local copy */

  if ((a_res = (struct a_res_str *)calloc(1, sizeof(struct a_res_str)))==NULL) {
    fprintf(stderr," [do_walign] Cannot allocate a_res");
    return NULL;
  }

#ifdef TFAST
  f_str->n10 = n10 = aatran(aa1,f_str->aa1x,n1,frame);
  aa1p = f_str->aa1x;
#else
  n10 = n1;
  aa1p = aa1;
#endif

  do_fasts(aa0, n0, aa1p, n10, ppst, f_str, &rst, &hoff, 1, f_str->maxsav_w);
  a_res->sw_score = rst.score[0];
  memcpy(&a_res->rst, &rst, sizeof(rst));

  /* the alignment portion takes advantage of the information left
     over in f_str after do_fasts is done.  in particular, it is
     easy to run a modified sconn() to produce the alignments.

     unfortunately, the alignment display routine wants to have
     things encoded as with bd_align and sw_align, so we need to do that.
  */

  /* unnecessary; do_fasts just did this */
  /*  kssort(f_str->vptr,f_str->nsave);  */

  /* at some point, we want one best score for each of the segments */


  a_res->min0 = a_res->min1 = a_res->max0 = a_res->max1 = 0;
  if (f_str->nsave <=0) {
    a_res->res = NULL;
    return a_res;
  }

  for ( ; f_str->nsave > 0; f_str->nsave--) 
    if (f_str->vptr[f_str->nsave-1]->score >0) break;

  if ((aa0t = (unsigned char *)calloc(n0+1,sizeof(unsigned char)))==NULL) {
    fprintf(stderr," cannot allocate aa0t %d\n",n0+1);
    exit(1);
  }

  /* copy aa0[] into f_str->aa0t[] */
  for (i=0; i<n0; i++) f_str->aa0t[i] = aa0t[i] = aa0[i];
  f_str->aa0t[i] = aa0t[i] = '\0';

  a_res->nres = sconn_a (aa0t,n0,aa1p,n10,f_str, a_res, ppst);

  free(aa0t);

  a_res->res = f_str->res;
  return a_res;
}

/* this version of sconn is modified to provide alignment information */
/* in addition, it needs to know whether a segment has been used before */

/* sconn_a fills in the res[nres] array, but this is passed implicitly
   through f_str->res[f_str->nres] */

int sconn_a (unsigned char *aa0, int n0,
	     const unsigned char *aa1, int n1,
	     struct f_struct *f_str,
	     struct a_res_str *a_res,
	     struct pstruct *ppst)
{
   int     i, si, cmpp (), n;
   unsigned char *aa0p;
   int sx, dx, doff, *aa0tip;

   struct savestr **v;
   struct slink *start, *sl, *sj, *so, *sarr;
   int     lstart, lstop, ltmp, plstart, tstart, plstop, ptstop, ptstart, tstop;

   int *res, nres, tres;

   double tatprob;

/*	sort the score left to right in lib pos */

   v = f_str->vptr;
   n = f_str->nsave;
   sarr = f_str->sarr;

   /* set things up in case nothing fits */
   if (n <=0 || v[0]->score <= 0) return 0;

   if (v[0]->score < 0) {
     sarr[0].vp = v[0];
     sarr[0].score = v[0]->score;
     sarr[0].next = NULL;
     sarr[0].prev = NULL;
     start = &sarr[0];
   }
   else {

     krsort (v, n);	/* sort from left to right in library */

     start = NULL;

     /*	for each alignment, see if it fits */


     for (i = 0, si = 0; i < n; i++) {
       /*	if the score is less than the join threshold, skip it */

       if (v[i]->score < 0) continue;

       lstart = v[i]->start;
       lstop = v[i]->stop;
       tstart = lstart - v[i]->dp + f_str->noff;
       tstop = lstop - v[i]->dp + f_str->noff;

       /*	put the alignment in the group */

       sarr[si].vp = v[i];
       sarr[si].score = v[i]->score;
       sarr[si].next = NULL;
       sarr[si].prev = NULL;
       sarr[si].tat = NULL;

       sarr[si].tatprob = 
	 calc_tatusov(NULL, &sarr[si], aa0, n0, aa1, n1, 
		      ppst->pam2[0], ppst->nsq, f_str, 
		      ppst->pseudocts, 1, ppst->zsflag);
       sarr[si].tat = sarr[si].newtat;


       /* 	if it fits, then increase the score */
       /* start points to a sorted (by total score) list of candidate
	  overlaps */

       for (sl = start; sl != NULL; sl = sl->next) { 
	 plstart = sl->vp->start;
	 plstop = sl->vp->stop;
	 ptstart = plstart - sl->vp->dp + f_str->noff;
	 ptstop = plstop - sl->vp->dp + f_str->noff;
#ifndef FASTM
	 if (plstart > lstop && (ptstop < tstart || ptstart > tstop)) {
#else
         if (plstop > lstart && ptstart > tstop) {
#endif
	   /* alignment always uses probabilistic scoring ... */
	   /*   sarr[si].score = sl->score + v[i]->score;
		sarr[si].prev = sl;
		break; */		/* quit as soon as the alignment has been added */

	   tatprob = calc_tatusov(sl, &sarr[si], aa0, n0, aa1, n1, 
				  ppst->pam2[0], ppst->nsq, f_str, 
				  ppst->pseudocts, 1, ppst->zsflag);
	   /* if our tatprob gets worse when we add this, forget it */
	   if(tatprob > sarr[si].tatprob) {
	     free(sarr[si].newtat->probs); /* get rid of new tat struct */
	     free(sarr[si].newtat);
	     continue; /* reuse this sarr[si] */
	   } else {
	     sarr[si].tatprob = tatprob;
	     free(sarr[si].tat->probs); /* get rid of old tat struct */
	     free(sarr[si].tat);
	     sarr[si].tat = sarr[si].newtat;
	     sarr[si].prev = sl;
	     sarr[si].score = sl->score + v[i]->score;
	     /*
	       fprintf(stderr,"sconn %d added %d/%d getting %d; si: %d, tat: %g\n",
	       i,v[i]->start, v[i]->score,sarr[si].score,si, tatprob);
	     */
	     break;
	   }
	 }
       }

       /* now recalculate the list of best scores */
       if (start == NULL)
	 start = &sarr[si];	/* put the first one in the list */
       else
	 for (sj = start, so = NULL; sj != NULL; sj = sj->next) {
	   /* if (sarr[si].score > sj->score) { */ /* new score better than old */
	   if ( sarr[si].tatprob < sj->tatprob ||
		((sarr[si].tatprob == sj->tatprob) && sarr[si].score > sj->score) ) {
	     sarr[si].next = sj;		/* next best after new score */
	     if (so != NULL)
	       so->next = &sarr[si];	/* prev_best->next points to best */
	     else  start = &sarr[si];	/* start points to best */
	     break;			/* stop looking */
	   }
	   so = sj;		/* previous candidate best */
	 }
       si++;				/* increment to next alignment */
     }
   }

   for (i = 0 ; i < si ; i++) {
     free(sarr[i].tat->probs);
     free(sarr[i].tat);
   }

   res = f_str->res;
   tres = nres = 0;
   aa0p = aa0;
   aa0tip = f_str->aa0ti;	/* point to temporary index */
   a_res->min1 = start->vp->start;
   a_res->min0 = 0;

   sx=start->vp->start-start->vp->dp+f_str->noff;
   f_str->aa0t_off = sx - f_str->aa0b[sx];

   for (sj = start; sj != NULL; sj = sj->prev ) {
     doff = (int)(aa0p-aa0) - (sj->vp->start-sj->vp->dp+f_str->noff);
     
     /* fprintf(stderr,"doff: %3d\n",doff); */
     
     for (dx=sj->vp->start,sx=sj->vp->start-sj->vp->dp+f_str->noff;
	  dx <= sj->vp->stop; dx++) {
       *aa0tip++ = f_str->aa0i[sx];	/* save index */
       *aa0p++ = f_str->aa0t[sx++];	/* save sequence at index */
       tres++;
       res[nres++] = 0;
     }
     sj->vp->dp -= doff;
     if (sj->prev != NULL) {
       if (sj->prev->vp->start - sj->vp->stop - 1 > 0 )
	 tres += res[nres++] = (sj->prev->vp->start - sj->vp->stop - 1);
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
     a_res->max1 = sj->vp->stop+1;
     a_res->max0 = a_res->max1 - sj->vp->dp + f_str->noff;
   }
   *aa0p = '\0';	/* be sure to terminate the string */

   /*
   fprintf(stderr,"(%3d - %3d):(%3d - %3d)\n",
	   a_res->min0,a_res->max0,a_res->min1,a_res->max1);
   */
   
   /* now replace f_str->aa0t with aa0
      (f_str->aa0t is permanent, aa0 is not)*/
   for (i=0; i<n0; i++) f_str->aa0t[i] = aa0[i];

   return tres;
}

/* for fasts (and fastf), pre_cons needs to set up f_str as well as do
   necessary translations - for right now, simply do do_walign */

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
  aln->qlfact= 1;
  aln->llfact = aln->llmult = 3;
  if (frame > 3) aln->llrev = 1;
  else aln->llrev = 0;
  aln->frame = 0;
#else	/* FASTS */
  aln->llfact = aln->llmult = aln->qlfact = 1;
  aln->llrev = aln->qlrev = 0;
  aln->frame = 0;
#endif
}

void aaptrshuffle(unsigned char *res, int n) {

  int i, j;
  unsigned char tmp;

  for( i = n; --i; ) {

    /* j = nrand(i); if (i == j) continue; */ /* shuffle */
    j = (n - 1) - i; if (i <= j ) break; /* reverse */

    tmp = res[i];
    res[i] = res[j];
    res[j] = tmp;
  }
}

void aa0shuffle(unsigned char *aa0, int n0, struct f_struct *f_str) {

  int i;
  int j;

  for(i = 0 ; i < f_str->nm0 ; i++) { /* for each fragment */

    aaptrshuffle(&(aa0[f_str->nmoff[i]]), 
		 f_str->nmoff[i+1] - f_str->nmoff[i] - 1 );

  }

}
