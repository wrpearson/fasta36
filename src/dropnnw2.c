/* $Id: dropnnw2.c $ */
/* $Revision: 1140 $  */

/* copyright (c) 1996, 2007, 2014 by William R. Pearson and The Rector &
   Visitors of the Univeristy of Virginia */

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

/* 4-April-2007 - convert to global alignment */

/* 17-Aug-2006 - removed globals *sapp/last - alignment should be thread safe */

/* 12-Oct-2005 - converted to use a_res and aln for alignment coordinates */

/* 4-Nov-2004 - Diagonal Altivec Smith-Waterman included */

/* 14-May-2003 - modified to return alignment start at 0, rather than
   1, for begin:end alignments

   25-Feb-2003 - modified to support Altivec parallel Smith-Waterman

   22-Sep-2003 - removed Altivec support at request of Sencel lawyers
*/

/* this code uses an implementation of the Smith-Waterman algorithm
   designed by Phil Green, U. of Washington, that is 1.5 - 2X faster
   than my Miller and Myers implementation. */

/* the shortcuts used in this program prevent it from calculating scores
   that are less than the gap penalty for the first residue in a gap. As
   a result this code cannot be used with very large gap penalties, or
   with very short sequences, and probably should not be used with prss3.
*/

/* version 3.2 fixes a subtle bug that was encountered while running
   do_walign() interspersed with do_work().  This happens only with -m
   9 and pvcomplib.  The fix was to more explicitly zero-out ss[] at
   the beginning of do_work.
*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <math.h>

#include "defs.h"
#include "param.h"

static char *verstr="6.0 April 2007";

#include "dropgsw2.h"

#define DROP_INTERN
#include "drop_func.h"

#ifdef SW_SSE2
#ifdef GLOBAL_GLOBAL
#include "global_sse2.h"
#define GLOBAL_BYTE global_sse2_byte
#define GLOBAL_WORD global_sse2_word
#else
#include "glocal_sse2.h"
#define GLOBAL_BYTE glocal_sse2_byte
#define GLOBAL_WORD glocal_sse2_word
#endif
#endif

struct swstr {int H, E;};

extern void init_karlin(const unsigned char *aa0, int n0, struct pstruct *ppst,
			double *aa0_f, double **kp);
extern int do_karlin(const unsigned char *aa1, int n1,
		     int **pam2, const struct pstruct *ppst,
		     double *aa0_f, double *kar_p, double *lambda, double *H);

extern int
NW_ALIGN(int IW, const unsigned char *B,
	 int M, int N,
	 int **W, int G, int H, int *res, int *nres
	 );

static int
FGLOBAL_ALIGN(int *pwaa, const unsigned char *aa1,
	      int n0, int n1,
	      int GG,int HH,
	      struct swstr *ss);

extern struct a_res_str *
nsw_malign (int ***pam2p, int pam_ix, int n0,
	    const unsigned char *aa1, int n1,
	    int score_thresh, int max_res,
	    int gdelval, int ggapval, 
	    struct swstr *ss, 
	    struct a_res_str *cur_ares,
	    int (*fn_walign)
	    (
	     int **pam2p, int n0,
	     const unsigned char *aa1, int n1,
	     int q, int r,
	     struct swstr *ss,
	     struct a_res_str *a_res
	     ),
	    int do_rep
	    );

static 
void DISPLAY(const unsigned char *A, const unsigned char *B, 
	     int M, int N,
	     int *S, int AP, int BP, char *sq);

extern void aancpy(char *to, char *from, int count, struct pstruct *ppst);

/* initialize for Smith-Waterman optimal score */

void
init_work (unsigned char *aa0, int n0,
	   struct pstruct *ppst,
	   struct f_struct **f_arg)
{
  int maxn0, ip;
  int *pwaa_s, *pwaa_a;
  int e, f, i, j, l;
  int *res;
  struct f_struct *f_str;
  int **pam2p;
  struct swstr *ss;
  int nsq;

#if defined(SW_ALTIVEC) || defined(SW_SSE2)
  int data, bias, ceiling, gap;
  unsigned char *  pc;
  unsigned short * ps;
  int  overflow;

  int n_count;
  int col_len;
#endif

  if (ppst->ext_sq_set) {
    nsq = ppst->nsqx; ip = 1;
  }
  else {
    /* with memory mapped databases with lc chars, always nsqx */
    nsq = ppst->nsqx; ip = 0;
  }

  /* initialize range of length appropriate */

  if (ppst->n1_low == 0 ) {
    ppst->n1_low = (int)(0.75 * (float)n0 + 0.5);
  }

#if defined(GLOBAL_GLOBAL)
  if (ppst->n1_high == BIGNUM) {
    ppst->n1_high = (int)(1.33 * (float)n0 - 0.5);
  }
#endif

  /* allocate space for function globals */
  f_str = (struct f_struct *)calloc(1,sizeof(struct f_struct));

  if((ppst->zsflag%10) == 6) {
    f_str->kar_p = NULL;
    init_karlin(aa0, n0, ppst, &f_str->aa0_f[0], &f_str->kar_p);
  }
  
  /* allocate space for the scoring arrays */
  if ((ss = (struct swstr *) calloc (n0+2, sizeof (struct swstr)))
      == NULL) {
    fprintf (stderr, "cannot allocate ss array %3d\n", n0);
    exit (1);
  }
  ss++;

  f_str->ss = ss;

  /* initialize variable (-S) pam matrix */
  if ((f_str->waa_s= (int *)calloc((nsq+1)*(n0+1),sizeof(int))) == NULL) {
    fprintf(stderr,"cannot allocate waa_s array %3d\n",nsq*n0);
    exit(1);
  }

  /* initialize pam2p[1] pointers */
  if ((f_str->pam2p[1]= (int **)calloc((n0+1),sizeof(int *))) == NULL) {
    fprintf(stderr,"cannot allocate pam2p[1] array %3d\n",n0);
    exit(1);
  }

  pam2p = f_str->pam2p[1];
  if ((pam2p[0]=(int *)calloc((nsq+1)*(n0+1),sizeof(int))) == NULL) {
    fprintf(stderr,"cannot allocate pam2p[1][] array %3d\n",nsq*n0);
    exit(1);
  }

  for (i=1; i<n0; i++) {
    pam2p[i]= pam2p[0] + (i*(nsq+1));
  }

  /* initialize universal (alignment) matrix */
  if ((f_str->waa_a= (int *)calloc((nsq+1)*(n0+1),sizeof(int))) == NULL) {
    fprintf(stderr,"cannot allocate waa_a struct %3d\n",nsq*n0);
    exit(1);
  }
   
  /* initialize pam2p[0] pointers */
  if ((f_str->pam2p[0]= (int **)calloc((n0+1),sizeof(int *))) == NULL) {
    fprintf(stderr,"cannot allocate pam2p[1] array %3d\n",n0);
    exit(1);
  }

  pam2p = f_str->pam2p[0];
  if ((pam2p[0]=(int *)calloc((nsq+1)*(n0+1),sizeof(int))) == NULL) {
    fprintf(stderr,"cannot allocate pam2p[1][] array %3d\n",nsq*n0);
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

  pwaa_s = f_str->waa_s;
  pwaa_a = f_str->waa_a;
  if (ppst->pam_pssm) {
    for (e = 0; e <nsq; e++)	{	/* for each residue in the alphabet */
      for (f = 0; f < n0; f++) {	/* for each position in aa0 */
	*pwaa_s++ = f_str->pam2p[ip][f][e] = ppst->pam2p[ip][f][e];
	*pwaa_a++ = f_str->pam2p[0][f][e]  = ppst->pam2p[0][f][e];
      }
    }
  }
  else {	/* initialize scanning matrix */
    for (e = 0; e <nsq; e++)	/* for each residue in the alphabet */
      for (f = 0; f < n0; f++)	{	/* for each position in aa0 */
	*pwaa_s++ = f_str->pam2p[ip][f][e]= ppst->pam2[ip][aa0[f]][e];
	*pwaa_a++ = f_str->pam2p[0][f][e] = ppst->pam2[0][aa0[f]][e];
      }
  }

  /* these structures are used for producing alignments */

#if defined(SW_SSE2)
  /* First we allocate memory for the workspace - i.e. two rows for H and
   * one row for F.  We also need enough space to hold a temporary
   * scoring profile which will be query_length * 16 (sse2 word length).
   * Since this might be run on Linux or AIX too, we don't assume 
   * anything about the memory allocation but align it ourselves.
   */
  f_str->workspace_memory  = (void *)malloc(3*16*(MAXTST+MAXLIB+32)+256);
  f_str->workspace  = (void *)((((size_t)f_str->workspace_memory) + 255) & (~0xff));

  /* We always use a scoring profile for the SSE2 implementation, but the layout
   * is a bit strange.  The scoring profile is parallel to the query, but is
   * accessed in a stripped pattern.  The query is divided into equal length
   * segments.  The number of segments is equal to the number of elements
   * processed in the SSE2 register.  For 8-bit calculations, the query will
   * be divided into 16 equal length parts.  If the query is not long enough
   * to fill the last segment, it will be filled with neutral weights.  The
   * first element in the SSE register will hold a value from the first segment,
   * the second element of the SSE register will hold a value from the
   * second segment and so on.  So if the query length is 288, then each
   * segment will have a length of 18.  So the first 16 bytes will  have
   * the following weights: Q1, Q19, Q37, ... Q271; the next 16 bytes will
   * have the following weights: Q2, Q20, Q38, ... Q272; and so on until
   * all parts of all segments have been written.  The last seqment will
   * have the following weights: Q18, Q36, Q54, ... Q288.  This will be
   * done for the entire alphabet.
   */

  f_str->word_score_memory = (void *)malloc((n0 + 32) * sizeof(short) * (nsq + 1) + 256);
  f_str->byte_score_memory = (void *)malloc((n0 + 32) * sizeof(char) * (nsq + 1) + 256);

  f_str->word_score = (unsigned short *)((((size_t)f_str->word_score_memory) + 255) & (~0xff));
  f_str->byte_score = (unsigned char *)((((size_t)f_str->byte_score_memory) + 255) & (~0xff));

  overflow = 0;
  gap = -2 * ppst->ggapval;

  if (ppst->pam_pssm) {
    /* Use a position-specific scoring profile. 
     * This is essentially what we are going to construct anyway, but we'll
     * reorder it to suit sse2.
     */       
    bias = 127;
    ceiling = 0;
    for (i = 1; i < nsq ; i++) {
      for (j = 0; j < n0 ; j++) {
        data = ppst->pam2p[ip][j][i];
        if (data < bias) {
          bias = data;
        }
        if (data > ceiling) {
          ceiling = data;
        }
      }
    }
    bias += gap;
    if (bias > 0) {
      bias = 0;
    }


    /* Fill our specially organized byte- and word-size scoring arrays. */
    ps = f_str->word_score;
    col_len = (n0 + 7) / 8;
    n_count = (n0 + 7) & 0xfffffff8;
    for (f = 0; f < n_count; ++f) {
      *ps++ = 0;
    }
    for (f = 1; f < nsq ; f++) {
      for (e = 0; e < col_len; e++) {
        for (i = e; i < n_count; i += col_len) {
          if ( i < n0) { data = ppst->pam2p[ip][i][f] + gap;}
          else {data = 0;}
          *ps++ = (unsigned short)(data);
        }
      }
    }
    pc = f_str->byte_score;
    col_len = (n0 + 15) / 16;
    n_count = (n0 + 15) & 0xfffffff0;
    for (f = 0; f < n_count; ++f) {
      *pc++ = 0;
    }
    for (f = 1; f < nsq ; f++) {
      for (e = 0; e < col_len; e++) {
        for (i = e; i < n_count; i += col_len) {
          if ( i < n0 ) { data = ppst->pam2p[ip][i][f] + gap;}
          else {data = 0;}
          if (data > 255) {
            printf("Fatal error. data: %d bias: %d, position: %d/%d, "
                   "Score out of range for 8-bit SSE2 datatype.\n",
                   data, bias, f, e);
            exit(1);
          }
          *pc++ = (unsigned char)(data-bias);
        }
      }
    }
  } else {
    /* Classical simple substitution matrix */
    /* Find the bias to use in the substitution matrix */
    bias = 127;
    ceiling = 0;
    for (i = 1; i < nsq ; i++) {
      for (j = 1; j < nsq ; j++) {
        data = ppst->pam2[ip][i][j];
        if (data < bias) {
          bias = data;
        }
        if (data > ceiling) {
          ceiling = data;
        }
      }
    }
    bias += gap;
    if (bias > 0) {
      bias = 0;
    }

    /* Fill our specially organized byte- and word-size scoring arrays. */
    ps = f_str->word_score;
    col_len = (n0 + 7) / 8;
    n_count = (n0 + 7) & 0xfffffff8;
    for (f = 0; f < n_count; ++f) {
      *ps++ = 0;
    }
    for (f = 1; f < nsq ; f++) {
      for (e = 0; e < col_len; e++) {
        for (i = e; i < n_count; i += col_len) {
          if (i >= n0) {
            data = 0;
          } else {
            data = ppst->pam2[ip][aa0[i]][f] + gap;
          }
          *ps++ = (unsigned short)(data);
        }
      }
    }

    pc = f_str->byte_score;
    col_len = (n0 + 15) / 16;
    n_count = (n0 + 15) & 0xfffffff0;
    for (f = 0; f < n_count; ++f) {
      *pc++ = 0;
    }
    for (f = 1; f < nsq ; f++) {
      for (e = 0; e < col_len; e++) {
        for (i = e; i < n_count; i += col_len) {
          if (i >= n0) {
            data = 0;
          } else {
            data = ppst->pam2[ip][aa0[i]][f] + gap;
          }
          if (data > 255) {
            printf("Fatal error. data: %d bias: %d, position: %d/%d, "
                   "Score out of range for 8-bit SSE2 datatype.\n",
                   data, bias, f, e);
            exit(1);
          }
          *pc++ = (unsigned char)(data-bias);
        }
      }
    }
  }
       
  f_str->ceiling = (unsigned char) (ceiling + gap - bias);
  f_str->bias = (unsigned char) (-bias);
  f_str->alphabet_size = nsq;

  /* Some variable to keep track of how many 8-bit runs we need to rerun
   * in 16-bit accuracy. If there are too many reruns it can be faster
   * to use 16-bit alignments directly. 
   */
   
  /* We can only do 8-bit alignments if the scores were small enough. */
  f_str->try_8bit = (overflow == 0) ? 1 : 0;

  f_str->done_8bit  = 0;
  f_str->done_16bit = 0;
#endif /* SW_SSE2 */

  /* minimum allocation for alignment */
  f_str->max_res = max(3*n0/2,MIN_RES);

  *f_arg = f_str;
}

void close_work (const unsigned char *aa0, int n0,
		 struct pstruct *ppst,
		 struct f_struct **f_arg)
{
  struct f_struct *f_str;

  f_str = *f_arg;

  if (f_str != NULL) {
    if (f_str->kar_p !=NULL) free(f_str->kar_p);
    f_str->ss--;
    free(f_str->ss);
    free(f_str->waa_a);
    free(f_str->pam2p[0][0]);
    free(f_str->pam2p[0]);
    free(f_str->waa_s);
    free(f_str->pam2p[1][0]);
    free(f_str->pam2p[1]);

#if defined(SW_ALTIVEC) || defined(SW_SSE2)
    free(f_str->workspace_memory);
    free(f_str->word_score_memory);
    free(f_str->byte_score_memory);
#endif
    free(f_str);
    *f_arg = NULL;
  }
}


/* pstring1 is a message to the manager, currently 512 */
/*void get_param(struct pstruct *pstr,char *pstring1)*/
void
get_param (const struct pstruct *ppst,
	   char **pstring1, char *pstring2,
	   struct score_count_s *s_info)
{

  char pg_str[120];
  char psi_str[120];

#ifdef SW_SSE2
#if defined(GLOBAL_GLOBAL) 
  char *pg_desc = "Global/Global affine Needleman-Wunsch (SSE2, Michael Farrar 2010)";
#else
  char *pg_desc = "Global/Local affine Needleman-Wunsch (SSE2, Michael Farrar 2010)";
#endif
#else
#if defined(GLOBAL_GLOBAL) 
  char *pg_desc = "Global/Global affine Needleman-Wunsch (2007)";
#else
  char *pg_desc = "Global/Local affine Needleman-Wunsch (2007)";
#endif
#endif

  strncpy(pg_str, pg_desc,  sizeof(pg_str));

  if (ppst->pam_pssm) { strncpy(psi_str,"-PSI",sizeof(psi_str));}
  else { psi_str[0]='\0';}

  sprintf (pstring1[0], "%s (%s)", pg_str, verstr);
  sprintf (pstring1[1], 
#ifdef OLD_FASTA_GAP
	   "%s matrix%s (%d:%d)%s, gap-penalty: %d/%d",
#else
	   "%s matrix%s (%d:%d)%s, open/ext: %d/%d",
#endif
	   ppst->pam_name, psi_str, ppst->pam_h,ppst->pam_l, 
	   (ppst->ext_sq_set)?"xS":"\0", ppst->gdelval, ppst->ggapval);

   if (pstring2 != NULL) {
#ifdef OLD_FASTA_GAP
     sprintf(pstring2,"; pg_name_alg: %s\n; pg_ver_rel: %s\n; pg_matrix: %s (%d:%d)%s\n; pg_gap-pen: %d %d\n",
#else
     sprintf(pstring2,"; pg_name_alg: %s\n; pg_ver_rel: %s\n; pg_matrix: %s (%d:%d)%s\n; pg_open-ext: %d %d\n",
#endif
	     pg_str,verstr,psi_str,ppst->pam_h,ppst->pam_l, 
	     (ppst->ext_sq_set)?"xS":"\0",ppst->gdelval,ppst->ggapval);
   }
}

void do_work (const unsigned char *aa0, int n0,
	      const unsigned char *aa1, int n1,
	      int frame,
	      const struct pstruct *ppst, struct f_struct *f_str,
	      int qr_flg, int shuff_flg, struct rstruct *rst,
	      struct score_count_s *s_info)
{
  int     score;
  double lambda, H;
  int i;
  
  rst->valid_stat = 1;
  s_info->s_cnt[0]++;
  s_info->tot_scores++;

#if defined(SW_SSE2)

  score = OVERFLOW_SCORE;

  if (f_str->try_8bit) {
    score = GLOBAL_BYTE(n0,
                        f_str->byte_score,
                        aa1,
                        n1,
#ifndef OLD_FASTA_GAP
                        //-(ppst->gdelval + ppst->ggapval),
                        -ppst->gdelval,
#else
                        //-ppst->gdelval,
                        -(ppst->gdelval - ppst->ggapval),
#endif
                        -ppst->ggapval,
                        f_str->ceiling,
                        f_str->bias,
                        f_str);
      
    f_str->done_8bit++;
    
    /* The 8 bit version is roughly 50% faster than the 16 bit version,
     * so we are fine if less than about 1/3 of the runs have to
     * be rerun with 16 bits. If it is more, and we have tried at least
     * 500 sequences, we switch off the 8-bit mode.
     */
    if (score == OVERFLOW_SCORE) {
      f_str->done_16bit++;
      if(f_str->done_8bit>500 && (3*f_str->done_16bit)>(f_str->done_8bit))
        f_str->try_8bit = 0;
    }
  }
      
  if (score == OVERFLOW_SCORE) {
    /* Overflow, so we have to redo it in 16 bits. */
    score = GLOBAL_WORD(n0,
                        f_str->word_score,
                        aa1,
                        n1,
#ifndef OLD_FASTA_GAP
                        //-(ppst->gdelval + ppst->ggapval),
                        -ppst->gdelval,
#else
                        //-ppst->gdelval,
                        -(ppst->gdelval - ppst->ggapval),
#endif
                        -ppst->ggapval,
                        f_str->ceiling,
                        f_str);
  }
#else

  score = FGLOBAL_ALIGN(f_str->waa_s,aa1,n0,n1,
#ifdef OLD_FASTA_GAP
                       -(ppst->gdelval - ppst->ggapval),
#else
                       -ppst->gdelval,
#endif
                       -ppst->ggapval,f_str->ss);
#endif

  rst->score[0] = score;

  if(((ppst->zsflag % 10) == 6) &&
     (do_karlin(aa1, n1, ppst->pam2[0], ppst,f_str->aa0_f, 
		f_str->kar_p, &lambda, &H)>0)) {
    rst->comp = 1.0/lambda;
    rst->H = H;
  }
  else {rst->comp = rst->H = -1.0;}

}

/* nw_walign is the equivalent of sw_walign from dropnnw.c -- it is to
   be called from nw_malign (the equivalent of sw_malign */

int
nw_walign (int **pam2p, int n0,
	   const unsigned char *aa1, int n1,
	   int q, int r,
	   struct swstr *ss,
	   struct a_res_str *a_res
	   )
{
  const unsigned char *aa1p;
  register int i, j;
  register struct swstr *ssj;
  int e, f, h, p;
  int qr, t;
  int score;
  int cost, I, J, K, L;

  qr = q + r;

  score = -BIGNUM;
  J = n0-1; L = 0;	/* alignments are global in aa0[n0] */

  /* initialize 0th row */
  ss[0].H = 0;
  ss[0].E = t = -q;	/* must be re-initialized because it was
			   filled in the reverse direction in previous
			   invocations */
  /* must count from ss+1, ss[0].H = 0 */
  for (ssj=ss+1; ssj <= ss+n0 ; ssj++) {
    ssj->H = t = t - r;
    ssj->E = t - q;
  }

  aa1p = aa1;
  i = 0;
  t = -q;
  while (*aa1p) {
    p = ss[0].H;
#ifndef GLOBAL_GLOBAL	/* GLOBAL_LOCAL */
    ss[0].H = h = t = 0;
#else	
    ss[0].H = h = t = t - r;
#endif
    f = t - q;
    /* pwaa = waa + (*aa1p++ * n0); */
    /* ssj must start at ss+1, ss[0].H = 0,
       but j must go 0 .. n0-1 for pam2p[j]  */
    for (ssj = ss+1, j=0; j < n0; ssj++,j++) {
      if ((h =   h     - qr) > /* gap open from left best */
	  /* gap extend from left gapped */
	  (f =   f     - r)) f = h;	/* if better, use new gap opened */
      if ((h = ssj->H - qr) >	/* gap open from up best */
	  /* gap extend from up gap */
	  (e = ssj->E - r)) e = h;	/* if better, use new gap opened */
      h = p + pam2p[j][*aa1p];
      /* h = p + *pwaa++; */		/* diagonal match */
      if (h < f ) h = f;	/* left gap better, reset */
      if (h < e ) h = e;	/* up gap better, reset */
      p = ssj->H;		/* save previous best score */
      ssj->H = h;		/* save (new) up diag-matched */
      ssj->E = e;		/* save upper gap opened */
    }
#ifndef GLOBAL_GLOBAL
    if (h > score) {		/* ? new best score at the end of each row */
      score = h;		/* save best */
      I = i;			/* row */
    }
#endif
    /*
    fprintf(stderr," r %d - score: %d ssj[]: %d\n", i,score,
	    ss[(i <= n0) ? i-1 : n0-1].H);
    */
    aa1p++;	/* aa1p goes down the path graph, row by row */
    i++;	/* increment the row */
  }		/* done with forward pass */

#ifdef GLOBAL_GLOBAL
  cost = score = h;
  K = 0;
  I = n1 - 1;
#else
  /*   fprintf(stderr, " r: %d - score: %d\n", I, score); */

  /* to get the start point, go backwards */
  
  cost = -BIGNUM;
  K = 0;
  ss[n0].H = 0;
  t = -q;
  for (ssj=ss+n0-1; ssj>=ss; ssj--) {
    ssj->H = t = t - r;
    ssj->E= t - q;
  }

  t = 0;
  for (i=I; i>=0; i--) {
    p = ss[n0].H;
    ss[n0].H = h = t = t-r;
    f = t-q;
    for (ssj=ss+J, j= J-1; j>=0; ssj--, j--) {
      if ((h =   h     - qr) > /* gap open from left best */
	  /* gap extend from left gapped */
	  (f =   f     - r)) f = h;	/* if better, use new gap opened */
      if ((h = ssj->H - qr) >	/* gap open from up best */
	  /* gap extend from up gap */
	  (e = ssj->E - r)) e = h;	/* if better, use new gap opened */
      h = p + pam2p[j][aa1[i]];		/* diagonal match */
      if (h < f ) h = f;	/* left gap better, reset */
      if (h < e ) h = e;	/* up gap better, reset */
      p = ssj->H;		/* save previous best score */
      ssj->H = h;		/* save (new) up diag-matched */
      ssj->E = e;		/* save upper gap opened */
    }
    if (h > cost) {
      cost = h;
      K = i;
      if (cost >= score) goto found;
    }
  }
  /* at this point, ss[0].E has a very high value for good alignments */
 found:
#endif /* not GLOBAL_GLOBAL */

/*   fprintf(stderr," *** %d=%d: L: %3d-%3d/%3d; K: %3d-%3d/%3d\n",score,cost,L,J+1,n0,K,I+1,n1); */

/* in the f_str version, the *res array is already allocated at 4*n0/3 */

  a_res->n1 = n1;
  a_res->max0 = J+1; a_res->min0 = L; a_res->max1 = I+1; a_res->min1 = K;
  
/* this code no longer refers to aa0[], it uses pam2p[0][L] instead */
  NW_ALIGN(L,&aa1[K-1],J-L+1,I-K+1,pam2p,q,r,a_res->res,&a_res->nres);

/*  DISPLAY(&aa0[L-1],&aa1[K-1],J-L+1,I-K+1,res,L,K,ppst->sq); */

/* return *res and nres */

  return score;
}

#define gap(k)  ((k) <= 0 ? 0 : g+h*(k))	/* k-symbol indel cost */

static int
FGLOBAL_ALIGN(int *waa, const unsigned char *aa1,
	      int n0, int n1,
	      int q, int r,
	      struct swstr *ss)
{
  const unsigned char *aa1p;
  register int *pwaa;
  int i,j;
  struct swstr *ssj;
  int t;
  int e, f, h, p;
  int qr;
  int score;
  int ij, max_col, max_ij;

  /* q - gap open is positve */
  /* r - gap extend is positive */

  qr = q+r;

  score = -BIGNUM;

  /* initialize 0th row */
  ss[0].H = 0;
  ss[0].E = t = -q;
  for (ssj = ss+1; ssj <= ss+n0 ; ssj++) {
    ssj->H = t = t - r;
    ssj->E = t - q;
  }

  aa1p = aa1;
  t = -q;
  while (*aa1p) {
    p = ss[0].H;
#if defined(GLOBAL_GLOBAL)
    ss[0].H = h = t = t - r;
#else		/* GLOBAL_LOCAL */
    ss[0].H = h = t = 0;
#endif
    f = t - q;
    pwaa = waa + (*aa1p++ * n0);
    for (ssj = ss+1; ssj <= ss+n0; ssj++) {	/* go across query */
      if ((h =   h    - qr) > (f =   f    - r)) f = h;
      if ((h = ssj->H - qr) > (e = ssj->E - r)) e = h;
      h = p + *pwaa++;
      if (h < f) h = f;
      if (h < e) h = e;
      p = ssj->H;
      ssj->H = h;
      ssj->E = e;
    }
#if !defined(GLOBAL_GLOBAL)	/* GLOBAL_LOCAL */
    if (h > score) {
      score = h;	/* at end of query, update score */
    }
#endif
  }				/* done with forward pass */
#ifdef GLOBAL_GLOBAL
  score = h;
#endif
  return score;
}

void do_opt (const unsigned char *aa0, int n0,
	     const unsigned char *aa1, int n1,
	     int frame,
	     struct pstruct *ppst, struct f_struct *f_str,
	     struct rstruct *rst)
{
}

/* this do_walign simply calls nsw_malign using nw_walign it is
   modeled after the same do_walign code in dropgsw2.c 

   It makes no sense to use this strategy for GLOBAL_GLOBAL
   alignments, but hopefully they will take care of themselves.
*/

struct a_res_str *
do_walign (const unsigned char *aa0, int n0,
	   const unsigned char *aa1, int n1,
	   int frame, int repeat_thresh,
	   struct pstruct *ppst, 
	   struct f_struct *f_str, 
	   int *have_ares)
{
  struct a_res_str *a_res, *tmp_a_res;
  int a_res_index;

  *have_ares = 0x3;	/* set 0x2 bit to indicate local copy */

  if ((a_res = (struct a_res_str *)calloc(1, sizeof(struct a_res_str)))==NULL) {
    fprintf(stderr," [do_walign] Cannot allocate a_res");
    return NULL;
  }
  
  a_res = nsw_malign(f_str->pam2p, (ppst->ext_sq_set ? 1 : 0), n0, aa1, n1,
		     repeat_thresh, f_str->max_res,
		     -ppst->gdelval, -ppst->ggapval,
		     f_str->ss, a_res,
		     &nw_walign, ppst->do_rep
		     );

  a_res_index = 0;
  for (tmp_a_res=a_res; tmp_a_res; tmp_a_res = tmp_a_res->next) {
    tmp_a_res->index = a_res_index++;
  }

  return a_res;
}

void
pre_cons(const unsigned char *aa1, int n1, int frame, struct f_struct *f_str) {

#ifdef TFAST
  f_str->n10 = aatran(aa1,f_str->aa1x,n1,frame);
#endif

}

/* aln_func_vals - set up aln.qlfact, qlrev, llfact, llmult, frame, llrev */
/* call from calcons, calc_id, calc_code */
void 
aln_func_vals(int frame, struct a_struct *aln) {

  aln->llfact = aln->llmult = aln->qlfact = 1;
  aln->llrev = 0;
  if (frame > 0) aln->qlrev = 1;
  else aln->qlrev = 0;
  aln->frame = 0;
}
