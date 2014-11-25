/* $Id: dropgsw2.c $ */

/* copyright (c) 1996, 2014 by William R. Pearson and The Rector &
   Visitors of the University of Virginia */

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

static char *verstr="7.2 Nov 2010";

#include "dropgsw2.h"

#define DROP_INTERN
#include "drop_func.h"

#ifdef SW_ALTIVEC
#include "smith_waterman_altivec.h"
#endif
#ifdef SW_SSE2
#include "smith_waterman_sse2.h"
#endif

struct swstr {int H, E;};

extern void init_karlin(const unsigned char *aa0, int n0, struct pstruct *ppst,
			double *aa0_f, double **kp);
extern int do_karlin(const unsigned char *aa1, int n1,
		     int **pam2, const struct pstruct *ppst,
		     double *aa0_f, double *kar_p, double *lambda, double *H);

extern int sw_walign (int **pam2p, int n0,
		      const unsigned char *aa1, int n1,
		      int q, int r,
		      struct swstr *ss,
		      struct a_res_str *a_res
		      );

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

void
SIM(const unsigned char *A, /* seq1 indexed A[1..M] */
    const unsigned char *B, /* seq2 indexed B[1..N] */
    int M, int N,		/* len seq1, seq2 */
    struct pstruct *ppst,	/* parameters */
    int nseq,			/* nseq - number of different sequences */
    int mini_score,		/* cut-off score */
    int max_count,		/* number of alignments */
    struct a_res_str *a_res);	/* alignment result structure */

static int
FLOCAL_ALIGN(const unsigned char *aa0, const unsigned char *aa1,
	     int n0, int n1, int low, int up,
	     int **W, int GG,int HH, int MW,
	     struct f_struct *f_str);

extern void aancpy(char *to, char *from, int count, struct pstruct *ppst);

int same_seq(const unsigned char *aa0, int n0,
	     const unsigned char *aa1, int n1);

static int
prof_score(const unsigned char *aa1, int n0, int *pwaa_s);

/* initialize for Smith-Waterman optimal score */

void
init_work (unsigned char *aa0, int n0,
	   struct pstruct *ppst,
	   struct f_struct **f_arg)
{
  int ip;
  int *pwaa_s, *pwaa_a;
  int e, f, i, j;
  struct f_struct *f_str;
  int **pam2p;
  struct swstr *ss;
  int nsq;

#if defined(SW_ALTIVEC) || defined(SW_SSE2)
  int l, data, bias;
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
    /* set for lower-case for memory mapped DBs with lower case encoding */
    nsq = ppst->nsqx; ip = 0;
  }

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

   ss[n0].H = -1;	/* this is used as a sentinel - normally H >= 0 */
   ss[n0].E = 1;
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
       pwaa[n0..2n0-1] has pam scores for amino acid 1 (A)
       pwaa[2n0..3n0-1] has pam scores for amino acid 2 (R), ...

       thus: pwaa = f_str->waa_s + (*aa1p++)*n0; sets up pwaa so that
       *pwaa++ rapidly moves though the scores of the aa1p[] position
       without further indexing

       For a real sequence profile, pwaa[0..n0-1] vs ['A'] could have
       a different score in each position.
   */

   pwaa_s = f_str->waa_s;
   pwaa_a = f_str->waa_a;
   if (ppst->pam_pssm) {
     for (e = 0; e <=nsq; e++)	{	/* for each residue in the alphabet */
       for (f = 0; f < n0; f++) {	/* for each position in aa0 */
	 *pwaa_s++ = f_str->pam2p[ip][f][e] = ppst->pam2p[ip][f][e];
	 *pwaa_a++ = f_str->pam2p[0][f][e]  = ppst->pam2p[0][f][e];
       }
     }
   }
   else {	/* initialize scanning matrix */
     for (e = 0; e <=nsq; e++)	/* for each residue in the alphabet */
       for (f = 0; f < n0; f++)	{	/* for each position in aa0 */
	 *pwaa_s++ = f_str->pam2p[ip][f][e]= ppst->pam2[ip][aa0[f]][e];
	 *pwaa_a++ = f_str->pam2p[0][f][e] = ppst->pam2[0][aa0[f]][e];
       }
   }

#if defined(SW_ALTIVEC)

   /* First we allocate memory for the workspace - i.e. the single row
    * of storage for H/F. Since this might be run on Linux or AIX too,
    * we don't assume anything about the memory allocation but align
    * it ourselves.  We need two vectors (16 bytes each) per element,
    * and some padding space to make it cache-line aligned.

    * MAXTST+MAXLIB is longest allowed database sequence length...
    * this should be m_msg.max_tot, but m_msg is not available, but
    * ppst->maxlen has maxn, which is appropriate.
    */

     f_str->workspace_memory  = (void *)malloc(2*16*(ppst->maxlen+SEQ_PAD)+256);
     f_str->workspace  = (void *) ((((size_t) f_str->workspace_memory) + 255) & (~0xff));
   /* We always use a scoring profile in altivec, but the layout is a bit strange 
    * in order to optimize memory access order and thus cache efficiency.
    * Normally we first try 8-bit scoring in altivec, and if this leads to overflow
    * we recompute the score with 16-bit accuracy. Because of this we need to construct
    * two score profiles.
    * Since altivec always loads 16 bytes from aligned memory, corresponding to 8 or 16 
    * elements (for 16 and 8 bit scoring, respectively), we organize the scoring 
    * profile like this for 8-bit accuracy:
    *
    * 1. The profile starts on 256-byte aligned memory (cache line on G5 is 128 bytes).
    * 2. First we have the score for the full alphabet for the first 16 residues of
    *    the query, i.e. positions 0-15 are the scores for the first 16 query letters
    *    vs. the first in the alphabet, positions 16-31 the scores for the same 16
    *    query positions against alphabet letter two, etc.
    * 3. After alphabet_size*16bytes we start with the scores for residues 16-31 in
    *    the query, organized in the same way.
    * 4. At the end of the query sequence, we pad the scoring to the next 16-tuple
    *    with neutral scores.
    * 5. The total size of the profile is thus alphabet_size*N, where N is the 
    *    size of the query rounded up to the next 16-tuple.
    *
    * The word (16-bit) profile is identical, but scores are stored as 8-tuples.
    */

   f_str->word_score_memory = (void *)malloc(10*2*(nsq+2)*(n0+1+16)+256);
   f_str->byte_score_memory = (void *)malloc(10*(nsq+2)*(n0+1+16)+256);

   f_str->word_score = (unsigned short *) ((((size_t) f_str->word_score_memory) + 255) & (~0xff));
   f_str->byte_score = (unsigned char *) ((((size_t) f_str->byte_score_memory) + 255) & (~0xff));

   overflow = 0;

   if (ppst->pam_pssm) {
     /* Use a position-specific scoring profile. 
      * This is essentially what we are going to construct anyway, but we'll
      * reorder it to suit altivec.
      */       
     bias = 127;
     for(i = 1; i < nsq ; i++) {
	 for(j = 0; j < n0 ; j++) {
	     data = ppst->pam2p[ip][j][i];
	     if(data<bias) bias = data;
           }
       }

     /* Fill our specially organized byte- and word-size scoring arrays. */
     ps = f_str->word_score;
     for(f = 0; f<n0 ; f+=8) {
       /* e=0 */
       for(i=0 ; i<8 ; i++) {
	 *ps++ = (unsigned short) 0;
       }
       /* for each chunk of 8 residues in our query */
       for(e = 1; e<=nsq; e++) {
	 for(i=0 ; i<8 ; i++) {
	   l = f + i;
	   if(l<n0) {
	     data = ppst->pam2p[ip][l][e] - bias;
	   }
	   else {
	     data = 0;
	   }
	   *ps++ = (unsigned short)data;
	 }
       }
     }
     pc = f_str->byte_score;
     for(f = 0; f<n0 ; f+=16) {
       /* e=0 */
       for(i=0 ; i<16 ; i++) {
	 *pc++ = (unsigned char)0;
       }       
           
       for(e = 1; e<=nsq; e++) {
	 for(i=0 ; i<16 ; i++) {
	   l = f + i;
	   if(l<n0) {
	     data = ppst->pam2p[ip][l][e] - bias;
	   }
	   else {
	     data = 0;
	   }
	   if(data>255) {
	     /*
	     printf("Fatal error. data: %d bias: %d, position: %d/%d, Score out of range for 8-bit Altivec/VMX datatype.\n",data,bias,l,e);
	     exit(1);
	     */
	     overflow = 1;
	   }
	   *pc++ = (unsigned char)data;
	 }
       }
     }
   }
   else {
     /* Classical simple substitution matrix */
     /* Find the bias to use in the substitution matrix */
     bias = 127;
     for(i = 1; i < nsq ; i++) {
       for(j = 1; j < nsq ; j++) {
	 data = ppst->pam2[ip][i][j];
#ifdef DEBUG
	 if (data < -1000) {
	   fprintf(stderr, "*** low data: %d [%d][%d][%d]\n",data,ip,i,j);
	 }
#endif
	 if(data<bias) bias = data;
       }
     }
     /* Fill our specially organized byte- and word-size scoring arrays. */
     ps = f_str->word_score;
     for(f = 0; f<n0 ; f+=8) {
       /* e=0 */
       for(i=0 ; i<8 ; i++) {
	 *ps++ = (unsigned short) 0;
       }       
       /* for each chunk of 8 residues in our query */
       for(e = 1; e<=nsq; e++) {
	 for(i=0 ; i<8 ; i++) {
	   l = f + i;
	   if(l<n0) {
	     data = ppst->pam2[ip][aa0[l]][e] - bias;
	   }
	   else {
	     data = 0;
	   }
	   *ps++ = (unsigned short)data;
	 }
       }
     }
     pc = f_str->byte_score;
     for(f = 0; f<n0 ; f+=16) {
       /* e=0 */
       for(i=0 ; i<16 ; i++) {
	 *pc++ = (unsigned char)0;
       }
           
       for(e = 1; e<=nsq; e++) {
	 for(i=0 ; i<16 ; i++) {
	   l = f + i;
	   if (l<n0) {
	     data = ppst->pam2[ip][aa0[l]][e] - bias;
	   }
	   else {
	     data = 0;
	   }
	   if(data>255) {
	     /*
	     printf("Fatal error. Score out of range for 8-bit Altivec/VMX datatype.\n");
	     exit(1);
	     */
	     overflow = 1;
	   }
	   *pc++ = (unsigned char)data;
	 }
       }
     }
   }
       
   f_str->bias = (unsigned char) (-bias);
   f_str->alphabet_size = nsq+1;

   /* Some variable to keep track of how many 8-bit runs we need to rerun
    * in 16-bit accuracy. If there are too many reruns it can be faster
    * to use 16-bit alignments directly. 
    */
   
   /* We can only do 8-bit alignments if the scores were small enough. */
   if(overflow==0) f_str->try_8bit   = 1;
   else f_str->try_8bit   = 0;

   f_str->done_8bit  = 0;
   f_str->done_16bit = 0;
       
#endif /* SW_ALTIVEC */

#if defined(SW_SSE2)
   /* First we allocate memory for the workspace - i.e. two rows for H and
    * one row for F.  We also need enough space to hold a temporary
    * scoring profile which will be query_length * 16 (sse2 word length).
    * Since this might be run on Linux or AIX too, we don't assume 
    * anything about the memory allocation but align it ourselves.
    */
    f_str->workspace_memory  = (void *)malloc(3*16*(MAXTST+MAXLIB+32)+256);
    f_str->workspace  = (void *) ((((size_t) f_str->workspace_memory) + 255) & (~0xff));

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

    f_str->word_score_memory = (void *)malloc((n0 + 32) * sizeof (short) * (nsq + 1) + 256);
    f_str->byte_score_memory = (void *)malloc((n0 + 32) * sizeof (char) * (nsq + 1) + 256);

    f_str->word_score = (unsigned short *) ((((size_t) f_str->word_score_memory) + 255) & (~0xff));
    f_str->byte_score = (unsigned char *) ((((size_t) f_str->byte_score_memory) + 255) & (~0xff));

    overflow = 0;

    if (ppst->pam_pssm) {
        /* Use a position-specific scoring profile. 
        * This is essentially what we are going to construct anyway, but we'll
        * reorder it to suit sse2.
        */       
        bias = 127;
        for (i = 1; i < nsq ; i++) {
            for (j = 0; j < n0 ; j++) {
                data = ppst->pam2p[ip][j][i];
                if (data < bias) {
                    bias = data;
                }
            }
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
		  if ( i < n0) { data = ppst->pam2p[ip][i][f];}
		  else {data = 0;}
		  *ps++ = (unsigned short)data;
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
		  if ( i < n0 ) { data = ppst->pam2p[ip][i][f] - bias;}
		  else {data = 0 - bias;}
		  if (data > 255) {
		    printf("Fatal error. data: %d bias: %d, position: %d/%d, "
			   "Score out of range for 8-bit SSE2 datatype.\n",
			   data, bias, f, e);
		    exit(1);
		  }
		  *pc++ = (unsigned char)data;
		}
	    }
        }
    }
    else 
    {
        /* Classical simple substitution matrix */
        /* Find the bias to use in the substitution matrix */
        bias = 127;
        for (i = 1; i < nsq ; i++) {
            for (j = 1; j < nsq ; j++) {
                data = ppst->pam2[ip][i][j];
		if (data < -128) {
		  fprintf(stderr,"*** ERROR *** data out of range: %d[%d][%d,%d]\n",
			  data, ip, i, j);
		}
                if (data < bias) {
                    bias = data;
                }
            }
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
                        data = ppst->pam2[ip][aa0[i]][f];
                    }
                    *ps++ = (unsigned short)data;
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
                        data = -bias;
                    } else {
                        data = ppst->pam2[ip][aa0[i]][f] - bias;
                    }
                    if (data > 255) {
                        printf("Fatal error. data: %d bias: %d, position: %d/%d, "
                               "Score out of range for 8-bit SSE2 datatype.\n",
                               data, bias, f, e);
                        exit(1);
                    }
                    *pc++ = (unsigned char)data;
                }
            }
        }
    }
       
    f_str->bias = (unsigned char) (-bias);
    f_str->alphabet_size = nsq+1;

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
    f_str->max_res =  max(3*n0/2,MIN_RES);

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
/*void get_param(struct pstruct *ppst,char *pstring1)*/
void    get_param (const struct pstruct *ppst,
		   char **pstring1, char *pstring2,
		   struct score_count_s *s_cnt_info)
{
  char pg_str[120];
  char psi_str[120];

#if defined(SW_ALTIVEC)
  strncpy(pg_str,"Smith-Waterman (Altivec/VMX, Erik Lindahl 2004)",sizeof(pg_str));
#endif
#if defined(SW_SSE2)
  strncpy(pg_str,"Smith-Waterman (SSE2, Michael Farrar 2006)",sizeof(pg_str));
#endif
#if !defined(SW_ALTIVEC) && !defined(SW_SSE2)
  strncpy(pg_str,"Smith-Waterman (PGopt)",sizeof(pg_str));
#endif

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
     sprintf(pstring2,"; pg_name_alg: %s\n; pg_ver_rel: %s\n; pg_matrix: %s%s (%d:%d)%s\n; pg_gap-pen: %d %d\n",
#else
     sprintf(pstring2,"; pg_name_alg: %s\n; pg_ver_rel: %s\n; pg_matrix: %s%s (%d:%d)%s\n; pg_open-ext: %d %d\n",
#endif
	     pg_str,verstr,ppst->pam_name,psi_str,ppst->pam_h,ppst->pam_l, 
	     (ppst->ext_sq_set)?"xS":"\0",ppst->gdelval,ppst->ggapval);
   }
}

void do_work (const unsigned char *aa0, int n0,
	      const unsigned char *aa1, int n1,
	      int frame,
	      const struct pstruct *ppst, struct f_struct *f_str,
	      int qr_flg, int shuff_flg, struct rstruct *rst,
	      struct score_count_s *sc_info)
{
  int     score;
  double lambda, H;
  int i;
  
#ifdef LALIGN
  if (same_seq(aa0,n0,aa1,n1)) {
    rst->score[0] = prof_score(aa1, n0, f_str->waa_s);
    return;
  }
#endif

  rst->alg_info = 0;
  rst->valid_stat = 1;
  sc_info->s_cnt[0]++;
  sc_info->tot_scores++;

#ifdef SW_ALTIVEC
  if(f_str->try_8bit)
  {
      score = smith_waterman_altivec_byte(aa0,
                                          f_str->byte_score,
                                          n0,
                                          aa1,
                                          n1,
                                          f_str->bias,
#ifndef OLD_FASTA_GAP
                                          -(ppst->gdelval + ppst->ggapval),
#else
                                          -ppst->gdelval,
#endif
                                          -ppst->ggapval,
                                          f_str);
      
      f_str->done_8bit++;
      
      if(score>=255)
      {
          /* Overflow, so we have to redo it in 16 bits. */
          score = smith_waterman_altivec_word(aa0,
                                              f_str->word_score,
                                              n0,
                                              aa1,
                                              n1,
                                              f_str->bias,
#ifndef OLD_FASTA_GAP
                                              -(ppst->gdelval + ppst->ggapval),
#else
                                              -ppst->gdelval,
#endif
                                              -ppst->ggapval,
                                              f_str);
          
          /* The 8 bit version is roughly 50% faster than the 16 bit version,
           * so we are fine if less than about 1/3 of the runs have to
           * be rerun with 16 bits. If it is more, and we have tried at least
           * 500 sequences, we switch off the 8-bit mode.
           */
          f_str->done_16bit++;
          if(f_str->done_8bit>500 && (3*f_str->done_16bit)>(f_str->done_8bit))
              f_str->try_8bit = 0;
      }
  }
  else
  { 
      /* Just use the 16-bit altivec version directly */
      score = smith_waterman_altivec_word(aa0,
                                          f_str->word_score,
                                          n0,
                                          aa1,
                                          n1,
                                          f_str->bias,
#ifndef OLD_FASTA_GAP
                                          -(ppst->gdelval + ppst->ggapval),
#else
                                          -ppst->gdelval,
#endif
                                          -ppst->ggapval,
                                          f_str);
  }      

#endif /* not Altivec */

#if defined(SW_SSE2)

  if(f_str->try_8bit)
  {
      score = smith_waterman_sse2_byte(aa0,
                                       f_str->byte_score,
                                       n0,
                                       aa1,
                                       n1,
                                       f_str->bias,
#ifndef OLD_FASTA_GAP
                                       -(ppst->gdelval + ppst->ggapval),
#else
                                       -ppst->gdelval,
#endif
                                       -ppst->ggapval,
                                       f_str);
      
      f_str->done_8bit++;
      
      if(score>=255)
      {
          /* Overflow, so we have to redo it in 16 bits. */
          score = smith_waterman_sse2_word(aa0,
                                           f_str->word_score,
                                           n0,
                                           aa1,
                                           n1,
#ifndef OLD_FASTA_GAP
                                           -(ppst->gdelval + ppst->ggapval),
#else
                                           -ppst->gdelval,
#endif
                                           -ppst->ggapval,
                                           f_str);
          
          /* The 8 bit version is roughly 50% faster than the 16 bit version,
           * so we are fine if less than about 1/3 of the runs have to
           * be rerun with 16 bits. If it is more, and we have tried at least
           * 500 sequences, we switch off the 8-bit mode.
           */
          f_str->done_16bit++;
          if(f_str->done_8bit>500 && (3*f_str->done_16bit)>(f_str->done_8bit))
              f_str->try_8bit = 0;
      }
  }
  else
  { 
      /* Just use the 16-bit altivec version directly */
      score = smith_waterman_sse2_word(aa0,
                                       f_str->word_score,
                                       n0,
                                       aa1,
                                       n1,
#ifndef OLD_FASTA_GAP
                                       -(ppst->gdelval + ppst->ggapval),
#else
                                       -ppst->gdelval,
#endif
                                       -ppst->ggapval,
                                       f_str);
  }      
#endif

#if !defined(SW_ALTIVEC) && !defined(SW_SSE2)

  score = FLOCAL_ALIGN(aa0,aa1,n0,n1,0,0,
                       NULL,
#ifdef OLD_FASTA_GAP
                       -(ppst->gdelval - ppst->ggapval),
#else
                       -ppst->gdelval,
#endif
                       -ppst->ggapval,0,f_str);
#endif

  rst->score[0] = score;
  rst->score[1] = rst->score[2] = 0;

  if(((ppst->zsflag % 10) == 6) &&
     (do_karlin(aa1, n1, ppst->pam2[0], ppst,f_str->aa0_f, 
		f_str->kar_p, &lambda, &H)>0)) {
    rst->comp = 1.0/lambda;
    rst->H = H;
  }
  else {rst->comp = rst->H = -1.0;}

}

static int
FLOCAL_ALIGN(const unsigned char *aa0, const unsigned char *aa1,
	     int n0, int n1, int low, int up,
	     int **W, int GG,int HH, int MW,
	     struct f_struct *f_str) {

  register int *pwaa;
  register struct swstr *ssj;
  struct swstr *ss;
  register int h, e, f, p;
  int temp, score;
  int gap_ext, n_gap_init;

  const unsigned char *aa1p;
  ss = f_str->ss;
  ss[n0].H = -1;
  ss[n0].E = 1;

  n_gap_init = GG + HH;
  gap_ext = -HH;	/* GG, HH are both positive,
			   gap_ext penalty should be negative */

  score = 0;
  for (h=0; h<n0; h++) {	  /* initialize 0th row */
    ss[h].H = ss[h].E = 0;
  }
  
  aa1p=aa1;
  while (*aa1p) {		/* relies on aa1[n1]==0 for EOS flag */
    /* waa_s has the offsets for each residue in aa0 into pam2 */
    /* waa_s has complexity (-S) dependent scores */
    pwaa = f_str->waa_s + (*aa1p++)*n0;
    ssj = ss;

    e = f = h = p = 0;
  zero_f:	/* in this section left-gap f==0, and is never examined */

    while (1) {	/* build until h > n_gap_init (f < 0 until h > n_gap_init) */
      		/* bump through the pam[][]'s for each of the aa1[] matches to
	  	   aa0[], because of the way *pwaa is set up */

      h = p + *pwaa++;		/* increment diag value */
      p = ssj->H;		/* get next diag value */
      if ((e = ssj->E) > 0 ) {	/* >0 from up-gap */
	if (p == -1) goto next_row;	/* done, -1=ss[n0].H sentinel */
	if (h < e) h = e;	/* up-gap better than diag */
	else 
	  if (h > n_gap_init) {	/* we won't starting a new up-gap */
	    e += gap_ext;	/* but we might be extending one */
	    goto transition;	/* good h > n_gap_diag; scan f */
	  }
	e += gap_ext;		/* up-gap decreased */
	ssj->E =  (e > 0) ?  e : 0;	/* set to 0 if < 0 */
	ssj++->H = h;		/* diag match updated */
      }
      else {			/* up-gap (->E) is 0 */
	if ( h > 0) {		/* diag > 0 */
	  if (h > n_gap_init) {	/* we won't be starting a new up-gap */
	    e = 0;		/* and we won't be extending one */
	    goto transition;	/* good h > n_gap_diag; scan f */
	  }
	  ssj++->H = h;		/* update diag */
	}
	else ssj++->H = 0;	/* update diag to 0 */
      }
    }

    /* here h > n_gap_init and h > e, => the next f will be > 0 */
  transition:
#ifdef DEBUG
    if ( h > 10000) 
      fprintf(stderr,"h: %d ssj: %d\n",h, (int)(ssj-ss));
#endif
    if ( score < h ) score = h;	/* save best score, only when h > n_gap_init */

    temp = h - n_gap_init;	/* best score for starting a new gap */
    if ( f < temp ) f = temp;	/* start a left-gap? */
    if ( e < temp ) e = temp;	/* start an up-gap? */
    ssj->E = ( e > 0 ) ? e : 0;	/* update up-gap */
    ssj++->H = h;		/* update diag */
    e = 0;

    do {			/* stay here until f <= 0 */
      h = p + *pwaa++;		/* diag + match/mismatch */
      p = ssj->H;		/* save next (right) diag */

      if ( h < f ) h = f;	/* update diag using left gap */
      f += gap_ext;		/* update next left-gap */

      if ((e = ssj->E) > 0) {	/* good up gap */
	if (p == -1) goto next_row;	/* at the end of the row */
	if ( h < e ) h = e;	/* update diag using up-gap */
	else
	  if ( h > n_gap_init ) {
	    e += gap_ext;	/* update up gap */
	    goto transition;	/* good diag > n_gap_init, restart */
	  }
	e += gap_ext;		/* update up-gap */
	ssj->E = (e > 0) ? e : 0;	/* e must be >= 0 */
	ssj++->H = h;		/* update diag */
      }
      else {			/* up-gap <= 0 */
	if ( h > n_gap_init ) {
	  e = 0;
	  goto transition;	/* good diag > n_gap_init; restart */
	}
	ssj++->H = h;		/* update diag */
      }
    } while ( f > 0 );		/* while left gap f > 0  */
    goto zero_f;		/* otherwise, go to f==0 section */
  next_row:
    ;
  }		/* end while(*aap1) {} */

  return score;

}		/* here we should be all done */

void do_opt (const unsigned char *aa0, int n0,
	     const unsigned char *aa1, int n1,
	     int frame,
	     struct pstruct *ppst, struct f_struct *f_str,
	     struct rstruct *rst)
{
}

struct a_res_str *
do_walign (const unsigned char *aa0, int n0,
	   const unsigned char *aa1, int n1,
	   int frame, int repeat_thresh,
	   struct pstruct *ppst, 
	   struct f_struct *f_str, 
	   int *have_ares)
{
  int a_res_index;
  struct a_res_str *a_res, *tmp_a_res;

  *have_ares = 0x3;	/* set 0x2 bit to indicate local copy */

  if ((a_res = (struct a_res_str *)calloc(1, sizeof(struct a_res_str)))==NULL) {
    fprintf(stderr," [do_walign] Cannot allocate a_res");
    return NULL;
  }

#ifndef LALIGN
  a_res = nsw_malign(f_str->pam2p, (ppst->ext_sq_set ? 1 : 0), n0, aa1, n1,
		     repeat_thresh, f_str->max_res,
		     -ppst->gdelval, -ppst->ggapval,
		     f_str->ss, a_res,
		     &sw_walign, ppst->do_rep);

#else	/* LALIGN */
  if (!ppst->show_ident && same_seq(aa0, n0, aa1, n1)) ppst->nseq = 1;
  else ppst->nseq = 2;

  SIM(aa0-1, aa1-1, n0, n1, ppst, ppst->nseq, repeat_thresh, ppst->max_repeat, a_res);
#endif

  /* set a_res->index for alignments */

  a_res_index = 0;
  for (tmp_a_res=a_res; tmp_a_res; tmp_a_res = tmp_a_res->next) {
    tmp_a_res->index = a_res_index++;
  }

  return a_res;
}

/*
#define XTERNAL
#include "upam.h"

void
print_seq_prof(unsigned char *A, int M,
	       unsigned char *B, int N,
	       int **w, int iw, int dir) {
  char c_max;
  int i_max, j_max, i,j;

  char *c_dir="LRlr";

  for (i=1; i<=min(60,M); i++) {
    fprintf(stderr,"%c",aa[A[i]]);
  }
  fprintf(stderr, - %d\n,M);

  for (i=0; i<min(60,M); i++) {
    i_max = -1;
    for (j=1; j<21; j++) {
      if (w[iw+i][j]> i_max) {
	i_max = w[iw+i][j]; 
	j_max = j;
      }
    }
    fprintf(stderr,"%c",aa[j_max]);
  }
  fputc(':',stderr);

  for (i=1; i<=min(60,N); i++) {
    fprintf(stderr,"%c",aa[B[i]]);
  }

  fprintf(stderr," -%c: %d,%d\n",c_dir[dir],M,N);
}
*/

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

/* calculate the 100% identical score */
int
prof_score(const unsigned char *aa1p, int n0, int *pwaa_s)
{
  int sum=0;

  while (*aa1p) {
    sum += pwaa_s[(*aa1p++)*n0];
    pwaa_s++;
  }
  return sum;
}

int same_seq(const unsigned char *aa0, int n0,
	     const unsigned char *aa1, int n1)
{
  const unsigned char *ap0, *ap1;
  int cnt=0;

  if (n0 != n1) return 0;

  ap0 = aa0;
  ap1 = aa1;
  
  while ( *ap0 && *ap0++ == *ap1++ ) {cnt++;}
  if (cnt != n0) return 0;
  return 1;
}
