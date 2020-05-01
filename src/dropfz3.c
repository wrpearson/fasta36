/* $Id: dropfz2.c 1280 2014-08-21 00:47:55Z wrp $ */

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

/* 18-Sept-2006 - removed static global variables for alignment */

/* 2002/06/23 finally correctly implement fix to translate 'N' to 'X' */

/* 1999/11/29 modification by Z. Zhang to translate DNA 'N' as 'X' */

/* implements an improved version of the fasty algorithm, see:

   W. R. Pearson, T. Wood, Z. Zhang, A W. Miller (1997) "Comparison of
   DNA sequences with protein sequences" Genomics 46:24-36

*/
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <ctype.h>

#include "defs.h"
#include "param.h"
#define XTERNAL
#include "upam.h"
#include "uascii.h"

#define NT_N 16

/* globals for fasta */
#define MAXWINDOW 64

#ifndef MAXSAV
#define MAXSAV 10
#endif

#ifndef ALLOCN0
static char *verstr="3.8 June 2014";
#else
static char *verstr="3.8an0 Jul 2014";
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

struct update_code_str {
  int p_op_idx;
  int p_op_cnt;
  int btop_enc;
  int show_code;
  int cigar_order;
  int show_ext;
  char *op_map;
};

#ifdef TFAST
static char *ori_code = "-x/=\\+*";	/* FASTX */
static char *cigar_code = "DXFMRI*";
#else
static char *ori_code = "+x/=\\-*";	/* TFASTX */
static char *cigar_code = "IXFMRD*";
#endif

unsigned long adler32(unsigned long, const unsigned char *, unsigned int);
void kpsort (struct savestr **v, int n);
extern void *init_stack(int, int);
extern void push_stack(void *, void *);
extern void *pop_stack(void *);
extern void *free_stack(void *);
extern struct domfeat_data * init_domfeat_data(const struct annot_str *annot_p);

#define SGW1 100
#define SGW2 300
struct smgl_str {
  int C[SGW1+1][SGW2+1];
  int st[SGW1+1][SGW2+1];
  int D[SGW2+7], I[SGW2+1];
};

struct sx_s {int C1, C2, C3, I1, I2, I3, flag; };

struct wgt { int  iii, ii, iv;};
struct wgtc {char c2, c3, c4, c5;};

typedef struct st_s { int C, I, D;} *st_ptr;

struct f_struct {
  struct dstruct *diag;
  int frame;
  int ndo;
  int noff;
  int hmask;			/* hash constants */
  int *pamh1;			/* pam based array */
  int *pamh2;			/* pam based kfact array */
  int *link, *harr;		/* hash arrays */
  int kshft;			/* shift width */
  int nsav, lowscor;		/* number of saved runs, worst saved run */
#ifndef TFAST
  unsigned char *aa0x, *aa0v;	/* aa0x - 111122223333 */
#else
  unsigned char *aa1x, *aa1v;	/* aa1x - 111122223333 */
#endif				/* aa1v - computed codons */
  struct sx_s *cur;
  int cur_sp_size;
  struct wgt **weight0;
  struct wgt **weight1;
  struct wgtc **weight_c;
  int *waa;
  int *res;
  int max_res;
  st_ptr up, down, tp;
  struct smgl_str smgl_s;
};

#define DROP_INTERN
#include "drop_func.h"

static int dmatchz(const unsigned char *aa0, int n0,
		   const unsigned char *aa1, int n1,
		   const unsigned char *aa1v,
		   int hoff, int window, 
		   int **pam2, int gdelval, int ggapval, int gshift,
		   struct f_struct *f_str);

int shscore(unsigned char *aa0, int n0, int **pam2);
int saatran(const unsigned char *ntseq, unsigned char *aaseq, int maxs, int frame);
extern int ELK_to_s(double E_join, int n0, int n1, double Lambda, double K, double H);

int savemax (struct dstruct *, int, 
	     struct savestr *vmax, struct savestr **lowmax);

int spam (const unsigned char *aa0, const unsigned char *aa1, 
	  struct savestr *dmax, int **pam2,
	  struct f_struct *f_str);
int sconn (struct savestr **v, int n,int cgap, int pgap, struct f_struct *f_str);
int lx_band(const unsigned char *prot_seq, int len_prot,
	    const unsigned char *dna_prot_seq, int len_dna_prot,
	    int **pam_matrix, int gopen, int gext,
	    int gshift, int start_diag, int width, struct f_struct *f_str);

struct a_res_str *
merge_ares_chains(struct a_res_str *cur_ares, 
		  struct a_res_str *tmpl_ares,
		  int score_ix, const char *msg);

static struct update_code_str *
init_update_data(int show_code);

static void
sprintf_code(char *tmp_str, struct update_code_str *, int op_idx, int op_cnt);

static void 
update_code(struct dyn_string_str *align_code_dyn,
	    struct update_code_str *update_data, int op, 
	    int sim_code, unsigned char sp0, unsigned char sp1);

static void
close_update_data(struct dyn_string_str *align_code_dyn,
		  struct update_code_str *update_data);

extern void w_abort (char *p, char *p1);
extern void aagetmap(char *to, int n);

/* initialize for fasta */
/* modified 30-August-1999 by Zheng Zhang to work with an extended alphabet */
/* Assume naa=47, and wgts[47][23] matches both upper and lower case
amoino acids with another amino acid.  And also assume the DNA letter
does not have upper/lower case difference.  If you also allow DNA
sequence to be upper/lower case letters, more needs be changed. Not
only here, but also in the alignment code, the way that pack a codon
into a number between 0-63 need be changed. */

/* modified so that if **weightci==NULL, do not fiddle with characters */

/* modified 3-Aug-2010 for NCBIstdaa alphabet, which requires MAXUC
   28, MAXLC 56, so we must have 58, not 47, entries */

void
init_weights(struct wgt ***weighti, struct wgtc ***weightci,
	     int **wgts, int gshift, int gsubs, int naa)
{
  int i, j, do_wgtc=0;
  int aa, b, a, x, y, z;
  int *wwt, e;
  struct wgt **weight;
  struct wgtc **weightc;
  char aacmap[64];
  int temp[MAXLC+1][64]; /*change*/
  char le[MAXLC+1][64];

  if (naa > MAXLC) {
    fprintf(stderr,"*** ERROR [%s:%d] - compilation problem naa(%d) > MAXLC(%d) ***\n",
	    __FILE__, __LINE__, naa, MAXLC);
  }

  if ((*weighti=(struct wgt **)calloc((size_t)(naa+1),sizeof(struct wgt *)))
      ==NULL) {
    fprintf(stderr,"*** ERROR [%s:%d] - cannot allocate weights array: %d\n",
	    __FILE__, __LINE__, naa);
    exit(1);
  }

  weight = *weighti;
  /* allocate weight[aa 0..MAXLC] */
  for (aa=0; aa <= naa; aa++) {
    if ((weight[aa]=(struct wgt *)calloc((size_t)256,sizeof(struct wgt)))
	==NULL) {
      fprintf(stderr,"*** ERROR [%s:%d] - cannot allocate weight[]: %d/%d\n",
	      __FILE__, __LINE__, aa,naa);
      exit(1);
    }
  }

  /* allocate weightci[aa 0..MAXLC] */
  if (weightci !=NULL) {
    if ((*weightci=(struct wgtc **)calloc((size_t)(naa+1),
					  sizeof(struct wgtc *)))==NULL) {
      fprintf(stderr,"*** ERROR [%s:%d] -  cannot allocate weight_c array: %d\n",
	      __FILE__, __LINE__, naa);
      exit(1);
    }
    weightc = *weightci;

    for (aa=0; aa <= naa; aa++) {
      if ((weightc[aa]=(struct wgtc *)calloc((size_t)256,sizeof(struct wgtc)))
	  ==NULL) {
	fprintf(stderr,"*** ERROR [%s:%d] -  cannot allocate weightc[]: %d/%d\n",
		__FILE__, __LINE__, aa,naa);
	exit(1);
      }
    }
    do_wgtc = 1;
  }
  else do_wgtc = 0;

  aagetmap(aacmap,64);

  for (aa = 0; aa < naa; aa++) {	/* change*/
    wwt = wgts[aa];			/* pam matrix */
    for (i = 0; i < 64; i++) {	/* i iterates through the codons */
      x = -10000;		/* large negative */
      y = i;
      for (j = 0; j < 64; j++) {	/* j iterates through the codons */
	z = ((~i & j) | (i & ~j));
	b = 0;		/* score = 0 */
	if (z % 4) b-= gsubs;
	if (z /16) b-= gsubs;
	if ((z /4) % 4) b -= gsubs;   
	b += wwt[aascii[aacmap[j]]];  /* add the match score for char j*/
	if (b > x) {
	  x = b;		/* x has the score */
	  y = j;		/* y has the character  (codon index)*/
	}
      }
#ifdef DEBUG
      if (y < 0 || y > 63) printf("%d %d %d %d ",aa, i, x, y);
#endif
      temp[aa][i] = x;
      le[aa][i] = y;
    }
    /*            printf("\n"); */
  }

  for (aa= 0; aa < naa; aa++) { 
      wwt = temp[aa];
      for (i = 0; i < 256; i++) {
          for (x=-100,b = 0; b < 4; b++) {
              z = (i/ (1 << ((b+1)*2)))*(1<<(b*2))+(i%(1<<(b*2)));
	      if (x < (e=wwt[z])) {
		  x = e;
		  if (do_wgtc) weightc[aa][i].c4 = aacmap[le[aa][z]];
	      }
          }
          weight[aa][i].iv=x-gshift;
          weight[aa][i].iii = wwt[i%64];

	  if (do_wgtc) {
	    weightc[aa][i].c5 = aacmap[le[aa][i%64]];
	    weightc[aa][i].c3 = aacmap[i%64];
	  }
          x = i %16;
          for (y = -100, b = 0; b < 3; b++) {
              z = ((x >> (b*2)) << (b*2+2)) + (x % (1 << (b*2))); 
              for (a = 0; a < 4; a++) {
		  if ((e =wwt[z+(a<<(b*2))]) > y) {
		      y = e;
		      if (do_wgtc) 
			weightc[aa][i].c2 = aacmap[le[aa][z+(a<<(b*2))]];
		  }
              }
          }
          weight[aa][i].ii = y-gshift;
      }
  }
  /*106=CGGG*/
  for (aa = 0; aa < naa; aa++) {
    weight[aa][106].iii = wgts[aa][23]; /* is 23 the code for 'X'?*/
    weight[aa][106].iv = weight[aa][106].ii = weight[aa][106].iii-gshift;
    if (do_wgtc) {
      weightc[aa][106].c5 = weightc[aa][106].c4 = weightc[aa][106].c3
	= weightc[aa][106].c2 = 'X';
    }
  }
}

void
free_weights(struct wgt ***weighti0, struct wgt ***weighti1, 
	     struct wgtc ***weightci, int naa)
{
  int aa;
  struct wgt **weight0;
  struct wgt **weight1;
  struct wgtc **weightc;

  weight0 = *weighti0;
  weight1 = *weighti1;
  weightc = *weightci;

  for (aa=0; aa < naa; aa++) {free(weight0[aa]);}
  for (aa=0; aa < naa; aa++) {free(weight1[aa]);}
  for (aa=0; aa < naa; aa++) {free(weightc[aa]);}

  free(weight0);
  free(weight1);
  free(weightc);
}

static void
pre_com(const unsigned char *aa0, int n0, unsigned char *aa0v) {
  int dnav, i;
  dnav = (hnt[aa0[0]]<<2) + hnt[aa0[1]];
  for (i=2; i<n0; i++) {
    dnav = ((dnav<<2)+hnt[aa0[i]])&255;
    if (aa0[i] == NT_N  || aa0[i-1]==NT_N || aa0[i-2] == NT_N)  {
      aa0v[i-2] = 106;
    }
    else {
      if (dnav == 106/*CGGG*/) {dnav = 42/*AGGG*/;}
      aa0v[i-2]=dnav;
    }
  }
}

static void
pre_com_r(const unsigned char *aa0, int n0, unsigned char *aa0v) {
  int dnav, i, ir;
  dnav = ((3-hnt[aa0[n0-1]])<<2) + 3-hnt[aa0[n0-2]];
  for (i=2, ir=n0-3; i<n0; i++,ir--) {
    dnav = ((dnav<<2)+3-hnt[aa0[ir]])&255; 
    if (aa0[ir] == NT_N || aa0[ir+1]==NT_N || aa0[ir+2] == NT_N)  {
      aa0v[i-2] = 106;
    }
    else {
      if (dnav == 106) dnav = 42;
      aa0v[i-2]=dnav;
    }
  }
}

void
init_work (unsigned char *aa0, int n0, 
	   struct pstruct *ppst,
	   struct f_struct **f_arg)
{
   int mhv, phv;
   int hmax;
   int i0, hv;
   int pamfact;
   int btemp;
   struct f_struct *f_str;
   struct bdstr *bss;
   /* these used to be globals, but do not need to be */
   int ktup, fact, kt1, lkt;

   int maxn0;
   int *pwaa;
   int i, j, q;
   struct swstr *ss, *r_ss;
   int *waa;
   int *res;
   int nsq, ip, *hsq, naat;
#ifndef TFAST
   int last_n0, itemp, dnav;
   unsigned char *fd, *fs, *aa0x, *aa0v;
   int n0x, n0x3;
#endif

   if (nt[NT_N] != 'N') {
     fprintf(stderr,"*** ERROR [%s:%d] -  nt[NT_N] (%d) != 'X' (%c) - recompile\n",
	     __FILE__, __LINE__, NT_N,nt[NT_N]);
     exit(1);
   }     

   if (ppst->ext_sq_set) {
     nsq = ppst->nsqx; ip = 1;
     hsq = ppst->hsqx;
   }
   else {
     nsq = ppst->nsqx; ip = 0;
     hsq = ppst->hsq;
   }

   f_str = (struct f_struct *)calloc(1,sizeof(struct f_struct));

   if (!ppst->param_u.fa.use_E_thresholds) {
     btemp = 2 * ppst->param_u.fa.bestoff / 3 +
       n0 / ppst->param_u.fa.bestscale +
       ppst->param_u.fa.bkfact *
       (ppst->param_u.fa.bktup - ppst->param_u.fa.ktup);
     btemp = min (btemp, ppst->param_u.fa.bestmax);
     if (btemp > 3 * n0) btemp = 3 * shscore(aa0,n0,ppst->pam2[0]) / 5;

     ppst->param_u.fa.cgap = btemp + ppst->param_u.fa.bestoff / 3;
     if (ppst->param_u.fa.optcut_set != 1) {
#ifndef TFAST
       ppst->param_u.fa.optcut = (btemp*5)/4;
#else
       ppst->param_u.fa.optcut = (btemp*4)/3;
#endif
     }
   }

#ifdef OLD_FASTA_GAP
   ppst->param_u.fa.pgap = ppst->gdelval + ppst->ggapval;
#else
   ppst->param_u.fa.pgap = ppst->gdelval + 2*ppst->ggapval;
#endif
   pamfact = ppst->param_u.fa.pamfact;
   ktup = ppst->param_u.fa.ktup;
   fact = ppst->param_u.fa.scfact * ktup;

#ifndef TFAST
   /* before hashing, we must set up some space and translate the sequence */

   maxn0 = n0 + 2;
   if ((aa0x =(unsigned char *)calloc((size_t)maxn0,
					     sizeof(unsigned char)))
       == NULL) {
     fprintf (stderr, "*** ERROR [%s:%d] - cannot allocate aa0x array %d\n",
	      __FILE__, __LINE__, maxn0);
     exit (1);
   }
   aa0x++;
   f_str->aa0x = aa0x;


   if ((aa0v =(unsigned char *)calloc((size_t)maxn0,
					     sizeof(unsigned char)))
       == NULL) {
     fprintf (stderr, "*** ERROR [%s:%d] - cannot allocate aa0v array %d\n",
	      __FILE__, __LINE__, maxn0);
     exit (1);
   }
   aa0v++;
   f_str->aa0v = aa0v;

   /* make a precomputed codon number series */
   pre_com(aa0, n0, aa0v);

   last_n0 = 0;
   for (itemp=0; itemp<3; itemp++) {
     n0x=saatran(aa0,&aa0x[last_n0],n0,itemp);
     /*         for (i=0; i<n0x; i++) {
	   fprintf(stderr,"%c",aa[aa0x[last_n0+i]]);
	   if ((i%60)==59) fprintf(stderr,"\n");
	   }
	   fprintf(stderr,"\n");
	   */
     last_n0 += n0x+1;
   }

   /*     fprintf(stderr,"\n"); */
   n0x = n0;
   n0x3 = n0x/3;

   /* now switch aa0 and aa0x for hashing functions */
   fs = aa0;
   aa0 = aa0x;
   aa0x = fs;
#endif

   /* naat must always be MAXLC because library can have LC aa residues */
   /*
   if (ppst->ext_sq_set) naat = MAXLC;
   else naat = MAXUC;
   */
   naat = MAXLC;

   init_weights(&f_str->weight0, NULL,
		ppst->pam2[ip],-ppst->gshift,-ppst->gsubs,naat);
   init_weights(&f_str->weight1, &f_str->weight_c,
		ppst->pam2[0],-ppst->gshift,-ppst->gsubs,naat);

   if (pamfact == -1)
      pamfact = 0;
   else if (pamfact == -2)
      pamfact = 1;

   for (i0 = 1, mhv = -1; i0 < ppst->nsq; i0++)
      if (hsq[i0] < NMAP && hsq[i0] > mhv)
	 mhv = ppst->hsq[i0];

   if (mhv <= 0) {
      fprintf (stderr, "*** ERROR [%s:%d] -  maximum hsq <=0 %d\n",
	       __FILE__, __LINE__, mhv);
      exit (1);
   }

   for (f_str->kshft = 0; mhv > 0; mhv /= 2) f_str->kshft++;

/*      kshft = 2;	*/
   kt1 = ktup - 1;
   hv = 1;
   for (i0 = 0; i0 < ktup; i0++)
      hv = hv << f_str->kshft;
   hmax = hv;
   f_str->hmask = (hmax >> f_str->kshft) - 1;

   if ((f_str->harr = (int *) calloc (hmax, sizeof (int))) == NULL) {
     fprintf (stderr, "*** ERROR [%s:%d] - cannot allocate hash array [%d]\n",
	      __FILE__, __LINE__, hmax);
     exit (1);
   }
   if ((f_str->pamh1 = (int *) calloc (ppst->nsq+1, sizeof (int))) == NULL) {
     fprintf (stderr, "*** ERROR [%s:%d] - cannot allocate pamh1 array [%d]\n",
	      __FILE__, __LINE__, ppst->nsq+1);
     exit (1);
   }
   if ((f_str->pamh2 = (int *)calloc (hmax, sizeof (int))) == NULL) {
     fprintf (stderr, "*** ERROR [%s:%d] -  cannot allocate pamh2 array [%d]\n",
	      __FILE__, __LINE__, hmax);
     exit (1);
   }
   if ((f_str->link = (int *) calloc (n0, sizeof (int))) == NULL) {
     fprintf (stderr, "*** ERROR [%s:%d] - cannot allocate hash link array [%d]",
	      __FILE__, __LINE__, n0);
     exit (1);
   }

   for (i0 = 0; i0 < hmax; i0++)
      f_str->harr[i0] = -1;
   for (i0 = 0; i0 < n0; i0++)
      f_str->link[i0] = -1;

   /* encode the aa0 array */
   phv = hv = 0;
   lkt = kt1;
   for (i0 = 0; i0 < min(n0,lkt); i0++) {
     if (hsq[aa0[i0]] >= NMAP) {
       hv=phv=0; lkt = i0+ktup; continue;
     }
     hv = (hv << f_str->kshft) + ppst->hsq[aa0[i0]];
     phv += ppst->pam2[ip][aa0[i0]][aa0[i0]] * ktup;
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
     hv = ((hv & f_str->hmask) << f_str->kshft) + ppst->hsq[aa0[i0]];
     f_str->link[i0] = f_str->harr[hv];
     f_str->harr[hv] = i0;
     if (pamfact) {
       f_str->pamh2[hv] = (phv += ppst->pam2[ip][aa0[i0]][aa0[i0]] * ktup);
       if (hsq[aa0[i0-kt1]] < NMAP) 
	 phv -= ppst->pam2[ip][aa0[i0 - kt1]][aa0[i0 - kt1]] * ktup;
     }
     else f_str->pamh2[hv] = fact * ktup;
   }

/* this has been modified from 0..<ppst->nsq to 1..<=ppst->nsq because the
   pam2[0][0] is now undefined for consistency with blast
*/

   if (pamfact)
      for (i0 = 1; i0 < ppst->nsq; i0++)
	 f_str->pamh1[i0] = ppst->pam2[ip][i0][i0] * ktup;
   else
      for (i0 = 1; i0 < ppst->nsq; i0++)
	 f_str->pamh1[i0] = fact;

   f_str->ndo = 0;	/* used to save time on diagonals with long queries */


#ifndef ALLOCN0
   if ((f_str->diag = (struct dstruct *) calloc ((size_t)MAXDIAG,
						 sizeof (struct dstruct)))==NULL) {
      fprintf (stderr,"*** ERROR [%s:%d] - cannot allocate diagonal arrays: %lu\n",
	       __FILE__, __LINE__, MAXDIAG *sizeof (struct dstruct));
      exit (1);
     };
#else
   if ((f_str->diag = (struct dstruct *) calloc ((size_t)n0,
					      sizeof (struct dstruct)))==NULL) {
      fprintf (stderr,"*** ERROR [%s:%d] - cannot allocate diagonal arrays: %ld\n",
	       __FILE__, __LINE__, (long)n0*sizeof (struct dstruct));
      exit (1);
     };
#endif

#ifndef TFAST
   /* done hashing, now switch aa0, aa0x back */
   fs = aa0;
   aa0 = aa0x;
   aa0x = fs;
#else
   if ((f_str->aa1x =(unsigned char *)calloc((size_t)ppst->maxlen+4,
					     sizeof(unsigned char)))
       == NULL) {
     fprintf (stderr, "*** ERROR [%s:%d] - cannot allocate aa1x array %d\n",
	      __FILE__, __LINE__, ppst->maxlen+4);
     exit (1);
   }
   f_str->aa1x++;

   if ((f_str->aa1v =(unsigned char *)calloc((size_t)ppst->maxlen+4,
					     sizeof(unsigned char))) == NULL) {
     fprintf (stderr, "*** ERROR [%s:%d] - cannot allocate aa1v array %d\n",
	      __FILE__, __LINE__, ppst->maxlen+4);
     exit (1);
   }
   f_str->aa1v++;

#endif

   if ((waa= (int *)malloc (sizeof(int)*(nsq+1)*n0)) == NULL) {
     fprintf(stderr,"*** ERROR [%s:%d] - cannot allocate waa struct %3d\n",
	     __FILE__, __LINE__, nsq*n0);
     exit(1);
   }

   pwaa = waa;
   for (i=0; i<nsq; i++) {
     for (j=0;j<n0; j++) {
       *pwaa = ppst->pam2[ip][i][aa0[j]];
       pwaa++;
     }
   }
   f_str->waa = waa;

#ifndef TFAST
   maxn0 = max(2*n0,MIN_RES);
#else
   maxn0 = max(4*n0,MIN_RES);
#endif
   if ((res = (int *)calloc((size_t)maxn0,sizeof(int)))==NULL) {
     fprintf(stderr,"*** ERROR [%s:%d] - cannot allocate alignment results array %d\n",
	     __FILE__, __LINE__, maxn0);
     exit(1);
   }
   f_str->res = res;
   f_str->max_res = maxn0;

   *f_arg = f_str;
}

/* pstring1 is a message to the manager, currently 512 */
/* pstring2 is the same information, but in a markx==10 format */
void
get_param (const struct pstruct *ppst,
	   char **pstring1, char *pstring2,
	   struct score_count_s *s_info)
{
  char options_str1[128];
  char options_str2[128];
#ifndef TFAST
  char *pg_str="FASTY";
#else
  char *pg_str="TFASTY";
#endif

  if (!ppst->param_u.fa.use_E_thresholds) {
    sprintf(options_str1,"join: %d (%.3g), opt: %d (%.3g)",
	    ppst->param_u.fa.cgap, (double)s_info->s_cnt[0]/(double)s_info->tot_scores,
	    ppst->param_u.fa.optcut, (double)s_info->s_cnt[2]/(double)s_info->tot_scores);
    sprintf(options_str2,"pg_join: %d (%.3g)\n; pg_optcut: %d (%.3g)",
	    ppst->param_u.fa.cgap, (double)s_info->s_cnt[0]/(double)s_info->tot_scores,
	    ppst->param_u.fa.optcut, (double)s_info->s_cnt[2]/(double)s_info->tot_scores);
  }
  else {
    sprintf(options_str1,"E-join: %.2g (%.3g), E-opt: %.2g (%.3g)",
	    ppst->param_u.fa.E_join, (double)s_info->s_cnt[0]/(double)s_info->tot_scores,
	    ppst->param_u.fa.E_band_opt, (double)s_info->s_cnt[2]/(double)s_info->tot_scores);
    sprintf(options_str2,"pg_join_E(): %.2g (%.3g)\n; pg_optcut_E(): %.2g (%.3g)",
	    ppst->param_u.fa.E_join, (double)s_info->s_cnt[0]/(double)s_info->tot_scores,
	    ppst->param_u.fa.E_band_opt, (double)s_info->s_cnt[2]/(double)s_info->tot_scores);
  }

  if (!ppst->param_u.fa.optflag)
    sprintf (pstring1[0], "%s (%s)",pg_str, verstr);
  else 
    sprintf (pstring1[0], "%s (%s) [optimized]",pg_str, verstr);

  sprintf (pstring1[1], 
#ifdef OLD_FASTA_GAP
	   "%s matrix (%d:%d)%s, gap-pen: %3d/%3d, shift: %3d, subs: %3d\n ktup: %d, %s, width: %3d",
#else
	   "%s matrix (%d:%d)%s, open/ext: %3d/%3d, shift: %3d, subs: %3d\n ktup: %d, %s, width: %3d",
#endif
	   ppst->pam_name, ppst->pam_h,ppst->pam_l,
	   (ppst->ext_sq_set) ? "xS":"\0",
	   ppst->gdelval, ppst->ggapval,
	   ppst->gshift,ppst->gsubs,
	   ppst->param_u.fa.ktup, options_str1, ppst->param_u.fa.optwid);

  if (ppst->param_u.fa.iniflag) strcat(pstring1[0]," init1");

  if (pstring2 != NULL) {
#ifdef OLD_FASTA_GAP
    sprintf (pstring2, "; pg_name_alg: %s\n; pg_ver_rel: %s\n; pg_matrix: %s (%d:%d)%s\n\
; pg_gap-pen: %d %d\n; pg_ktup: %d\n; %s\n",
#else
	     sprintf (pstring2, "; pg_name_alg: %s\n; pg_ver_rel: %s\n; pg_matrix: %s (%d:%d)%s\n\
; pg_open-ext: %d %d\n; pg_ktup: %d\n; %s\n",
#endif
		      pg_str,verstr,ppst->pam_name, ppst->pam_h,ppst->pam_l,
		      (ppst->ext_sq_set) ? "xS":"\0",  ppst->gdelval,
		      ppst->ggapval,ppst->param_u.fa.ktup,options_str2);
	     }
  }

void
close_work (const unsigned char *aa0, int n0,
	    struct pstruct *ppst,
	    struct f_struct **f_arg)
{
  struct f_struct *f_str;
  int naat;

  f_str = *f_arg;

  if (f_str != NULL) {
   if (ppst->ext_sq_set) naat = MAXLC;
   else naat = MAXUC;
    free_weights(&f_str->weight0,&f_str->weight1,&f_str->weight_c,naat);
    free(f_str->cur);
#ifndef TFAST
    f_str->aa0v--;
    free(f_str->aa0v);
    f_str->aa0x--;
    free(f_str->aa0x);
#else	/* TFAST */
    f_str->aa1x--;
    free(f_str->aa1x);
    f_str->aa1v--;
    free(f_str->aa1v);
#endif
    free(f_str->res);
    free(f_str->waa);
    free(f_str->diag);
    free(f_str->link);
    free(f_str->pamh2); 
    free(f_str->pamh1);
    free(f_str->harr);
    free(f_str);
    *f_arg = NULL;
  }
}

void do_fastz (const unsigned char *aa0, int n0,
	       const unsigned char *aa1, int n1,
	       const unsigned char *aa1v,
	       const struct pstruct *ppst, struct f_struct *f_str,
	       struct rstruct *rst, int *hoff, int shuff_flg,
	       struct score_count_s *s_info)
{
   int     nd;		/* diagonal array size */
   int     lhval;
   int     kfact;
   int i;
   register struct dstruct *dptr;
   struct savestr vmax[MAXSAV];	/* best matches saved for one sequence */
   struct savestr *vptr[MAXSAV];
   struct savestr *lowmax;
   int lowscor;
   register int tscor;
   int xdebug = 0;

#ifndef ALLOCN0
   register struct dstruct *diagp;
#else
   register int dpos;
   int     lposn0;
#endif
   struct dstruct *dpmax;
   register int lpos;
   int     tpos;
   struct savestr *vmptr;
   int     scor, tmp;
   int     im, ib, nsave;
   int ktup, kt1, ip, lkt, ktup_sq;
   const int *hsq;
   int c_gap, opt_cut;
#ifndef TFAST
   int n0x31, n0x32;
   n0x31 = (n0-2)/3;
   n0x32 = n0x31+1+(n0-n0x31-1)/2;
#else
   unsigned char *fs, *fd;
   int n1x31, n1x32, last_n1, itemp;
   n1x31 = (n1-2)/3;
   n1x32 = n1x31+1+(n1-n1x31-1)/2;
#endif

  if (ppst->ext_sq_set) {
    ip = 1;
    hsq = ppst->hsqx;
  }
  else {
    ip = 0;
    hsq = ppst->hsq;
  }

  ktup = ppst->param_u.fa.ktup;
  ktup_sq = ktup*ktup;
  if (ktup == 1) ktup_sq *= 2;

  kt1 = ktup-1;

  if (n1 < ktup) {
    rst->score[0] = rst->score[1] = rst->score[2] = 0;
    return;
  }

  if (n0+n1+1 >= MAXDIAG) {
    fprintf(stderr,"*** ERROR [%s:%d] - n0,n1 too large > %d: %d, %d\n",
	    __FILE__, __LINE__, n0,n1, MAXDIAG);
    rst->score[0] = rst->score[1] = rst->score[2] = -1;
    return;
  }

  if (ppst->param_u.fa.use_E_thresholds) {
    c_gap = ELK_to_s(ppst->param_u.fa.E_join*ktup_sq*2.5, n0, n1, ppst->pLambda, ppst->pK, ppst->pH);
    opt_cut = ELK_to_s(ppst->param_u.fa.E_band_opt*ktup_sq*2.0, n0, n1, ppst->pLambda, ppst->pK, ppst->pH);
    rst->valid_stat = 0;
  }
  else {
    c_gap = ppst->param_u.fa.cgap;
    opt_cut = ppst->param_u.fa.optcut;
    rst->valid_stat = 1;
  }

  f_str->noff = n0 - 1;

#ifdef ALLOCN0
   nd = n0;
#endif

#ifndef ALLOCN0
   nd = n0 + n1;
#endif

   dpmax = &f_str->diag[nd];
   for (dptr = &f_str->diag[f_str->ndo]; dptr < dpmax;)
   {
      dptr->stop = -1;
      dptr->dmax = NULL;
      dptr++->score = 0;
   }

   for (vmptr = vmax; vmptr < &vmax[MAXSAV]; vmptr++)
      vmptr->score = 0;
   lowmax = vmax;
   lowscor = 0;

  if (n1 > 1000 && aa1[0]==23 && aa1[100]==23 &&
      aa1[1400]==23 && aa1[1401]!=23) {
    xdebug = 1;
  }
  else xdebug = 0;

   /* start hashing */
   lhval = 0;
   lkt = kt1;
   for (lpos = 0; (lpos < lkt || hsq[aa1[lpos]]>=NMAP) && lpos<n1; lpos++) {
     /* restart lhval calculation */
     if (hsq[aa1[lpos]]>=NMAP) {
       lhval = 0; lkt=lpos+ktup;
       continue;
     }
     lhval = ((lhval & f_str->hmask) << f_str->kshft) + ppst->hsq[aa1[lpos]];
   }

#ifndef ALLOCN0
   diagp = &f_str->diag[f_str->noff + lkt];
   for (; lpos < n1; lpos++, diagp++) {
     if (hsq[aa1[lpos]]>=NMAP) {
       lpos++ ; diagp++;
       while (lpos < n1 && hsq[aa1[lpos]]>=NMAP) {lpos++; diagp++;}
       if (lpos >= n1) break;
       lhval = 0;
     }
     lhval = ((lhval & f_str->hmask) << f_str->kshft) + ppst->hsq[aa1[lpos]];
     for (tpos = f_str->harr[lhval]; tpos >= 0; tpos = f_str->link[tpos]) {
       if ((tscor = (dptr = &diagp[-tpos])->stop) >= 0) {
#else
   lposn0 = f_str->noff + lpos;
   for (; lpos < n1; lpos++, lposn0++) {
     if (hsq[aa1[lpos]]>=NMAP) {lhval = 0; goto loopl;}
     lhval = ((lhval & f_str->hmask) << f_str->kshft) + ppst->hsq[aa1[lpos]];
     for (tpos = f_str->harr[lhval]; tpos >= 0; tpos = f_str->link[tpos]) {
       dpos = lposn0 - tpos;
       if ((tscor = (dptr = &f_str->diag[dpos % nd])->stop) >= 0) {
#endif
	 tscor += ktup;
	 if ((tscor -= lpos) <= 0) {
	   scor = dptr->score;
	   if ((tscor += (kfact = f_str->pamh2[lhval])) < 0 && lowscor < scor) {
#ifdef ALLOCN0
	     lowscor = savemax (dptr, dpos, vmax, &lowmax);
#else
	     lowscor = savemax (dptr, dptr- f_str->diag, vmax, &lowmax);
#endif
	   }
	   if ((tscor += scor) >= kfact) {
	     dptr->score = tscor;
	     dptr->stop = lpos;
	   }
	   else {
	     dptr->score = kfact;
	     dptr->start = (dptr->stop = lpos) - kt1;
	   }
	 }
	 else {
	   dptr->score += f_str->pamh1[aa0[tpos]];
	   dptr->stop = lpos;
	 }
       }
       else {
	 dptr->score = f_str->pamh2[lhval];
	 dptr->start = (dptr->stop = lpos) - kt1;
       }
     }				/* end tpos */

#ifdef ALLOCN0
      /* reinitialize diag structure */
   loopl:
     if ((dptr = &f_str->diag[lpos % nd])->score > lowscor)
       lowscor = savemax (dptr, lpos, vmax, &lowmax);
     dptr->stop = -1;
     dptr->dmax = NULL;
     dptr->score = 0;
#endif
   }	/* end lpos */

#ifdef ALLOCN0
   for (tpos = 0, dpos = f_str->noff + n1 - 1; tpos < n0; tpos++, dpos--) {
     if ((dptr = &f_str->diag[dpos % nd])->score > lowscor)
       lowscor = savemax (dptr, dpos, vmax, &lowmax, f_str);
   }
#else
   for (dptr = f_str->diag; dptr < dpmax;) {
     if (dptr->score > lowscor) savemax (dptr, dptr - f_str->diag, vmax, &lowmax);
     dptr->stop = -1;
     dptr->dmax = NULL;
     dptr++->score = 0;
   }
   f_str->ndo = nd;
#endif

   for (nsave = 0, vmptr = vmax; vmptr < &vmax[MAXSAV]; vmptr++) {
     if (vmptr->score > 0) {
       vmptr->score = spam (aa0, aa1, vmptr, ppst->pam2[ip], f_str);
       vptr[nsave++] = vmptr;
     }
   }

   if (nsave <= 0) {
     rst->score[0] = rst->score[1] = rst->score[2] = 0;
     return;
   }
       
#ifndef TFAST
   /* FASTX code here to modify the start, stop points for 
      the three phases of the translated protein sequence
   */

   for (ib=0; ib<nsave; ib++) {
     if (f_str->noff-vptr[ib]->dp+vptr[ib]->start >= n0x32)
       vptr[ib]->dp += n0x32;
     if (f_str->noff-vptr[ib]->dp +vptr[ib]->start >= n0x31)
       vptr[ib]->dp += n0x31;
   }
#else
   /* TFAST code here to modify the start, stop points for 
      the three phases of the translated protein sequence
      TFAST modifies library start points, rather than 
      query start points
   */

   for (ib=0; ib<nsave; ib++) {
     if (vptr[ib]->start >= n1x32) {
       vptr[ib]->start -= n1x32;
       vptr[ib]->stop -= n1x32;
       vptr[ib]->dp -= n1x32;
     }
     if (vptr[ib]->start >= n1x31) {
       vptr[ib]->start -= n1x31;
       vptr[ib]->stop -= n1x31;
       vptr[ib]->dp -= n1x31;
     }
   }
#endif /* TFAST */

   scor = sconn (vptr, nsave, c_gap, 
		 ppst->param_u.fa.pgap, f_str);

   for (vmptr=vptr[0],ib=1; ib<nsave; ib++)
     if (vptr[ib]->score > vmptr->score) vmptr=vptr[ib];

/*  kssort (vptr, nsave); */

   rst->score[1] = vmptr->score;
   rst->score[0] = max (scor, vmptr->score);
   rst->score[2] = rst->score[0];		/* initn */

   s_info->tot_scores++;
   if (rst->score[0] > c_gap) { s_info->s_cnt[0]++;}
#ifndef TFAST
   *hoff=f_str->noff - vmptr->dp;
#else /* TFAST */
   *hoff=vmptr->dp-f_str->noff;
#endif /* TFAST */
   if (ppst->param_u.fa.optflag) {
     if (/* shuff_flg || */ rst->score[0] > opt_cut) {
       s_info->s_cnt[2]++;
       rst->valid_stat = 1;
       rst->score[2] = dmatchz(aa0, n0,aa1,n1, aa1v,
			       *hoff,ppst->param_u.fa.optwid,
			       ppst->pam2[ip],
			       ppst->gdelval,ppst->ggapval,ppst->gshift,
			       f_str);
     }
   }
}

void do_work (const unsigned char *aa0, int n0,
	      const unsigned char *aa1, int n1,
	      int frame,
	      const struct pstruct *ppst,
	      struct f_struct *f_str,
	      int qr_flg, int shuff_flg, struct rstruct *rst,
	      struct score_count_s *s_info)
{
  int hoff;
  int last_n1, itx, dnav, n10, i, ir;
  unsigned char *aa1x;

  rst->escore = 1.0;
  rst->segnum = rst->seglen = 1;
  rst->valid_stat = 0;

  if (n1 < ppst->param_u.fa.ktup) {
    rst->score[0] = rst->score[1] = rst->score[2] = 0;
    return;
  }

#ifndef TFAST
  do_fastz (f_str->aa0x, n0, aa1, n1, f_str->aa0v, ppst, f_str, rst, &hoff, shuff_flg, s_info);
#else
   /* make a precomputed codon number series */

  if (frame == 0) {
    pre_com(aa1, n1, f_str->aa1v);
  }
  else {
    pre_com_r(aa1, n1, f_str->aa1v);
  }

  /* make translated sequence */
  last_n1 = 0;
  aa1x = f_str->aa1x;
#ifdef DEBUG
  if (frame > 1) {
    fprintf(stderr, "*** ERROR [%s:%d] - fz_walign - frame: %d - out of range [0,1]\n",
	    __FILE__, __LINE__, frame);
  }
#endif

  for (itx= frame*3; itx< frame*3+3; itx++) {
    n10  = saatran(aa1,&aa1x[last_n1],n1,itx);
    last_n1 += n10+1;
  }
  n10 = last_n1-1;

  do_fastz (aa0, n0, f_str->aa1x, n10, f_str->aa1v, ppst, f_str, rst, &hoff, shuff_flg, s_info);
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
  int optflag, tscore, hoff;
  int last_n1, itx, n10, i, ir;
  unsigned char *aa1x;
  struct score_count_s s_info = {0, 0, 0, 0};

  optflag = ppst->param_u.fa.optflag;
  ppst->param_u.fa.optflag = 1;

#ifndef TFAST
  do_fastz (f_str->aa0x, n0, aa1, n1, f_str->aa0v, ppst, f_str, rst, &hoff, 0, &s_info);
#else
   /* make a precomputed codon number series */

  if (frame == 0) {
    pre_com(aa1, n1, f_str->aa1v);
  }
  else {
    pre_com_r(aa1, n1, f_str->aa1v);
  }

  /* make translated sequence */
  last_n1 = 0;
  aa1x = f_str->aa1x;
  for (itx= frame*3; itx< frame*3+3; itx++) {
    n10  = saatran(aa1,&aa1x[last_n1],n1,itx);
    last_n1 += n10+1;
  }
  n10 = last_n1-1;

  do_fastz (aa0, n0, f_str->aa1x, n10, f_str->aa1v, ppst, f_str, rst, &hoff, 0, &s_info );
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
      if (vmptr != *lowmax) return (*lowmax)->score;
   }
   else {
      i = (*lowmax)->score = dptr->score;
      (*lowmax)->dp = dpos;
      (*lowmax)->start = dptr->start;
      (*lowmax)->stop = dptr->stop;
      dptr->dmax = *lowmax;
   }

   for (vmptr = vmax; vmptr < vmax+MAXSAV; vmptr++) {
     if (vmptr->score < i) {
       i = vmptr->score;
       *lowmax = vmptr;
     }
   }
   return i;
}

int spam (const unsigned char *aa0,
	  const unsigned char *aa1,
	  struct savestr *dmax, int **pam2,
	  struct f_struct *f_str)
{
   int     lpos;
   int     tot, mtot;
   struct {
     int     start, stop, score;
   } curv, maxv;
   const unsigned char *aa0p, *aa1p;

   aa1p = &aa1[lpos = dmax->start];
   aa0p = &aa0[lpos - dmax->dp + f_str->noff];
   curv.start = lpos;

   tot = curv.score = maxv.score = 0;
   for (; lpos <= dmax->stop; lpos++) {
     tot += pam2[*aa0p++][*aa1p++];
     if (tot > curv.score) {
       curv.stop = lpos;
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

/*	if (maxv.start != dmax->start || maxv.stop != dmax->stop)
		printf(" new region: %3d %3d %3d %3d\n",maxv.start,
			dmax->start,maxv.stop,dmax->stop);
*/
   dmax->start = maxv.start;
   dmax->stop = maxv.stop;

   return maxv.score;
}

#define XFACT 10

int sconn (struct savestr **v, int n, 
       int cgap, int pgap, struct f_struct *f_str)
{
  int     i, si;
  struct slink {
    int     score;
    struct savestr *vp;
    struct slink *next;
  }      *start, *sl, *sj, *so, sarr[MAXSAV];
  int     lstart, tstart, plstop, ptstop;

/*	sort the score left to right in lib pos */

   kpsort (v, n);

   start = NULL;

/*	for the remaining runs, see if they fit */

   for (i = 0, si = 0; i < n; i++)
   {

/*	if the score is less than the gap penalty, it never helps */
      if (v[i]->score < cgap)
	 continue;
      lstart = v[i]->start;
      tstart = lstart - v[i]->dp + f_str->noff;

/*	put the run in the group */
      sarr[si].vp = v[i];
      sarr[si].score = v[i]->score;
      sarr[si].next = NULL;

/* 	if it fits, then increase the score */
      for (sl = start; sl != NULL; sl = sl->next)
      {
	 plstop = sl->vp->stop;
	 ptstop = plstop - sl->vp->dp + f_str->noff;
	 if (plstop < lstart+XFACT && ptstop < tstart+XFACT) {
	   sarr[si].score = sl->score + v[i]->score + pgap;
	   break;
	 }
      }

/*	now recalculate where the score fits */
      if (start == NULL)
	 start = &sarr[si];
      else
	 for (sj = start, so = NULL; sj != NULL; sj = sj->next)
	 {
	    if (sarr[si].score > sj->score)
	    {
	       sarr[si].next = sj;
	       if (so != NULL)
		  so->next = &sarr[si];
	       else
		  start = &sarr[si];
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

static int
dmatchz(const unsigned char *aa0, int n0,
	const unsigned char *aa1, int n1,
	const unsigned char *aa1v,
	int hoff, int window, 
	int **pam2, int gdelval, int ggapval, int gshift,
	struct f_struct *f_str)
{

   hoff -= window/2;

#ifndef TFAST
   return lx_band(aa1,n1,f_str->aa0v,n0-2, 
		  pam2,
#ifdef OLD_FASTA_GAP
		  -(gdelval - ggapval),
#else
		  -gdelval,
#endif
		  -ggapval,-gshift,
		  hoff,window,f_str);
#else
   return lx_band(aa0,n0,aa1v,n1-2, 
		  pam2,
#ifdef OLD_FASTA_GAP
		  -(gdelval - ggapval),
#else
		  -gdelval,
#endif
		  -ggapval,-gshift,
		  hoff,window,f_str);
#endif
}

static void
init_row(struct sx_s *row, int sp) {
  int i;
  for (i = 0; i < sp; i++) {
      row[i].C1 = row[i].I1 = 0;
      row[i].C2 = row[i].I2 = 0;
      row[i].C3 = row[i].I3 = 0;
      row[i].flag = 0;
  }
}

int lx_band(const unsigned char *prot_seq,  /* array with protein sequence numbers*/
	    int len_prot,    /* length of prot. seq */
	    const unsigned char *dna_prot_seq, /* translated DNA sequence numbers*/
	    int len_dna_prot,   /* length trans. seq. */
	    int **pam_matrix,   /* scoring matrix */
	    int gopen, int gext, /* gap open, gap extend penalties */
	    int gshift,         /* frame-shift penalty */
	    int start_diag,     /* start diagonal of band */
	    int width,         /* width for band alignment */
	    struct f_struct *f_str)
{
  void *ckalloc();
  int i, j, bd, bd1, x1, x2, sp, p1=0, p2=0, end_prot;
  struct sx_s *last, *tmp;
  int sc, del, best = 0, cd,ci, e1, e2, e3, cd1, cd2, cd3, f, gg;
  const unsigned char *dp;
  register struct sx_s *ap, *aq;
  struct wgt *wt, *ww;
  int aa, b, a,x,y,z;

  sp = width+7;
  gg = gopen+gext;
  /* sp = sp/3+1; */

  if (f_str->cur == NULL ) {
    f_str->cur_sp_size = sp;
    f_str->cur = (struct sx_s *) ckalloc(sizeof(struct sx_s)*sp);
  }
  else if (f_str->cur_sp_size != sp) {
    free(f_str->cur);
    f_str->cur = (struct sx_s *) ckalloc(sizeof(struct sx_s)*sp);
    f_str->cur_sp_size = sp;
  }

  init_row(f_str->cur, sp);

  /*
  if (start_diag %3 !=0) start_diag = start_diag/3-1;
  else start_diag = start_diag/3;
  if (width % 3 != 0) width = width/3+1;
  else width = width /3;
  */

  x1 = start_diag; 		/* x1 = lower bound of DNA */
  x2 = 1;               /* the amount of position shift from last row*/

  end_prot = max(0,-width-start_diag) + (len_dna_prot+5)/3 + width;
  end_prot = min(end_prot,len_prot);

  /* i counts through protein sequence, x1 through DNAp */

  for (i = max(0, -width-start_diag), x1+=i; i < len_prot; i++, x1++) {
      bd = min(x1+width, (len_dna_prot+2)/3);	/* upper bound of band */
      bd1 = max(0,x1);	                /* lower bound of band */
      wt = f_str->weight0[prot_seq[i]];
      del = 1-x1;   /*adjustment*/
      bd += del; 
      bd1 +=del;

      ap = &f_str->cur[bd1]; aq = ap+1;
      e1 = f_str->cur[bd1-1].C3; e2 = ap->C1; cd1 = cd2= cd3= 0;
      for (dp = &dna_prot_seq[(bd1-del)*3]; ap < &f_str->cur[bd]; ap++) {
	  ww = &wt[(unsigned char) *dp++];
	  sc = max(max(e1+ww->iv, (e3=ap->C2)+ww->ii), e2+ww->iii);
	  if (cd1 > sc) sc = cd1;
	  cd1 -= gext;
	  if ((ci = aq->I1) > 0) {
	      if (sc < ci) { ap->C1 = ci; ap->I1 = ci-gext;}
	      else {
		  ap->C1 = sc;
		  sc -= gg;
		  if (sc > 0) {
		      if (sc > best) best =sc;
		      if (cd1 < sc) cd1 = sc;
		      ap->I1 = max(ci-gext, sc);
		  } else ap->I1 = ci-gext;
	      }
	  } else {
	      if (sc <= 0) {
		  ap->I1 = ap->C1 = 0;
	      } else {
		  ap->C1 = sc; sc-=gg;
		  if (sc >0) {
		      if (sc > best) best =sc;
		      if (cd1 < sc) cd1 = sc;
		      ap->I1 = sc;
		  } else ap->I1 = 0;
	      }
	  }
	  ww = &wt[(unsigned char) *dp++];
	  sc = max(max(e2+ww->iv, (e1=ap->C3)+ww->ii), e3+ww->iii);
	  if (cd2 > sc) sc = cd2;
	  cd2 -= gext;
	  if ((ci = aq->I2) > 0) {
	      if (sc < ci) { ap->C2 = ci; ap->I2 = ci-gext;}
	      else {
		  ap->C2 = sc;
		  sc -= gg;
		  if (sc > 0) {
		      if (sc > best) best =sc;
		      if (cd2 < sc) cd2 = sc;
		      ap->I2 = max(ci-gext, sc);
		  }
	      }
	  } else {
	      if (sc <= 0) {
		  ap->I2 = ap->C2 = 0;
	      } else {
		  ap->C2 = sc; sc-=gg;
		  if (sc >0) {
		      if (sc > best) best =sc;
		      if (cd2 < sc) cd2 = sc;
		      ap->I2 = sc;
		  } else ap->I2 = 0;
	      }
	  }
	  ww = &wt[(unsigned char)*dp++];
	  sc = max(max(e3+ww->iv, (e2=aq->C1)+ww->ii), e1+ww->iii);
	  if (cd3 > sc) sc = cd3;
	  cd3 -= gext;
	  if ((ci = aq++->I3) > 0) {
	      if (sc < ci) { ap->C3 = ci; ap->I3 = ci-gext;}
	      else {
		  ap->C3 = sc;
		  sc -= gg;
		  if (sc > 0) {
		      if (sc > best) best =sc;
		      if (cd3 < sc) cd3 = sc;
		      ap->I3 = max(ci-gext, sc);
		  }
	      }
	  } else {
	      if (sc <= 0) {
		  ap->I3 = ap->C3 = 0;
	      } else {
		  ap->C3 = sc; sc-=gg;
		  if (sc >0) {
		      if (sc > best) best =sc;
		      if (cd3 < sc) cd3 = sc;
		      ap->I3 = sc;
		  } else ap->I3 = 0;
	      }
	  }
      }
  }
  /*  printf("The best score is %d\n", best); */
  return best+gg;
}

/* ckalloc - allocate space; check for success */
void *ckalloc(size_t amount)
{
	void *p;

	if ((p = (void *)malloc( (size_t)amount)) == NULL)
		w_abort("Ran out of memory.","");
	return(p);
}

/* calculate the 100% identical score */
int
shscore(unsigned char *aa0, int n0, int **pam2)
{
  int i, sum;
  for (i=0,sum=0; i<n0; i++)
    sum += pam2[aa0[i]][aa0[i]];
  return sum;
}

#define WIDTH 60

typedef struct mat *match_ptr;

typedef struct mat {
	int i, j, l;
	match_ptr next;
} match_node;

typedef struct { int i,j;} state;
typedef state *state_ptr;


void *ckalloc();
static match_ptr small_global(), global();
static int local_align(), find_best();
static void init_row2(), init_ROW();

int
pro_dna(const unsigned char *prot_seq,  /* array with prot. seq. numbers*/
	int len_prot,    /* length of prot. seq */
	const unsigned char *dna_prot_seq, /* trans. DNA seq. numbers*/
	int len_dna_prot,   /* length trans. seq. */
	int **pam_matrix,   /* scoring matrix */
	int gopen, int gext, /* gap open, gap extend penalties */
	int gshift,         /* frame-shift penalty */
	struct f_struct *f_str,
	int max_res,
	struct a_res_str *a_res) /* alignment info */
{
  match_ptr align, ap, aq;
  int x, y, ex, ey, i, score;
  int *alignment;

  f_str->up = (st_ptr) ckalloc(sizeof(struct st_s)*(len_dna_prot+10));
  f_str->down = (st_ptr) ckalloc(sizeof(struct st_s)*(len_dna_prot+10));
  f_str->tp = (st_ptr) ckalloc(sizeof(struct st_s)*(len_dna_prot+10));

  /*local alignment find the best local alignment x and y
    is the starting position of the best local alignment
    and ex ey is the ending position */

  score= local_align(&x, &y, &ex, &ey, 
		     pam_matrix, gopen, gext,
		     dna_prot_seq, len_dna_prot,
		     prot_seq, len_prot, f_str);

  f_str->up += 3; f_str->down += 3; f_str->tp += 3;

  /* x, y - start in prot, dna_prot */
  a_res->min0 = x;	/* prot */
  a_res->min1 = y;	/* DNA */
  a_res->max0 = ex;	/* prot */
  a_res->max1 = ey;	/* DNA */

  align = global(x, y, ex, ey, 
		 pam_matrix, gopen, gext,
		 dna_prot_seq, prot_seq,
		 0, 0, f_str);

  alignment = a_res->res;

  for (ap = align, i= 0; ap; i++) {
    if (i < max_res) alignment[i] = ap->l;
    aq = ap->next; free(ap); ap = aq;
  }
  if (i >= max_res)
    fprintf(stderr,"*** ERROR [%s:%d] - alignment truncated: %d >= %d***\n",
	    __FILE__, __LINE__, i, max_res);

  /*   up = &up[-3]; down = &down[-3]; tp = &tp[-3]; */
  free(&f_str->up[-3]); free(&f_str->tp[-3]); free(&f_str->down[-3]);

  a_res->nres = i;
  return score;
}

static void
swap(void **a, void **b)
{
    void *t = *a;
    *a = *b;   *b = t;
}

/*
   local alignment find the best local alignment x and y
   is the starting position of the best local alignment
   and ex ey is the ending position 
*/
static int
local_align(int *x, int *y, int *ex, int *ey, 
	    int **wgts, int gop, int gext,
	    const unsigned char *dnap, int ld,
	    const unsigned char *pro, int lp,
	    struct f_struct *f_str)
{
  int i, j,  score, x1,x2,x3,x4, e1 = 0, e2 = 0, e3,
    sc, del,  e, best = 0,  cd, ci, c;
  struct wgt *wt, *ww;
  state_ptr cur_st, last_st, cur_i_st;
  st_ptr cur, last;
  const unsigned char *dp;
  int *cur_d_st, *st_up;

  /*      
	  Array rowiC stores the best scores of alignment ending at a position
	  Arrays rowiD and rowiI store the best scores of alignment ending
	  at a position with a deletion or insrtion
	  Arrays sti stores the starting position of the best alignment whose
	  score stored in the corresponding row array.
	  The program stores two rows to complete the computation, same is
	  for the global alignment routine.
  */


  st_up = (int *) ckalloc(sizeof(int)*(ld+10));
  init_row2(st_up, ld+5);

  ld += 2;

  init_ROW(f_str->up, ld+1);
  init_ROW(f_str->down, ld+1);
  cur = f_str->up+1;
  last = f_str->down+1; 

  cur_st = (state_ptr) ckalloc(sizeof(state)*(ld+1));
  last_st = (state_ptr) ckalloc(sizeof(state)*(ld+1));
  cur_i_st = (state_ptr) ckalloc(sizeof(state)*(ld+1));
  cur_d_st = st_up; 
  dp = dnap-2;
  for (i = 0; i < lp; i++) {
    wt = f_str->weight1[pro[i]];  e2 =0; e1 = last[0].C;
    for (j = 0; j < 2; j++) {
      cur_st[j].i = i+1;
      cur_st[j].j = j+1;
    }
    for (j = 2; j < ld; j++) {
      ww = &wt[(unsigned char) dp[j]];
      del = -1;
      if (j >= 3) {
	sc = 0;
	e3 = e2; e2 = e1;
	e1 = last[j-2].C; 
	if ((e=e2+ww->iii) > sc) {sc = e; del = 3;}
	if ((e=e1+ww->ii) > sc) {sc = e; del = 2;}
	if ((e = e3+ww->iv) > sc) {sc = e; del = 4;} 
      } else {
	sc = e2  = 0;
	if (ww->iii > 0) {sc = ww->iii; del = 3;}
      }
      if (sc < (ci=last[j].I)) {
	sc = ci; del = 0;
      }
      if (sc < (cd=cur[j].D)) {
	sc = cd; del = 5;
      }
      cur[j].C = sc;
      e = sc  - gop;
      if (e > cd) {
	cur[j+3].D = e-gext;
	cur_d_st[j+3] = 3;
      } else {
	cur[j+3].D = cd-gext;
	cur_d_st[j+3] = cur_d_st[j]+3;
      }
      switch(del) {
      case 5:
	c = cur_d_st[j];
	cur_st[j].i = cur_st[j-c].i;
	cur_st[j].j = cur_st[j-c].j;
	break;
      case 0:
	cur_st[j].i = cur_i_st[j].i;
	cur_st[j].j = cur_i_st[j].j;
	break;
      case 2:
      case 3:
      case 4:
	if (i) {
	  if (j-del >= 0) {
	    cur_st[j].i = last_st[j-del].i;
	    cur_st[j].j = last_st[j-del].j;
	  } else {
	    cur_st[j].i = i;
	    cur_st[j].j = 0;
	  }
	} else {
	  cur_st[j].i = 0;
	  cur_st[j].j = max(0, j-del+1);
	}
	break;
      case -1:
	cur_st[j].i = i+1;
	cur_st[j].j = j+1;
	break;
      }
      if (e > ci) {
	cur[j].I  = e -gext;
	cur_i_st[j].i = cur_st[j].i;
	cur_i_st[j].j = cur_st[j].j;
      } else {
	cur[j].I  = ci- gext;
      }
      if (sc > best) {
	x1 = cur_st[j].i;
	x2 = cur_st[j].j;
	best =sc;
	x3 = i;
	x4 = j;
      }
    }
    swap((void *)&last, (void *)&cur);
    swap((void *)&cur_st, (void *)&last_st);
  }
  /*	printf("The best score is %d\n", best);*/
  *x = x1; *y = x2; *ex = x3; *ey = x4;
  free(cur_st); free(last_st); free(cur_i_st); 
  free(st_up);
  return best;
}

/* 
   Both global_up and global_down do linear space score only global 
   alignments on subsequence pro[x]...pro[ex], and dna[y]...dna[ey].
   global_up do the algorithm upwards, from row x towards row y.
   global_down do the algorithm downwards, from row y towards x.
*/

static void
global_up(st_ptr *row1, st_ptr *row2, 
	  int x, int y, int ex, int ey, 
	  int **wgts, int gop, int gext,
	  unsigned char *dnap, unsigned char *pro, 
	  int N, struct f_struct *f_str)
{
  int i, j, k, sc, e, e1, e2, e3, t, ci, cd, score;
  struct wgt *wt, *ww;
  st_ptr cur, last;

  cur = *row1; last = *row2;
  sc = -gop;
  for (j = 0; j <= ey-y+1; j++) {
    if (j % 3 == 0) {last[j].C = sc; sc -= gext; last[j].I = sc-gop;}
    else { last[j].I = last[j].C = -10000;}
  }  
  last[0].C = 0; cur[0].D = cur[1].D = cur[2].D = -10000;
  last[0].D = last[1].D = last[2].D = -10000;
  if (N) last[0].I = -gext;
  for (i = 1; i <= ex-x+1; i++) {
    wt = f_str->weight1[pro[i+x-1]]; e1 = -10000; e2 = last[0].C;
    for (j = 0; j <= ey-y+1; j++) {
      t = j+y;
      sc = -10000; 
      ww = &wt[(unsigned char) dnap[t-3]]; 
      if (j < 4) {
	if (j == 3) {
	  sc = e2+ww->iii;
	} else if (j == 2) {
	  sc = e2 + ww->ii;
	}
      } else {
	e3 = e2; e2 = e1;
	e1 = last[j-2].C;
	sc = max(e2+ww->iii, max(e1+ww->ii, e3+ww->iv));
      }
      sc = max(sc, max(ci=last[j].I, cd = cur[j].D));
      cur[j].C = sc;
      cur[j+3].D = max(cd, sc-gop)-gext;
      cur[j].I = max(ci, sc-gop)-gext;
    }
    swap((void *)&last, (void *)&cur);
  }
  /*printf("global up score =%d\n", last[ey-y+1].C);*/
  for (i = 0; i <= ey-y+1; i++) last[i].I = cur[i].I;	
  if (*row1 != last) swap((void *)row1, (void *)row2);
}

static void
global_down(st_ptr *row1, st_ptr *row2,
	    int x, int y, int ex, int ey, 
	    int **wgts, int gop, int gext,
	    unsigned char *dnap, unsigned char *pro,
	    int N, struct f_struct *f_str)
{
  int i, j, k, sc, del, *tmp, e,  t, e1,e2,e3, ci,cd, score;
  struct wgt *wt, *w1, *w2, *w3;
  st_ptr cur, last;

  cur = (*row1); last = *row2;
  sc = -gop;
  for (j = ey-y+1; j >= 0; j--) {
    if ((ey-y+1-j) % 3) {last[j].C = sc; sc-=gext; last[j].I = sc-gop;}
    else  last[j].I =  last[j].C = -10000;
    cur[j].I = -10000;
  } 
  last[ey-y+1].C = 0;
  if (N) last[ey-y+1].I = -gext;
  cur[ey-y+1].D = cur[ey-y].D = cur[ey-y-1].D = -10000;
  last[ey-y+1].D = last[ey-y].D = last[ey-y-1].D = -10000;
  for (i = ex-x; i >= 0; i--) {
    wt = f_str->weight1[pro[i+x]]; e2 = last[ey-y+1].C; 
    e1 = -10000;
    w3 = &wt[(unsigned char) dnap[ey]]; 
    w2 = &wt[(unsigned char) dnap[ey-1]];
    for (j = ey-y+1; j >= 0; j--) {
      t = j+y;
      w1 = &wt[(unsigned char) dnap[t-1]];
      sc = -10000;
      if (t+3 > ey) {
	if (t+2 == ey) {
	  sc = e2+w2->iii; 
	} else if (t+1 == ey) {
	  sc = e2+w1->ii;
	}
      } else {
	e3 = e2; e2 = e1;
	e1 = last[j+2].C;
	sc = max(e2+w2->iii, max(e1+w1->ii,e3+w3->iv)) ;
      }
      if (sc < (cd= cur[j].D)) {
	sc = cd; 
	cur[j-3].D = cd-gext;
      } else cur[j-3].D =max(cd, sc-gop)-gext;
      if (sc < (ci= last[j].I)) {
	sc = ci;
	cur[j].I = ci - gext;
      } else cur[j].I = max(sc-gop,ci)-gext;
      cur[j].C = sc;
      w3 = w2; w2 = w1;
    }
    swap((void *)&last, (void *)&cur);
  }
  for (i = 0; i <= ey-y+1; i++) last[i].I = cur[i].I;
  if (*row1 != last) swap((void *)row1, (void *)row2);
}

static void
init_row2(int *row, int ld) {
  int i;
  for (i = 0; i < ld; i++) row[i] = 0;
}

static void init_ROW(st_ptr row, int ld) {
  int i;
  for (i = 0; i < ld; i++) row[i].I = row[i].D = row[i].C = 0;
}

static match_ptr
combine(match_ptr x1, match_ptr x2, int st) {
  match_ptr x;

  if (x1 == NULL) return x2;
  for (x = x1; x->next; x = x->next);
  x->next = x2;
  if (st) {
    for (x = x2; x; x = x->next) {
      x->j++;
      if (x->l == 3 || x->l == 4) break;
    }
    x->l--;
  }
  return x1;
}

/*
   global use the two upwards and downwards score only linear
   space global alignment subroutine to recursively build the
   alignment.
*/

match_ptr
global(int x, int y, int ex, int ey,
       int **wgts, int gop, int gext,
       unsigned char *dnap, unsigned char *pro, int N1, int N2,
       struct f_struct *f_str)
{
  int m;
  int m1, m2;
  match_ptr x1, x2, mm1, mm2;

  /*printf("%d %d %d %d %d %d\n", x,y, ex, ey, N1, N2);*/
  /*
    if the space required is limited, we can do a quadratic space
    algorithm to find the alignment.
  */

  if (ex <= x) {
    mm1  = NULL;
    for (m = y+3; m <= ey; m+=3) {
      x1 = (match_ptr) ckalloc(sizeof(match_node));
      x1->l = 5; x1->next = mm1; 
      if (mm1== NULL) mm2 = x1;
      mm1 = x1;
    }
    if (ex == x) {
      if ((ey-y) % 3 != 0) {
	x1  = (match_ptr) ckalloc(sizeof(match_node));
	x1->l = ((ey-y) % 3) +1; x1->next = NULL;
	if (mm1) mm2->next = x1; else mm1 = x1;
      } else mm2->l = 4;
    }
    return mm1;
  }
  if (ey <= y) {
    mm1  = NULL;
    for (m = x; m <= ex; m++) {
      x1 = (match_ptr) ckalloc(sizeof(match_node));
      x1->l = 0; x1->next = mm1; mm1 = x1;
    }
    return mm1;
  }
  if (ex -x < SGW1 && ey-y < SGW2) 
    return small_global(x,y,ex,ey,wgts, gop, gext,  dnap, pro, N1, N2,f_str);
  m = (x+ex)/2;
  /*     
	 Do the score only global alignment from row x to row m, m is
	 the middle row of x and ex. Store the information of row m in
	 upC, upD, and upI.
  */
  global_up(&f_str->up, &f_str->tp,  x, y, m, ey,
	    wgts, gop, gext,
	    dnap, pro, N1, f_str);
  /* 
     Do the score only global alignment downwards from row ex
     to row m+1, store information of row m+1 in downC downI and downD
  */
  global_down(&f_str->down, &f_str->tp, m+1, y, ex, ey,
	      wgts, gop, gext,
	      dnap, pro, N2, f_str);

  /*
    Use this information for row m and m+1 to find the crossing
    point of the best alignment with the middle row. The crossing
    point is given by m1 and m2. Then we recursively call global
    itself to compute alignments in two smaller regions found by
    the crossing point and combine the two alignments to form a
    whole alignment. Return that alignment.
  */
  if (find_best(f_str->up, f_str->down, &m1, &m2, ey-y+1, y, gop)) {
    x1 = global(x, y, m, m1, wgts, gop, gext, dnap, pro, N1, 0, f_str);
    x2 = global(m+1, m2, ex, ey, wgts, gop, gext, dnap, pro, 0, N2, f_str);
    if (m1 == m2) x1 = combine(x1,x2,1);
    else x1 = combine(x1, x2,0);
  } else {
    x1 = global(x, y, m-1, m1, wgts, gop, gext, dnap, pro, N1, 1, f_str);
    x2 = global(m+2, m2, ex, ey, wgts, gop, gext, dnap, pro, 1, N2, f_str);
    mm1 = (match_ptr) ckalloc(sizeof(match_node));
    mm1->i = m; mm1->l = 0; mm1->j = m1;
    mm2 = (match_ptr) ckalloc(sizeof(match_node));
    mm2->i = m+1; mm2->l = 0; mm2->j = m1;
    mm1->next = mm2; mm2->next = x2;
    x1 = combine(x1, mm1, 0);
  }
  return x1;
}

static int
find_best(st_ptr up, st_ptr down, int *m1, int *m2, int ld, int y, int gop) {

  int i, best = -1000, j = 0, s1, s2, s3, s4, st;

  for (i = 1; i < ld; i++) {
    s2 = up[i].C + down[i].C;
    s4 = up[i].I + down[i].I + gop;
    if (best < s2) {
      best = s2; j = i; st = 1;
    }
    if (best < s4) {
      best = s4; j = i; st = 0;
    }
  }
  *m1 = j-1+y;
  *m2 = j+y;
  /*printf("score=%d\n", best);*/
  return st;
} 

/*
   An alignment is represented as a linked list whose element
   is of type match_node. Each element represent an edge in the
   path of the alignment graph. The fields of match_node are
   l ---  gives the type of the edge.
   i, j --- give the end position.
*/

static match_ptr
small_global(int x, int y, int ex, int ey, 
	     int **wgts, int gop, int gext,
	     unsigned char *dnap, unsigned char *pro,
	     int N1, int N2, struct f_struct *f_str) {

  /* int C[SGW1+1][SGW2+1], st[SGW1+1][SGW2+1], D[SGW2+7], I[SGW2+1]; */
  int i, j, e, sc, score, del, k, t,  ci, cd;
  int *cI, *cD, *cC, *lC, *cst, e2, e3, e4;
  match_ptr mp, first;
  struct wgt *wt, *ww;

  /*printf("small_global %d %d %d %d\n", x, y, ex, ey);*/
  sc = -gop-gext; f_str->smgl_s.C[0][0] = 0; 

  cI = f_str->smgl_s.I;
  if (N1) cI[0] = -gext; else cI[0] = sc;

  for (j = 1; j <= ey-y+1; j++) {
    if (j % 3== 0) {
      f_str->smgl_s.C[0][j] = sc;
      sc -= gext;
      cI[j] = sc-gop;
    }
    else {
      cI[j] = f_str->smgl_s.C[0][j] = -10000;
    }
    f_str->smgl_s.st[0][j] = 5;
  }

  lC = &f_str->smgl_s.C[0][0];
  cD = f_str->smgl_s.D; cD[0] = cD[1] = cD[2] = -10000;
  for (i = 1; i <= ex-x+1; i++) {
    cC = &f_str->smgl_s.C[i][0];	
    wt = f_str->weight1[pro[i+x-1]]; cst = &f_str->smgl_s.st[i][0];
    for (j = 0; j <=ey-y+1; j++) {
      ci = cI[j];
      cd= cD[j];
      t = j+y;
      ww = &wt[(unsigned char) dnap[t-3]];
      if (j >= 4) {
	sc = lC[j-3]+ww->iii; e2 = lC[j-2]+ww->ii;  
	e4 = lC[j-4]+ww->iv; del = 3;
	if (e2 > sc) { sc = e2; del = 2;}
	if (e4 >= sc) { sc = e4; del = 4;}
      } else {
	if (j == 3) {
	  sc = lC[0]+ww->iii; del =3;
	} else if (j == 2) {
	  sc = lC[0]+ww->ii; del = 2;
	} else {sc = -10000; del = 0;}
      }
      if (sc < ci) {
	sc = ci; del = 0; 
      }
      if (sc <= cd) {
	sc = cd;
	del = 5;
      }
      cC[j] = sc;
      sc -= gop;
      if (sc <= cd) {
	del += 10;
	cD[j+3] = cd - gext;
      } else cD[j+3] = sc -gext;
      if (sc < ci) {
	del += 20;
	cI[j] = ci-gext;
      } else cI[j] = sc-gext;
      *(cst++) = del;
    }
    lC = cC;
  }
  /*printf("small global score =%d\n", f_str->smgl_s.C[ex-x+1][ey-y+1]);*/
  if (N2 && cC[ey-y+1] <  ci+gop) f_str->smgl_s.st[ex-x+1][ey-y+1] =0;
  first = NULL; e = 1;
  for (i = ex+1, j = ey+1; i > x || j > y; i--) {
    mp = (match_ptr) ckalloc(sizeof(match_node));
    mp->i = i-1;
    k  = (t=f_str->smgl_s.st[i-x][j-y])%10;
    mp->j = j-1;
    if (e == 5 && (t/10)%2 == 1) k = 5;
    if (e == 0 && (t/20)== 1) k = 0;
    if (k == 5) { j -= 3; i++; e=5;}
    else {j -= k;if (k==0) e= 0; else e = 1;}
    mp->l = k;
    mp->next = first;
    first = mp;
  }

  /*	for (i = 0; i <= ex-x; i++) {
	for (j = 0; j <= ey-y; j++) 
	printf("%d ", C[i][j]);
	printf("\n");
	}
  */
  return first;	
}

#define XTERNAL
#include "upam.h"

void
display_alig(int *a, unsigned char *dna, unsigned char *pro,
	     int length, int ld, struct f_struct *f_str)
{
  int len = 0, i, j, x, y, lines, k, iaa;
  static char line1[100], line2[100], line3[100],
    tmp[10] = "         ", *st;
  char *dna1, c1, c2, c3;

  line1[0] = line2[0] = line3[0] = '\0'; x= a[0]; y = a[1]-3;

  printf("\n%5d\n%5d", y+3, x);
  for (len = 0, j = 2, lines = 0; j < length; j++) {
    i = a[j];
    line3[len] = ' ';
    switch (i) {
    case 3: 
      y += 3;
      line2[len] = NCBIstdaa[iaa=pro[x++]];
      line1[len] = f_str->weight_c[iaa][(unsigned char) dna[y]].c5;
      if (line1[len] != f_str->weight_c[iaa][(unsigned char) dna[y]].c3)
	line3[len] = f_str->weight_c[iaa][(unsigned char) dna[y]].c3;
      break;
    case 2:
      y += 2;
      line1[len] = '\\';
      line2[len++] = ' ';
      line2[len] = NCBIstdaa[iaa=pro[x++]];
      line1[len] = f_str->weight_c[iaa][(unsigned char) dna[y]].c2;
      line3[len] = f_str->weight_c[iaa][(unsigned char) dna[y]].c3;
      break;
    case 4:
      y += 4;
      line1[len] = '/';
      line2[len++] = ' ';
      line2[len] = NCBIstdaa[iaa=pro[x++]];
      line1[len] = f_str->weight_c[iaa][(unsigned char) dna[y]].c4;
      line3[len] = f_str->weight_c[iaa][(unsigned char) dna[y]].c3;
      break;
    case 5:
      y += 3;
      line1[len] = f_str->weight_c[0][(unsigned char) dna[y]].c3;
      line2[len] = '-';
      break;
    case 0:
      line1[len] = '-';
      line2[len] = NCBIstdaa[pro[x++]];
      break;
    }
    len++;
    line1[len] = line2[len]  = line3[len]  = '\0'; 
    if (len >= WIDTH) {
      for (k = 10; k <= WIDTH; k+=10) 
	printf("    .    :");
      if (k-5 < WIDTH) printf("    .");
      c1 = line1[WIDTH]; c2 = line2[WIDTH]; c3 = line3[WIDTH];
      line1[WIDTH] = line2[WIDTH] = line3[WIDTH] = '\0';
      printf("\n     %s\n     %s\n     %s\n", line1, line3, line2);
      line1[WIDTH] = c1; line2[WIDTH] = c2;
      strncpy(line1, &line1[WIDTH], sizeof(line1)-1);
      strncpy(line2, &line2[WIDTH], sizeof(line2)-1);
      strncpy(line3, &line3[WIDTH], sizeof(line3)-1);
      len = len - WIDTH;
      printf("\n%5d\n%5d", y+3, x);
    }
  }
  for (k = 10; k < len; k+=10) 
    printf("    .    :");
  if (k-5 < len) printf("    .");
  printf("\n     %s\n     %s\n     %s\n", line1, line3, line2);
}


/* alignment store the operation that align the protein and dna sequence.
   The code of the number in the array is as follows:
   0:     delete of an amino acid.
   2:     frame shift, 2 nucleotides match with an amino acid
   3:     match an  amino acid with a codon
   4:     the other type of frame shift
   5:     delete of a codon
   

   Also the first two element of the array stores the starting point 
   in the protein and dna sequences in the local alignment.

   Display looks like where WIDTH is assumed to be divisible by 10.

    0    .    :    .    :    .    :    .    :    .    :    .    :
     AACE/N\PLK\G\HK\Y/LWA\S\C\E/P\PRIRZ/G\HK\Y/LWA\S\C\E/P\PRIRZ
          I S   G S  V F   N R Q L A     G S  V F   N R Q L A    
     AACE P P-- G HK Y TWA A C E P P---- G HK Y TWA A C E P P----

   60    .    :    .    :    .    :    .    :    .    :    .    :
     /G\HK\Y/LWA\S\C\E/P\PRIRZ/G\HK\Y/LWA\S\C\E/P\PRIRZ/G\HK\Y/LW
      G S  V F   N R Q L A     G S  V F   N R Q L A     G S  V F 
      G HK Y TWA A C E P P---- G HK Y TWA A C E P P---- G HK Y TW

For frame shift, the middle row show the letter in the original sequence,
and the letter in the top row is the amino acid that is chose by the 
alignment (translated codon chosen from 4 nucleotides, or 2+1).
*/

/* fatal - print message and die */
void
fatal(msg)
     char *msg;
{
  fprintf(stderr, "%s\n", msg);
  exit(1);
}

/* 10-Feb-2010 - fz_walign modified to ensure that the final alignment
   overlaps the initial lz_band() region.  In earlier versions, the
   final alignment (using pam2p[0]) might have been outside the band
   region */

void
fz_walign (const unsigned char *aa0, int n0,
	   const unsigned char *aa1, int n1,
	   int frame, int max_res,
	   struct pstruct *ppst, 
	   struct f_struct *f_str, 
	   struct a_res_str *a_res,
	   int score_thresh)
{
  int score;
  int i, last_n1, itemp, n10;
  int hoff, nt_min, nt_max, n_nt, n_aa, w_fact;
  int l_min, l_max, window;
  unsigned char *fs, *fd;
  /*
  unsigned char *aa1_min_s, aa1_max_s;
  */
  unsigned char *local_aa1;
  int optflag_s;
  int itx;
  unsigned char *aa1x;
  struct score_count_s s_info = {0,0,0,0};;
  
#ifndef TFAST
  do_fastz (f_str->aa0x, n0, aa1, n1, f_str->aa0v, ppst, f_str, &a_res->rst, &hoff, 1, &s_info);
#else
  /* make translated sequence */
  last_n1 = 0;
  aa1x = f_str->aa1x;
  for (itx= frame*3; itx< frame*3+3; itx++) {
    n10  = saatran(aa1,&aa1x[last_n1],n1,itx);
    last_n1 += n10+1;
  }
  n10 = last_n1-1;

  /* do_fastz (lz_band) also needs a pre-computed number series */
  if (frame == 0) {
    pre_com(aa1, n1, f_str->aa1v);
  }
  else {
    pre_com_r(aa1, n1, f_str->aa1v);
  }

  do_fastz (aa0, n0, f_str->aa1x, n10, f_str->aa1v, ppst, f_str, &a_res->rst, &hoff, 1, &s_info);
#endif

  if (a_res->rst.score[ppst->score_ix] <= score_thresh) {
    a_res->sw_score = 0;
    a_res->n1 = n1;
    return;
  }

#ifndef TFAST
  window = min(n1, ppst->param_u.fa.optwid);
  l_min = max(0, -window - hoff);
  l_max = min(n1, n0-hoff+window);
  local_aa1 = (unsigned char *)aa1;
  if (l_min > 0 || l_max < n1 -1 ) {
    local_aa1 = (unsigned char *)calloc(l_max - l_min + 2,sizeof(unsigned char));
    local_aa1++;
    memcpy(local_aa1,aa1+l_min, l_max - l_min);
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
#else
  window = min(n0, ppst->param_u.fa.optwid);
  if (frame==0) {
    l_min = max(0,(hoff-window)*3);
    l_max = min((hoff+window+n0)*3,n1);

    local_aa1 = (unsigned char *)aa1;
    if (l_min > 0 || l_max < n1 -1 ) {
      local_aa1 = (unsigned char *)calloc(l_max - l_min + 2,sizeof(unsigned char));
      local_aa1++;
      memcpy(local_aa1,aa1+l_min, l_max - l_min);
    }

    /*
    if (l_min > 0) {
      aa1_min_s = aa1[l_min-1];
      aa1[l_min-1] = '\0';
    }
    if (l_max < n1-1) {
      aa1_max_s = aa1[l_max];
      aa1[l_max] = '\0';
    }
    */
    /* re-do precomputed codon number series for limited region */
    pre_com(local_aa1, l_max - l_min, f_str->aa1v);
  }
  else {
    /* things are more complicated here because the mapping of hoff is
       with respect to the reversed aa1 */

    l_max = n1 - max(0,(hoff-window)*3);
    l_min = n1 - min((hoff+window+n0)*3,n1);

    local_aa1 = (unsigned char *)aa1;
    if (l_min > 0 || l_max < n1 -1 ) {
      local_aa1 = (unsigned char *)calloc(l_max - l_min + 2,sizeof(unsigned char));
      local_aa1++;
      memcpy(local_aa1,aa1+l_min, l_max - l_min);
    }

    /*
    if (l_min > 0) {
      aa1_min_s = aa1[l_min-1];
      aa1[l_min-1] = '\0';
    }
    if (l_max < n1-1) {
      aa1_max_s = aa1[l_max];
      aa1[l_max] = '\0';
    }
    */

    pre_com_r(local_aa1, l_max - l_min, f_str->aa1v);
  }
#endif

  a_res->sw_score = 
    pro_dna(
#ifndef TFAST
	    aa1+l_min, l_max - l_min,
	    f_str->aa0v, n0-2,
#else
	    aa0, n0,
	    f_str->aa1v, l_max - l_min-2,
#endif
	    ppst->pam2[0],
	    -ppst->gdelval, -ppst->ggapval, -ppst->gshift,
	    f_str, f_str->max_res, a_res);

  if (l_min > 0 || l_max < n1 - 1) { free(--local_aa1); }
  /*
  if (l_min > 0) {
    aa1[l_min-1] = aa1_min_s;
  }
  if (l_max < n1 - 1) {
    aa1[l_max] = aa1_max_s;
  }
  */
#ifndef TFAST
  a_res->min0 += l_min;
  a_res->max0 += l_min;
#else
  if (frame==1) {
    a_res->min1 += n1 - l_max;
    a_res->max1 += n1 - l_max;
  }
  else {
    a_res->min1 += l_min;
    a_res->max1 += l_min;
  }    
#endif

  /* display_alig(f_str->res,f_str->aa0v+2,aa1,*nres,n0-2,f_str); */
}

/*
  fz_malign is a recursive interface to fz_walign() that is called
  from do_walign(). fz_malign() first does an alignment, then checks
  to see if the score is greater than the threshold. If so, it tries
  doing a left and right alignment.

  In this implementation, the DNA sequence is preserved as DNA for
  TFAST, so that it can be sub-setted and translated correctly.  Thus,
  the translation required for f_str->aa1x and f_str->aa1v is done at
  each recursive level (in fz_walign).
 */

struct a_res_str *
fz_malign (const unsigned char *aa0, int n0,
	   const unsigned char *aa1, int n1,
	   int frame, 
	   int score_thresh, int max_res,
	   struct pstruct *ppst,
	   struct f_struct *f_str,
	   struct a_res_str *cur_ares,
	   int first_align)
{
  struct a_res_str *tmpl_ares, *tmpr_ares, *this_ares;
  struct a_res_str *mtmpl_ares, *mtmpr_ares, *mt_next;
  int sq_start, sq_end, sq_save;
  int hoff, score_ix;
  int min_alen;
  struct rstruct rst;
  unsigned char *local_aa1;
  /*   char save_res; */
  int iphase, i;
  unsigned char *fd;
  int max_sub_score = -1;

  score_ix = ppst->score_ix;

#ifdef TFAST  
  min_alen = min(n0,MIN_LOCAL_LEN)*3;	/* n0 in aa, min_alen in nt */
#else 
  min_alen = min(n0/3,MIN_LOCAL_LEN);	/* no in nt, min_alen in aa */
#endif

  /* now we need alignment storage - get it */
  if ((cur_ares->res = (int *)calloc((size_t)max_res,sizeof(int)))==NULL) {
    fprintf(stderr,"*** ERROR [%s:%d] - cannot allocate alignment results array %d\n",
	    __FILE__, __LINE__, max_res);
    exit(1);
  }

  cur_ares->next = NULL;

  fz_walign(aa0, n0, aa1, n1, frame, max_res, ppst, f_str, cur_ares, (first_align ? 1 : score_thresh));

  /* in cur_ares, min0,max0 are always protein, min1,max1 are always
     DNA, but n0 could be protein or DNA, depending on
     FASTY/TFASTY */

  if (!ppst->do_rep || cur_ares->rst.score[ppst->score_ix] <= score_thresh) { 
    return cur_ares;
  }

  /* have a score >= threshold - try left and right */

  /* in code below, cur_ares->min0/max0 always refers to aa
     cur_ares->min1/max1 always refers to nt

     however, things are more complex because if frame==1, then
     offsets are from the end (n1), not the beginning.  There is no
     frame==1 for fasty, only for TFASTY
  */
  cur_ares->v_start = sq_start = 0;
#ifdef TFAST
  if (frame == 0) {sq_end = cur_ares->min1-1;}   /* aa1[sq_start --> sq_end] */
  else {sq_end = n1 - cur_ares->max1;}
  sq_save = sq_end;
#else
  sq_save = sq_end = cur_ares->min0;
#endif
  cur_ares->v_len = sq_end - sq_start;

  if (cur_ares->v_len >= min_alen) { /* try the left  */
      /* allocate a_res */	
      tmpl_ares = (struct a_res_str *)calloc(1, sizeof(struct a_res_str));
      local_aa1 = (unsigned char *)calloc(cur_ares->v_len+2,sizeof(unsigned char));
      local_aa1++;
      memcpy(local_aa1, aa1, cur_ares->v_len);

      /* 
      save_res = aa1[sq_save];
      aa1[sq_save] = '\0';
      */
      tmpl_ares = fz_malign(aa0, n0, local_aa1, cur_ares->v_len,
			    frame, score_thresh, max_res, 
			    ppst, f_str, tmpl_ares,0);

      free(--local_aa1);
      /*       aa1[sq_save] = save_res; */

      if (tmpl_ares->rst.score[ppst->score_ix] > score_thresh) {
	max_sub_score = tmpl_ares->rst.score[ppst->score_ix];
#ifdef TFAST
	if (frame == 1) {
	  for (this_ares = tmpl_ares; this_ares; this_ares = this_ares->next) {
	    this_ares->v_start += n1 - sq_end;
	    this_ares->min1 += n1 - sq_end;
	    this_ares->max1 += n1 - sq_end;
	  }
	}
#endif
      }
      else {
	if (tmpl_ares->res) free(tmpl_ares->res);
	free(tmpl_ares);
	tmpl_ares=NULL;
      }
    }
    else {tmpl_ares = NULL;}

  /* now the right end */
  /* for fasty -- max positions refer to the aa,codon, not the next
     residue, so they must be incremented */

    sq_end = n1;
#if TFAST
    if (frame == 0) {sq_start = cur_ares->max1+1;}
    else {sq_start = n1 - cur_ares->min1;}
#else
    sq_start = cur_ares->max0+1;
#endif
    sq_save = sq_start-1;
    cur_ares->v_len = sq_end - sq_start;

  if (cur_ares->v_len >= min_alen) { /* try the right  */
      /* allocate a_res */	
      tmpr_ares = (struct a_res_str *)calloc(1, sizeof(struct a_res_str));

      /* find boundaries */
      local_aa1 = (unsigned char *)calloc(cur_ares->v_len+2,sizeof(unsigned char));
      local_aa1++;
      memcpy(local_aa1,aa1+sq_start,cur_ares->v_len);
      /*
      save_res = aa1[sq_save];
      aa1[sq_save] = '\0';
      */

      tmpr_ares = fz_malign(aa0, n0, 
			    local_aa1, cur_ares->v_len,
			    frame,
			    score_thresh, max_res,
			    ppst, f_str, tmpr_ares,0);
      free(--local_aa1);
      /* 
      aa1[sq_save] = save_res;
      */

      if (tmpr_ares->rst.score[ppst->score_ix] > score_thresh) {
	/* adjust the left boundary */
	for (this_ares = tmpr_ares; this_ares; this_ares = this_ares->next) {
#ifndef TFAST
	  this_ares->min0 += sq_start;
	  this_ares->max0 += sq_start;
#else
	  if (frame == 0) {
	    this_ares->v_start += sq_start;
	    this_ares->min1 += sq_start;
	    this_ares->max1 += sq_start;
	  }
#endif
	}
	if (tmpr_ares->rst.score[ppst->score_ix] > max_sub_score) {
	  max_sub_score = tmpr_ares->rst.score[ppst->score_ix];
	}
      }
      else {
	if (tmpr_ares->res) free(tmpr_ares->res);
	free(tmpr_ares);
	tmpr_ares=NULL;
      }
    }
    else {tmpr_ares = NULL;}

    if (max_sub_score <= score_thresh) return cur_ares;

    cur_ares = merge_ares_chains(cur_ares, tmpl_ares, score_ix, "left");
    cur_ares = merge_ares_chains(cur_ares, tmpr_ares, score_ix, "right");

    return cur_ares;
}

/* do_walign() can be called with aa0,n0 as nt (FASTY) or 
   aa0,n0 as aa (TFASTY).  if aa0 is nt, then f_str->aa0x,y have the
   translations already.  if aa0 is aa, then f_str->aa1x,y must be
   generated.
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
  int hoff, use_E_thresholds_s, optflag_s, optcut_s, optwid_s, score;
  int last_n1, itx, itt, n10, iphase;
  unsigned char *fs, *fd;
  struct rstruct rst;
#ifdef DEBUG
  unsigned long adler32_crc;
#endif

  *have_ares = 0x3;	/* set 0x2 bit to indicate local copy */

  if ((a_res = (struct a_res_str *)calloc(1, sizeof(struct a_res_str)))==NULL) {
    fprintf(stderr,"*** ERROR [%s:%d] - cannot allocate a_res [%lu]",
	    __FILE__, __LINE__, sizeof(struct a_res_str));
    return NULL;
  }

#ifdef DEBUG
  adler32_crc = adler32(1L,aa1,n1);
#endif

  f_str->frame = frame;	/* need frame for later pre_cons() in calcons() */

  use_E_thresholds_s = ppst->param_u.fa.use_E_thresholds;
  optflag_s = ppst->param_u.fa.optflag;
  optcut_s = ppst->param_u.fa.optcut;
  optwid_s = ppst->param_u.fa.optwid;
  ppst->param_u.fa.use_E_thresholds = 0;
  ppst->param_u.fa.optflag = 1;
  if (!ppst->param_u.fa.optwid_set) {
    ppst->param_u.fa.optwid *= 2;
  }

  a_res = fz_malign(aa0, n0, aa1, n1, frame,
		    repeat_thresh, f_str->max_res,
		    ppst, f_str, a_res, 1);

#ifdef DEBUG
  if (adler32(1L,aa1,n1) != adler32_crc) {
    fprintf(stderr,"*** ERROR [%s:%d] adler32_crc mismatch n1: %d\n",__FILE__, __LINE__, n1);
  }
#endif

  a_res_index = 0;
  for (tmp_a_res=a_res; tmp_a_res; tmp_a_res = tmp_a_res->next) {
    tmp_a_res->index = a_res_index++;
  }

  ppst->param_u.fa.use_E_thresholds = use_E_thresholds_s;
  ppst->param_u.fa.optflag = optflag_s;
  ppst->param_u.fa.optwid = optwid_s;
  return a_res;
}

void
pre_cons(const unsigned char *aa1, int n1, int frame, struct f_struct *f_str) {

#ifdef TFAST
  int i, last_n1, itemp, n10;
  unsigned char *fs, *fd;
  int itx;

  /* make a precomputed codon number series */
  if (frame==0) {
    pre_com(aa1, n1, f_str->aa1v);
  }
  else { /* must do things backwards */
    pre_com_r(aa1, n1, f_str->aa1v);
  }
#endif
}

/* aln_func_vals - set up aln.qlfact, qlrev, llfact, llmult, frame, llrev */
/* call from calcons, calc_id, calc_code */
void 
aln_func_vals(int frame, struct a_struct *aln) {

#ifndef TFAST
  aln->llrev = 0;
  aln->llfact = 1;
  aln->llmult = 1;
  aln->qlfact = 3;
  aln->frame = frame;
  if (frame > 0) aln->qlrev = 1;
  else aln->qlrev = 0;
  aln->llrev = 0;
#else	/* TFASTX */
  aln->qlfact = 1;
  aln->qlrev = 0;
  aln->llfact = 3;
  aln->llmult = 1;
  aln->frame = frame;
  if (frame > 0) aln->llrev = 1;
  else aln->llrev = 0;
  aln->qlrev=0;
#endif	/* TFASTX */
}

#include "structs.h"
#include "a_mark.h"

extern int align_type(int score, char sp0, char sp1, int nt_align, struct a_struct *aln, int pam_x_id_sim);

extern void
process_annot_match(int *itmp, int *pam2aa0v, 
		    long ip, long ia, char *sp1, char *sp1a, const unsigned char *sq,
		    struct annot_entry *annot_arr_p, int n_annots, char **ann_comment,
		    void *annot_stack, int *have_push_features, int *v_delta,
		    int *d_score_p, int *d_ident_p, int *d_alen_p, int *d_gaplen_p,
		    struct domfeat_data **left_domain_head_p,
		    struct domfeat_data *left_domain_p,
		    long *left_end_p, int init_score);

extern int
next_annot_match(int *itmp, int *pam2aa0v, 
		 long ip, long ia, char *sp1, char *sp1a, const unsigned char *sq,
		 int i_annot, int n_annot, struct annot_entry **annot_arr, char **ann_comment,
		 void *annot_stack, int *have_push_features, int *v_delta,
		 int *d_score_p, int *d_ident_p, int *d_alen_p, int *d_gaplen_p,
		 struct domfeat_data **left_domain_head_p,
		 struct domfeat_data *left_domain_p,
		 long *left_domain_end, int init_score);

extern void
close_annot_match (int ia, void *annot_stack, int *have_push_features,
		   int *d_score_p, int *d_ident_p, int *d_alen_p, int *d_gaplen_p,
		   struct domfeat_data **left_domain_p,
		   long *left_end_p, int init_score);

extern void
comment_var(long i0, char sp0, long i1, char sp1, char o_sp1, char sim_char,
	    const char *ann_comment, struct dyn_string_str *annot_var_dyn,
	    int target, int d_type);

void
display_push_features(void *annot_stack, struct dyn_string_str *annot_var_dyn,
		      long i0_pos, char sp0, long i1_pos, char sp1, char sym, 
		      int score, double comp, int sw_score, int n0, int n1,
		      void *pstat_void, int d_type);

#define DP_FULL_FMT 1	/* Region: score: bits: id: ... */
#define Q_TARGET 0
#define L_TARGET 1

int seq_pos(int pos, int rev,int off);

/* values of calc_func_mode */
#define CALC_CONS 1
#define CALC_CODE 2
#define CALC_ID   3
#define CALC_ID_DOM   4

/* add_annot_code: adds annotation codes to struct dyn_string_str ann_code_dyn */
void
add_annot_code(int have_ann, char sp0, char sp1, 
	       char ann_aa1_i1,
	       long q_off_pos, long l_off_pos, char sim_sym_code,
	       struct dyn_string_str *ann_code_dyn)
{
  char ann_ch0, ann_ch1;
  char tmp_astr[MAX_STR];

  ann_ch0 = ann_ch1 = '\0';

  if (have_ann && ann_aa1_i1 != ' ') {
    ann_ch0 = 'X';
    ann_ch1 = ann_aa1_i1;
  }
  else {return;}

  if (!(ann_ch1 == '['  || ann_ch1 == ']')) {
    sprintf(tmp_astr, "|%c%c:%ld%c%c%ld%c",
	    ann_ch0,ann_ch1, q_off_pos+1,sp0,
	    sim_sym_code, l_off_pos+1,sp1);
    dyn_strcat(ann_code_dyn, tmp_astr);
  }
}

/* universal alignment code builder for calc_cons_a(), calc_code(), and calc_id() */
/* see cal_cons2.c/calc_cons_u() for strategy */

int
calc_cons_u( /* inputs */
	    const unsigned char *aa0, int n0,
	    const unsigned char *aa1, int n1,
	    struct a_res_str *a_res,	/* alignment encoding */
	    struct pstruct *ppst,
	    struct f_struct *f_str,
	    void *pstat_void,
	    /* annotation stuff */
	    const unsigned char *ann_arr,
	    const unsigned char *aa0a, const struct annot_str *annot0_p,
	    const unsigned char *aa1a, const struct annot_str *annot1_p,
	    int calc_func_mode, 	/* CALC_CONS, CALC_CODE, CALC_ID */
	    int display_code, 	/* used only by CALC_CODE */
	    /* outputs */
	    int *nc,
	    char *seqc0, char *seqc1, char *seqca, int *cumm_seq_score,
	    char *seqc0a, char *seqc1a,
	    struct a_struct *aln,
	    int *score_delta, 
	    struct dyn_string_str *annot_var_dyn,
	    struct dyn_string_str *align_code_dyn)
{
  int i0, i1, i, j;
  int lenc, not_c, itmp, ngap_p, ngap_d, nfs;
  int *i_spa;
  char *sp0_p, *sp0a_p, *sp1_p, *sp1a_p, *spa_p, t_spa;
  char sp0_c, sp1_c, spa_c;	/* used for CALC_ID, CALC_CODE */
  char sp0a_c, sp1a_c;		/* used for CALC_CODE */

  struct update_code_str *update_data_p;

  const unsigned char *sq;
  unsigned char aap;

  const unsigned char *ap0, *ap1;
  const unsigned char *ap1a;	/* ap1 always points to protein, and
				   only protein has annotations */
  const struct annot_str *annotp_p;	/* protein annotations from annot_str */
  int comment_target;

  int *rp, *rpmax;
  int have_ann;

  /* variables for variant changes/region scores */
  char tmp_str[MAX_LSTR];
  void *annot_stack;
  int have_push_features, prev_match, *have_push_features_p;

  char *sim_sym = aln_map_sym[MX_ACC];
  struct annot_entry **s_annotp_arr_p;
  int  i1_annot, v_delta, v_tmp;
  long i0_offset, i1_offset;

  long i1_left_end;
  int show_code, annot_fmt, start_flag;

  int d1_score, d1_ident, d1_alen, d1_gaplen;
  struct domfeat_data *left_domain_list1, *left_domain_head1;

  char *ann_comment;

  *score_delta = 0;
  d1_score = d1_ident = d1_alen = d1_gaplen = 0;
  i1_left_end = -1;
  left_domain_head1 = left_domain_list1 = NULL;

  NULL_dyn_string(annot_var_dyn);

  if (ppst->ext_sq_set) {sq = ppst->sqx;}
  else {sq = ppst->sq;}

#ifndef TFAST	/* FASTX */
  comment_target = 1;
  aln->amin1 = aln->smin1 = a_res->min0;	/* prot */
  aln->amin0 = aln->smin0 = a_res->min1;	/* DNA */

  i0_offset = aln->q_offset;
  i1_offset = aln->l_offset;

  ap0 = f_str->aa0v;	/* translated DNA */
  ap1 = aa1;		/* protein */

  ap1a = aa1a;
  annotp_p = annot1_p;

  if (calc_func_mode == CALC_CONS) {
    have_ann = (seqc0a !=NULL && aa1a != NULL);
    sp0_p = seqc0;	/* translated DNA */
    sp1_p = seqc1;	/* protein */
    spa_p = seqca;
    sp1a_p = seqc1a;	/* protein library can have annotation */
    sp0a_p = seqc0a;	/* sp0a is always ' ' - no translated annotation */
    annot_fmt = DP_FULL_FMT;
  }
  else if (calc_func_mode == CALC_ID || calc_func_mode == CALC_ID_DOM) {
    /* does not require aa0a/aa1a, only for variants */
    have_ann = ((annotp_p && annotp_p->n_annot > 0) || (annot0_p && annot0_p->n_annot > 0));
    sp0_p = &sp0_c;
    sp1_p = &sp1_c;
    spa_p = &spa_c;
    sp0a_p = &sp0a_c;
    sp1a_p = &sp1a_c;
    annot_fmt = 3;
  }
  else if (calc_func_mode == CALC_CODE) {
    spa_p = &spa_c;
    sp0_p = &sp0_c;
    sp1_p = &sp1_c;
    sp0a_p = &sp0a_c;
    sp1a_p = &sp1a_c;

    show_code = (display_code & (SHOW_CODE_MASK+SHOW_CODE_EXT));	/* see defs.h; SHOW_CODE_ALIGN=2,_CIGAR=3,_CIGAR_EXT=4 */
    annot_fmt = 2;
    if (display_code & SHOW_ANNOT_FULL) {
      annot_fmt = 1;
    }
    /* have_ann encodes number of sequences annotated */
    have_ann = 0;
    if ((annotp_p && annotp_p->n_annot > 0) || (aa1a != NULL)) { have_ann |= 2;}
    update_data_p = init_update_data(show_code);
  }
  else {
    fprintf(stderr,"*** ERROR [%s:%d] --- calc_cons_u() invalid calc_func_mode: %d\n",
	    __FILE__, __LINE__, calc_func_mode);
    exit(1);
  }

#else		/* TFASTX */
  comment_target = 0;
  aln->amin0 = aln->smin0 = a_res->min0;	/* DNA */
  aln->amin1 = aln->smin1 = a_res->min1;	/* prot */

  i1_offset = aln->q_offset;
  i0_offset = aln->l_offset;

  ap1 = aa0;		/* aa0 is protein */
  /* with fx_malign(), there is no guarantee that we have a valid f_str->aa1y, so make one */
  pre_cons(aa1,n1,aln->frame, f_str);
  ap0 = f_str->aa1v;	/* aa1 is DNA */
  ap1a = aa0a;
  annotp_p = annot0_p;

  have_ann = (seqc0a !=NULL && aa0a != NULL);
  if (calc_func_mode == CALC_CONS) {
    sp1_p = seqc0;		/* sp1 points to protein query */
    sp0_p = seqc1;		/* sp0 points to DNA */
    spa_p = seqca;
    sp1a_p = seqc0a;	/* protein query can have annotation */
    sp0a_p = seqc1a;	/* sp0a is always ' ' - no translated annotation */
    annot_fmt = DP_FULL_FMT;
  }
  else if (calc_func_mode == CALC_ID || calc_func_mode == CALC_ID_DOM)  {
    have_ann = (annotp_p && annotp_p->n_annot > 0);
    spa_p = &spa_c;
    sp0_p = &sp1_c;
    sp1_p = &sp0_c;

    sp0a_p = &sp1a_c;
    sp1a_p = &sp0a_c;
    annot_fmt = 3;

    /* does not require aa0a/aa1a, only for variants */
  }
  else if (calc_func_mode == CALC_CODE) {
    spa_p = &spa_c;
    sp0_p = &sp1_c;
    sp1_p = &sp0_c;

    sp0a_p = &sp1a_c;
    sp1a_p = &sp0a_c;

    show_code = (display_code & (SHOW_CODE_MASK+SHOW_CODE_EXT));	/* see defs.h; SHOW_CODE_ALIGN=2,_CIGAR=3,_CIGAR_EXT=4 */
    annot_fmt = 2;
    if (display_code & SHOW_ANNOT_FULL) {
      annot_fmt = 1;
    }

    /* have_ann encodes number of sequences annotated */
    if ((annotp_p && annotp_p->n_annot > 0) || (ap1a != NULL)) { have_ann |= 1;}

    update_data_p = init_update_data(show_code);
  }
  else {
    fprintf(stderr,"*** ERROR [%s:%d] --- calc_cons_u() invalid calc_func_mode: %d\n",
	    __FILE__, __LINE__, calc_func_mode);
    exit(1);
  }
#endif
  if (cumm_seq_score) i_spa = cumm_seq_score;

  rp = a_res->res;
  rpmax = &a_res->res[a_res->nres];

  lenc = not_c = aln->nident = aln->nmismatch = aln->nsim = aln->npos = ngap_p = ngap_d = nfs= 0;
  i0 = a_res->min1 - 3;
  i1 = a_res->min0;

  v_delta = 0;
  i1_annot = 0;
  annot_stack = NULL;
  s_annotp_arr_p = NULL;
  have_push_features = prev_match = 0;
  if (have_ann) {
    have_push_features_p = &have_push_features;

    if (annotp_p && annotp_p->n_annot > 0) {
      annot_stack = init_stack(64,64);
      left_domain_list1=init_domfeat_data(annotp_p);

      s_annotp_arr_p = annotp_p->s_annot_arr_p;

      while (i1_annot < annotp_p->n_annot) {
	if (s_annotp_arr_p[i1_annot]->pos >= i1+i1_offset) {break;}
	if (s_annotp_arr_p[i1_annot]->end <= i1+i1_offset) {i1_annot++; continue;}

	if (s_annotp_arr_p[i1_annot]->label == '-') {
	  process_annot_match(&itmp, NULL, i1_offset+seq_pos(i1,aln->llrev,0), i0_offset + seq_pos(i0,aln->qlrev,0),
			      sp1_p, sp1a_p, sq, s_annotp_arr_p[i1_annot], annotp_p->n_annot,  &ann_comment, 
			      annot_stack, have_push_features_p, &v_delta,
			      &d1_score, &d1_ident, &d1_alen, &d1_gaplen,
			      &left_domain_head1, &left_domain_list1[i1_annot], &i1_left_end, 0);
	}
	i1_annot++;
      }
    }
  }

  while (rp < rpmax) {
    /*    fprintf(stderr,"%d %d %d (%c) %d (%c)\n"
	  ,(int)(rp-res),*rp,i0,sq[ap0[i0]],i1,sq[ap1[i1]]);
    */

    switch (*rp++) {
    case 3:		/* match */
      i0 += 3;

      *sp1_p = sq[aap=ap1[i1]];
      *sp0_p = f_str->weight_c[aap][ap0[i0]].c5;
      itmp = ppst->pam2[0][aap][pascii[*sp0_p]];


      if (have_ann) {
	if (calc_func_mode != CALC_ID  && calc_func_mode != CALC_ID_DOM) {
	  *sp1a_p = ann_arr[ap1a[i1]];
	  *sp0a_p = ' ';
	}
	if (s_annotp_arr_p) {
	  if (i1+i1_offset == s_annotp_arr_p[i1_annot]->pos || i1+i1_offset == i1_left_end) {
	    i1_annot = next_annot_match(&itmp, ppst->pam2[0][pascii[*sp0_p]],
#ifndef TFAST
					i1_offset+seq_pos(i1,aln->llrev,0),	/* annotated target (prot) coordinate */
					i0_offset+seq_pos(i0,aln->qlrev,0),
#else
					i1_offset+seq_pos(i1,aln->qlrev,0),
					i0_offset+seq_pos(i0,aln->llrev,0),
#endif
					sp1_p, sp1a_p, sq,
					i1_annot, annotp_p->n_annot, s_annotp_arr_p,
					&ann_comment, annot_stack, have_push_features_p, &v_delta,
					&d1_score, &d1_ident, &d1_alen, &d1_gaplen,
					&left_domain_head1, left_domain_list1, &i1_left_end,
					0);

	    /* must be out of the loop to capture the last value */
	    if (ppst->sq[ap1[i1]] != *sp1_p) {
	      t_spa = align_type(itmp, *sp0_p, *sp1_p, 0, NULL, ppst->pam_x_id_sim);

	      if (calc_func_mode != CALC_ID  && calc_func_mode != CALC_ID_DOM) {
	      comment_var(
#ifndef TFAST
			  i0_offset+seq_pos(i0,aln->qlrev,0), *sp0_p,
			  i1_offset+seq_pos(i1,aln->llrev,0), *sp1_p,
#else
			  i1_offset+seq_pos(i1,aln->qlrev,0), *sp0_p,
			  i0_offset+seq_pos(i0,aln->llrev,0), *sp1_p,
#endif
			  sq[ap1[i1]], sim_sym[t_spa], ann_comment,
			  annot_var_dyn,comment_target,annot_fmt);
	      }
	      else {
		sprintf(tmp_str,"%c%d%c;",ppst->sq[ap1[i1]],i1+1,*sp1_p);
		/*  SAFE_STRNCAT(annot_var_s,tmp_str,n_annot_var_s); */
		dyn_strcat(annot_var_dyn, tmp_str);
	      }
	    }
	  }
	  prev_match = 1;
	  d1_score += itmp;
	}
	if (calc_func_mode == CALC_CONS) {sp0a_p++; sp1a_p++;}
      }

      *spa_p = align_type(itmp, *sp0_p, *sp1_p, 0, aln, ppst->pam_x_id_sim);
      d1_alen++;
      if (*spa_p == M_IDENT) {d1_ident++;}

      if (cumm_seq_score) *i_spa++ = itmp;

      if (calc_func_mode == CALC_CODE) {
#ifndef TFAST
	update_code(align_code_dyn, update_data_p, 3, *spa_p, *sp0_p, *sp1_p);
#else
	update_code(align_code_dyn, update_data_p, 3, *spa_p, *sp1_p, *sp0_p);
#endif

	if (have_ann && have_push_features) {
    	  add_annot_code(have_ann, *sp0_p, *sp1_p, *sp1a_p,
#ifndef TFAST
			 i0_offset+seq_pos(i0,aln->qlrev,0),
			 i1_offset+seq_pos(i1,aln->llrev,0),
#else
			 i1_offset+seq_pos(i1,aln->qlrev,0),
			 i0_offset+seq_pos(i0,aln->llrev,0),
#endif
			 sim_sym[*spa_p], annot_var_dyn);
	}
      }

      if (have_ann && have_push_features &&  calc_func_mode != CALC_ID) {
	display_push_features(annot_stack, annot_var_dyn,
#ifndef TFAST
			      i0_offset+seq_pos(i0,aln->qlrev,0), *sp0_p,
			      i1_offset+seq_pos(i1,aln->llrev,0), *sp1_p,
#else
			      i1_offset+seq_pos(i1,aln->qlrev,0), *sp0_p,
			      i0_offset+seq_pos(i0,aln->llrev,0), *sp1_p,
#endif
			      sim_sym[*spa_p],
			      a_res->rst.score[ppst->score_ix], a_res->rst.comp, a_res->sw_score,
			      n0, n1, pstat_void, annot_fmt);
	have_push_features = 0;
      }

      i1++;
      if (calc_func_mode == CALC_CONS) {sp0_p++; sp1_p++; spa_p++;}
      lenc++;
      break;
    case 2:		/* frame shift +2, then match */
      nfs++;
      i0 += 2;
      *sp0_p = '/';
      *sp1_p = '-';
      *spa_p = M_DEL;

      if (calc_func_mode == CALC_CODE) {
#ifndef TFAST
	update_code(align_code_dyn, update_data_p, 2, *spa_p,*sp0_p,*sp1_p);
#else
	update_code(align_code_dyn, update_data_p, 2, *spa_p,*sp1_p,*sp0_p);
#endif
      }

      if (cumm_seq_score) *i_spa++ = ppst->gshift;

      if (calc_func_mode == CALC_CONS) {
	sp0_p++; sp1_p++; spa_p++;
	if (have_ann) {*sp0a_p++ = *sp1a_p++ = ' ';}
      }

      not_c++;

      *sp1_p = sq[aap=ap1[i1]];
      *sp0_p = f_str->weight_c[aap][ap0[i0]].c2;
      itmp = ppst->pam2[0][pascii[*sp0_p]][aap];

      if (cumm_seq_score) *i_spa++ = ppst->gshift;

      if (have_ann) {
	have_push_features = 0;
	if (calc_func_mode != CALC_ID && calc_func_mode != CALC_ID_DOM) {
	  *sp1a_p = ann_arr[ap1a[i1]];
	  *sp0a_p = ' ';
	}
	if (s_annotp_arr_p) {
	  if (i1 + i1_offset == s_annotp_arr_p[i1_annot]->pos || i1+i1_offset == i1_left_end) {
	    i1_annot = next_annot_match(&itmp, ppst->pam2[0][pascii[*sp0_p]],
#ifndef TFAST
					i1_offset+seq_pos(i1,aln->llrev,0),	/* annotated target (prot) coordinate */
					i0_offset+seq_pos(i0,aln->qlrev,0),
#else
					i1_offset+seq_pos(i1,aln->qlrev,0),
					i0_offset+seq_pos(i0,aln->llrev,0),
#endif
					sp1_p, sp1a_p, sq,
					i1_annot, annotp_p->n_annot, s_annotp_arr_p,
					&ann_comment, annot_stack, have_push_features_p, &v_delta,
					&d1_score, &d1_ident, &d1_alen, &d1_gaplen,
					&left_domain_head1, left_domain_list1, &i1_left_end,0);

	    /* must be out of the loop to capture the last value */
	    if (ppst->sq[ap1[i1]] != *sp1_p) {
	      t_spa = align_type(itmp, *sp0_p, *sp1_p, 0, NULL, ppst->pam_x_id_sim);

	      if (calc_func_mode != CALC_ID && calc_func_mode != CALC_ID_DOM) {
		comment_var(
#ifndef TFAST
			    i0_offset+seq_pos(i0,aln->qlrev,0), *sp0_p,
			    i1_offset+seq_pos(i1,aln->llrev,0), *sp1_p,
#else
			    i1_offset+seq_pos(i1,aln->qlrev,0), *sp0_p,
			    i0_offset+seq_pos(i0,aln->llrev,0), *sp1_p,
#endif
			    sq[ap1[i1]], sim_sym[t_spa], ann_comment,
			    annot_var_dyn,comment_target,annot_fmt);
	      }
	      else {
		sprintf(tmp_str,"%c%d%c;",sq[ap1[i1]],i1+1,*sp1_p);
		/*   SAFE_STRNCAT(annot_var_s,tmp_str,n_annot_var_s); */
		dyn_strcat(annot_var_dyn, tmp_str);
	      }
	    }
	  }
	  d1_score += ppst->gshift;
	  d1_score += itmp;
	  prev_match = 1;
	}
	if (calc_func_mode == CALC_CONS) {sp0a_p++; sp1a_p++;}
      }

      *spa_p = align_type(itmp, *sp0_p, *sp1_p, 0, aln, ppst->pam_x_id_sim);

      if (calc_func_mode == CALC_CODE) {
#ifndef TFAST
	update_code(align_code_dyn, update_data_p, 3, *spa_p,*sp0_p,*sp1_p);
#else
	update_code(align_code_dyn, update_data_p, 3, *spa_p,*sp1_p,*sp0_p);
#endif
      }

      d1_alen++;
      if (*spa_p == M_IDENT) {d1_ident++;}

      if (have_ann && have_push_features  &&  calc_func_mode != CALC_ID) {
	if (calc_func_mode == CALC_CODE) {
    	  add_annot_code(have_ann, *sp0_p, *sp1_p, *sp1a_p,
#ifndef TFAST
			 i0_offset+seq_pos(i0,aln->qlrev,0),
			 i1_offset+seq_pos(i1,aln->llrev,0),
#else
			 i1_offset+seq_pos(i1,aln->qlrev,0),
			 i0_offset+seq_pos(i0,aln->llrev,0),
#endif
			 sim_sym[*spa_p], annot_var_dyn);
	}
	display_push_features(annot_stack, annot_var_dyn,
#ifndef TFAST
			      i0_offset+seq_pos(i0,aln->qlrev,0), *sp0_p,
			      i1_offset+seq_pos(i1,aln->llrev,0), *sp1_p,
#else
			      i1_offset+seq_pos(i1,aln->qlrev,0), *sp0_p,
			      i0_offset+seq_pos(i0,aln->llrev,0), *sp1_p,
#endif
			      sim_sym[*spa_p],
			      a_res->rst.score[ppst->score_ix], a_res->rst.comp, a_res->sw_score,
			      n0, n1, pstat_void, annot_fmt);
	have_push_features = 0;
      }

      if (cumm_seq_score) *i_spa++ = itmp;

      i1++;
      if (calc_func_mode == CALC_CONS) {sp0_p++; sp1_p++; spa_p++;}
      lenc++;
      break;
    case 4:		/* frame shift, -1, then match */
      nfs++;
      i0 += 4;

      *sp0_p = '\\';
      *sp1_p = '-';
      *spa_p = M_DEL;

      if (calc_func_mode == CALC_CODE) {
#ifndef TFAST
        update_code(align_code_dyn, update_data_p, 4, *spa_p,*sp0_p,*sp1_p);
#else
        update_code(align_code_dyn, update_data_p, 4, *spa_p,*sp1_p,*sp0_p);
#endif
      }

      if (calc_func_mode == CALC_CONS) {sp0_p++; sp1_p++; spa_p++;}

      if (cumm_seq_score) *i_spa++ = ppst->gshift;

      if (have_ann && calc_func_mode == CALC_CONS) {*sp1a_p++ = *sp0a_p++ = ' ';}
      not_c++;

      *sp1_p = sq[aap=ap1[i1]];
      *sp0_p = f_str->weight_c[aap][ap0[i0]].c4;
      itmp = ppst->pam2[0][pascii[*sp0_p]][aap];

      if (have_ann) {
	if (calc_func_mode != CALC_ID  && calc_func_mode != CALC_ID_DOM) {
	  *sp1a_p = ann_arr[ap1a[i1]];
	  *sp0a_p = ' ';
	}
	if (s_annotp_arr_p && (i1+i1_offset == s_annotp_arr_p[i1_annot]->pos || i1+i1_offset == i1_left_end)) {
	  i1_annot = next_annot_match(&itmp, ppst->pam2[0][pascii[*sp0_p]],
#ifndef TFAST
					i1_offset+seq_pos(i1,aln->llrev,0),	/* annotated target (prot) coordinate */
					i0_offset+seq_pos(i0,aln->qlrev,0),
#else
					i1_offset+seq_pos(i1,aln->qlrev,0),
					i0_offset+seq_pos(i0,aln->llrev,0),
#endif
				      sp1_p, sp1a_p, sq,
				      i1_annot, annotp_p->n_annot, s_annotp_arr_p, &ann_comment,
				      annot_stack, have_push_features_p, &v_delta,
				      &d1_score, &d1_ident, &d1_alen, &d1_gaplen,
				      &left_domain_head1, left_domain_list1, &i1_left_end,0);

	  /* must be out of the loop to capture the last value */
	  if (ppst->sq[ap1[i1]] != *sp1_p) {
	    t_spa = align_type(itmp, *sp0_p, *sp1_p, 0, NULL, ppst->pam_x_id_sim);
	    if (calc_func_mode != CALC_ID  && calc_func_mode != CALC_ID_DOM) {
	      comment_var(
#ifndef TFAST
			    i0_offset+seq_pos(i0,aln->qlrev,0), *sp0_p,
			    i1_offset+seq_pos(i1,aln->llrev,0), *sp1_p,
#else
			    i1_offset+seq_pos(i1,aln->qlrev,0), *sp0_p,
			    i0_offset+seq_pos(i0,aln->llrev,0), *sp1_p,
#endif
			  sq[ap1[i1]], sim_sym[t_spa], ann_comment,
			  annot_var_dyn, comment_target, annot_fmt);
	    }
	    else {
	      sprintf(tmp_str,"%c%d%c;",ppst->sq[ap1[i1]],i1+1,*sp1_p);
	      /*   SAFE_STRNCAT(annot_var_s,tmp_str,n_annot_var_s); */
	      dyn_strcat(annot_var_dyn, tmp_str);
	    }
	  }
	  d1_score += ppst->gshift;
	  d1_score += itmp;
	  prev_match = 1;
	}
	if (calc_func_mode == CALC_CONS) {sp0a_p++; sp1a_p++;}
      }

      *spa_p = align_type(itmp, *sp0_p, *sp1_p, 0, aln, ppst->pam_x_id_sim);

      if (calc_func_mode == CALC_CODE) {
#ifndef TFAST
	update_code(align_code_dyn, update_data_p, 3, *spa_p,*sp0_p,*sp1_p);
#else
	update_code(align_code_dyn, update_data_p, 3, *spa_p,*sp1_p,*sp0_p);
#endif
      }

      d1_alen++;
      if (*spa_p == M_IDENT) {d1_ident++;}
      if (cumm_seq_score) *i_spa++ = itmp;

      if (have_ann && have_push_features &&  calc_func_mode != CALC_ID) {
	display_push_features(annot_stack, annot_var_dyn,
#ifndef TFAST
			      i0_offset+seq_pos(i0,aln->qlrev,0), *sp0_p,
			      i1_offset+seq_pos(i1,aln->llrev,0), *sp1_p,
#else
			      i1_offset+seq_pos(i1,aln->qlrev,0), *sp0_p,
			      i0_offset+seq_pos(i0,aln->llrev,0), *sp1_p,
#endif
			      sim_sym[*spa_p], 
			      a_res->rst.score[ppst->score_ix], a_res->rst.comp, a_res->sw_score,
			      n0, n1, pstat_void, annot_fmt);
	have_push_features = 0;
      }

      i1++;
      if (calc_func_mode == CALC_CONS) {sp0_p++; sp1_p++; spa_p++;}
      lenc++;
      break;
    case 5:		/* insertion in 0 */
      if (have_ann && calc_func_mode == CALC_CONS) {
	*sp1a_p++ = *sp0a_p++ = ' ';
      }

      if (cumm_seq_score) {
	if (prev_match) *i_spa = ppst->gdelval;
	*i_spa++ = ppst->ggapval;
      }

      if (prev_match) d1_score += ppst->gdelval;
      d1_score += ppst->ggapval;

      prev_match = 0;

      i0 += 3;
      *sp0_p = f_str->weight_c[0][ap0[i0]].c3;
      *sp1_p = '-';
      *spa_p = M_DEL;

      if (calc_func_mode == CALC_CODE) {
	*spa_p = 5;
#ifndef TFAST
	update_code(align_code_dyn, update_data_p, 5, *spa_p,*sp0_p,*sp1_p);
#else
	update_code(align_code_dyn, update_data_p, 5, *spa_p,*sp1_p,*sp0_p);
#endif
      }

      lenc++;
      ngap_p++;

      if (calc_func_mode == CALC_CONS) {sp0_p++; sp1_p++; spa_p++;}
      break;
    case 0:		/* insertion in 1 */

      *sp0_p = '-';
      *sp1_p = sq[ap1[i1]];
      *spa_p = M_DEL;

      if (calc_func_mode == CALC_CODE) {
	*spa_p = 5; /* indel code */
#ifndef TFAST
	update_code(align_code_dyn, update_data_p, 0, *spa_p,*sp0_p,*sp1_p);
#else
	update_code(align_code_dyn, update_data_p, 0, *spa_p,*sp1_p,*sp0_p);
#endif	
      }

      if (cumm_seq_score) {
	if (prev_match) *i_spa = ppst->gdelval;
	*i_spa++ += ppst->gdelval;
      }

      if (calc_func_mode == CALC_CONS) {sp0_p++; sp1_p++; spa_p++;}

      if (have_ann) {
	if (calc_func_mode != CALC_ID  && calc_func_mode != CALC_ID_DOM) {
	  *sp0a_p = ' ';
	  *sp1a_p = ann_arr[ap1a[i1]];
	}

	if (s_annotp_arr_p) {
	  /* coordiates are much more complex for next_annot_match,
	     and comment_var, because they may need to be reversed */

	  if (i1 + i1_offset == s_annotp_arr_p[i1_annot]->pos) {
	    i1_annot = next_annot_match(&itmp, ppst->pam2[0][ap0[i0]],
#ifndef TFAST
					i1_offset+seq_pos(i1,aln->llrev,0),	/* annotated target (prot) coordinate */
					i0_offset+seq_pos(i0,aln->qlrev,0),
#else
					i1_offset+seq_pos(i1,aln->qlrev,0),
					i0_offset+seq_pos(i0,aln->llrev,0),
#endif
					sp1_p, sp1a_p, sq, 
					i1_annot, annotp_p->n_annot, s_annotp_arr_p,
					&ann_comment, annot_stack, have_push_features_p, &v_delta,
					&d1_score, &d1_ident, &d1_alen,  &d1_gaplen,
					&left_domain_head1, left_domain_list1, &i1_left_end,0);

	  }

	  if (prev_match) d1_score += ppst->gdelval;
	  d1_score += ppst->ggapval;
	  d1_alen++;
	  d1_gaplen++;
	  prev_match = 0;
	}
	if (calc_func_mode == CALC_CONS) {sp0a_p++; sp1a_p++;}
      }

      if (have_ann && have_push_features  &&  calc_func_mode != CALC_ID) {
	display_push_features(annot_stack, annot_var_dyn,
#ifndef TFAST
			      i0_offset+seq_pos(i0,aln->qlrev,0), *sp0_p,
			      i1_offset+seq_pos(i1,aln->llrev,0), *sp1_p,
#else
			      i1_offset+seq_pos(i1,aln->qlrev,0), *sp0_p,
			      i0_offset+seq_pos(i0,aln->llrev,0), *sp1_p,
#endif
			      sim_sym[*spa_p], 
			      a_res->rst.score[ppst->score_ix], a_res->rst.comp, a_res->sw_score,
			      n0, n1, pstat_void, annot_fmt);
	have_push_features = 0;
      }

      i1++;
      lenc++;
      ngap_d++;
      break;
    }
  }

  /* done with alignment loop */

  if (calc_func_mode == CALC_CODE) {
    close_update_data(align_code_dyn, update_data_p);
  }

  if (have_ann) {
    if (calc_func_mode != CALC_ID  && calc_func_mode != CALC_ID_DOM) {*sp0a_p = *sp1a_p = '\0';}
    if (s_annotp_arr_p) {
      have_push_features = 0;

      if (s_annotp_arr_p && i1_left_end > 0) {
	close_annot_match(-1, annot_stack, have_push_features_p,
			  &d1_score, &d1_ident, &d1_alen,  &d1_gaplen,
			  &left_domain_head1, &i1_left_end,
			  0);
      }

      if (have_push_features  &&  calc_func_mode != CALC_ID) {
	display_push_features(annot_stack, annot_var_dyn,
#ifndef TFAST
			      i0_offset+seq_pos(i0,aln->qlrev,0), *sp0_p,
			      i1_offset+seq_pos(i1,aln->llrev,0), *sp1_p,
#else
			      i1_offset+seq_pos(i1,aln->qlrev,0), *sp0_p,
			      i0_offset+seq_pos(i0,aln->llrev,0), *sp1_p,
#endif
			      sim_sym[*spa_p], 
			      a_res->rst.score[ppst->score_ix], a_res->rst.comp, a_res->sw_score,
			      n0, n1, pstat_void, annot_fmt);
	have_push_features = 0;
      }
    }
    if (left_domain_list1) free(left_domain_list1);
    free_stack(annot_stack);
  }
  *spa_p = '\0';

#ifndef TFAST
  aln->amax0 = i0;
  aln->amax1 = i1;
  aln->ngap_q = ngap_d;
  aln->ngap_l = ngap_p;
#else
  aln->amax1 = i0;
  aln->amax0 = i1;
  aln->ngap_q = ngap_p;
  aln->ngap_l = ngap_d;
#endif
  aln->calc_last_set = 1;

  aln->nfs = nfs;

  *score_delta = v_delta;

  if (lenc < 0) lenc = 1;
  *nc = lenc;
/*	now we have the middle, get the right end */
  return lenc+not_c;
}

int
calc_cons_a(const unsigned char *aa0, int n0,
	    const unsigned char *aa1, int n1,
	    int *nc,
	    struct a_struct *aln,
	    struct a_res_str *a_res, 
	    struct pstruct *ppst,
	    char *seqc0, char *seqc1, char *seqca, int *cumm_seq_score,
	    const unsigned char *ann_arr,
	    const unsigned char *aa0a, const struct annot_str *annot0_p, char *seqc0a,
	    const unsigned char *aa1a, const struct annot_str *annot1_p, char *seqc1a,
	    int *score_delta,
	    struct dyn_string_str *annot_var_dyn,
	    struct f_struct *f_str,
	    void *pstat_void)
{
  return calc_cons_u(
		     aa0, n0, aa1, n1,
		     a_res, ppst, f_str, pstat_void,
		     ann_arr, aa0a, annot0_p, aa1a, annot1_p, CALC_CONS, 0,
		     nc, seqc0, seqc1, seqca, cumm_seq_score,
		     seqc0a, seqc1a, aln, score_delta, annot_var_dyn, NULL
		     );
}

void
calc_astruct(struct a_struct *aln_p, struct a_res_str *a_res_p, struct f_struct *f_str) {

  aln_p->calc_last_set = 0;

#ifndef TFAST	/* FASTX */
  aln_p->amin1 = a_res_p->min0;	/* prot */
  aln_p->amin0 = a_res_p->min1;	/* DNA */
  aln_p->amax1 = a_res_p->max0;	/* prot */
  aln_p->amax0 = a_res_p->max1;	/* DNA */
#else		/* TFASTX */
  aln_p->amin0 = a_res_p->min0;	/* DNA */
  aln_p->amin1 = a_res_p->min1;	/* prot */
  aln_p->amax0 = a_res_p->max0;	/* DNA */
  aln_p->amax1 = a_res_p->max1;	/* prot */
#endif
}

/* build an array of match/ins/del - length strings */

/* modified 10-June-2014 to distinguish matches from mismatches, op=1
   (previously unused) indicates an aligned non-identity */

/* op_codes are:  0 - aa insertion
   		  1 - (now) aligned non-identity
   		  2 - -1 frameshift
   		  3 - aligned identity
   		  4 - +1 frameshift
   		  5 - codon insertion
*/

static struct update_code_str *
init_update_data(int show_code) {

  struct update_code_str *update_data_p;

  if ((update_data_p = (struct update_code_str *)calloc(1,sizeof(struct update_code_str)))==NULL) {
    fprintf(stderr,"*** ERROR [%s:%d] - init_update_data(): cannot allocate update_code_str\n",
	      __FILE__, __LINE__);
    return NULL;
  }

  update_data_p->p_op_idx = -1;
  update_data_p->p_op_cnt = 0;
  update_data_p->show_code = show_code;
  update_data_p->btop_enc = 0;

  if ((show_code & SHOW_CODE_CIGAR) == SHOW_CODE_CIGAR) {
    update_data_p->op_map = cigar_code;
    update_data_p->cigar_order = 1;
  }
  else if ((show_code & SHOW_CODE_BTOP) == SHOW_CODE_BTOP) {
    update_data_p->op_map = ori_code;
    update_data_p->cigar_order = 0;
    update_data_p->btop_enc = 1;
  }    
  else {
    update_data_p->op_map = ori_code;
    update_data_p->cigar_order = 0;
  }

  if ((show_code & SHOW_CODE_EXT) == SHOW_CODE_EXT) {
    update_data_p->show_ext = 1;
  }
  else {
    update_data_p->show_ext = 0;
  }

  return update_data_p;
}

static void
close_update_data(struct dyn_string_str *align_code_dyn,
		  struct update_code_str *up_dp) {
  char tmp_cnt[MAX_SSTR];

  if (!up_dp) return;

  if (up_dp->p_op_cnt) {
    if (up_dp->btop_enc) {
      sprintf(tmp_cnt,"%d",up_dp->p_op_cnt);
      up_dp->p_op_cnt = 0;
    }
    else {
      sprintf_code(tmp_cnt,up_dp, up_dp->p_op_idx, up_dp->p_op_cnt);
    }
    dyn_strcat(align_code_dyn, tmp_cnt);
  }

  free(up_dp);
}

/* update_indel_code() has been modified to work more correctly with
   ggsearch/glsearch, which, because alignments can start with either
   insertions or deletions, can produce an initial code of "0=".  When
   that happens, it is ignored and no code is added.

   *align_code_dyn - alignment string (dynamic)
   op -- encoded operation, currently 0=match, 1-delete, 2-insert, 3-term-match, 4-mismatch
   op_cnt -- length of run
   show_code -- SHOW_CODE_CIGAR uses cigar_code, otherwise legacy
*/

/* update_indel_code() is called for insertions and deletions
   update_match_code() is called for every match
*/

static void
sprintf_code(char *tmp_str, struct update_code_str *up_dp, int op_idx, int op_cnt) {

  if (op_cnt == 0) return;

  if (up_dp->cigar_order) {
    sprintf(tmp_str,"%d%c",op_cnt,up_dp->op_map[op_idx]);
  }
  else {
    sprintf(tmp_str,"%c%d",up_dp->op_map[op_idx],op_cnt);
  }
}

/* only called for btop alignment encoding, for identity, update
   count, otherwise, print previous count and current difference.
   assumes that up_dp->p_op_cnt only tracks identity

   for fx/fz, op=0, 
*/

static void
sprintf_btop(char *tmp_str, 
	     struct update_code_str *up_dp, 
	     int op, int sim_code,
	     unsigned char sp0, unsigned char sp1)
{
  char local_str[MAX_SSTR];
  local_str[0]='\0';

  tmp_str[0] = '\0';

  /* only aligned identities update counts */
  if (op==3 && sim_code == M_IDENT) {
    if ((sp0 == '*' && (sp1 == '*' || toupper(sp1) == 'U'))
	|| (sp1 == '*' && (sp0 == '*' || toupper(sp0) == 'U'))) {
      if (up_dp->p_op_cnt > 0) {
	sprintf(tmp_str,"%d**",up_dp->p_op_cnt);
	up_dp->p_op_cnt = 0;
	return;
      }
    }
    else {
      up_dp->p_op_cnt++;
      return;
    }
  }
  else {
    if (up_dp->p_op_cnt > 0) {
      sprintf(local_str,"%d",up_dp->p_op_cnt);
    }
    up_dp->p_op_cnt = 0;
    sprintf(tmp_str,"%s%c%c",local_str,sp0,sp1);
  }
}

static void
update_code(struct dyn_string_str *align_code_dyn,
	    struct update_code_str *up_dp, int op, 
	    int sim_code,  unsigned char sp0, unsigned char sp1)
{
  char tmp_cnt[MAX_SSTR];
  tmp_cnt[0]='\0';

  if (up_dp->btop_enc) {
    sprintf_btop(tmp_cnt, up_dp, op, sim_code, sp0, sp1);
    dyn_strcat(align_code_dyn, tmp_cnt);
    return;
  }

  /* there are two kinds of "op's", one time and accumulating */
  /* op == 2, 4 are one-time: */

  switch (op) {
  case 2:
  case 4:
    sprintf_code(tmp_cnt,up_dp, up_dp->p_op_idx,up_dp->p_op_cnt);
    dyn_strcat(align_code_dyn, tmp_cnt);
    sprintf_code(tmp_cnt,up_dp, op, 1);
    dyn_strcat(align_code_dyn, tmp_cnt);
    up_dp->p_op_cnt = 0;
    break;
  case 0:
  case 5:
    if (op == up_dp->p_op_idx) {
      up_dp->p_op_cnt++;
    }
    else {
      sprintf_code(tmp_cnt,up_dp, up_dp->p_op_idx,up_dp->p_op_cnt);
      dyn_strcat(align_code_dyn, tmp_cnt);
      up_dp->p_op_idx = op;
      up_dp->p_op_cnt = 1;
    }
    break;
  case 1:
  case 3:
    if (sp0 != '*' && sp1 != '*') {	/* default case, not termination */
      if (up_dp->show_ext) {
	if (sim_code != M_IDENT) { op = 1;}
      }
    }
    else {	/* have a termination codon, output for !SHOW_CODE_CIGAR */
      if (!up_dp->cigar_order) {
	if (sp0 == '*' || sp1 == '*') { op = 6;}
      }
      else if (up_dp->show_ext && (sp0 != sp1)) { op = 1;}
    }

    if (up_dp->p_op_cnt == 0) {
      up_dp->p_op_idx = op;
      up_dp->p_op_cnt = 1;
    }
    else if (op != up_dp->p_op_idx) {
      sprintf_code(tmp_cnt,up_dp, up_dp->p_op_idx,up_dp->p_op_cnt);
      dyn_strcat(align_code_dyn, tmp_cnt);
      up_dp->p_op_idx = op;
      up_dp->p_op_cnt = 1;
    }
    else {
      up_dp->p_op_cnt++;
    }
    break;
  }
  return;
}

int calc_code(const unsigned char *aa0, int n0,
	      const unsigned char *aa1, int n1,
	      struct a_struct *aln,
	      struct a_res_str *a_res,
	      struct pstruct *ppst,
	      struct dyn_string_str *align_code_dyn,
	      const unsigned char *ann_arr, 
	      const unsigned char *aa0a,
	      const struct annot_str *annot0_p,
	      const unsigned char *aa1a,
	      const struct annot_str *annot1_p,
	      struct dyn_string_str *annot_code_dyn,
	      int *score_delta,
	      struct f_struct *f_str,
	      void *pstat_void,
	      int display_code)
{
  int nc;

  return calc_cons_u(
		     aa0, n0, aa1, n1,
		     a_res, ppst, f_str, pstat_void,
		     ann_arr, aa0a, annot0_p, aa1a, annot1_p, CALC_CODE, 
		     display_code,
		     &nc, NULL, NULL, NULL, NULL,
		     NULL, NULL, aln, score_delta, annot_code_dyn,
		     align_code_dyn
		     );
}

/* calc_id never looks at domains or features, only variation */

int calc_id(const unsigned char *aa0, int n0,
	    const unsigned char *aa1, int n1,
	    struct a_struct *aln, 
	    struct a_res_str *a_res,
	    struct pstruct *ppst,
	    const struct annot_str *annot0_p,
	    const struct annot_str *annot1_p,
	    int *score_delta,
	    struct dyn_string_str *annot_var_dyn,
	    struct f_struct *f_str)
{
  int nc;

  return calc_cons_u(
		     aa0, n0, aa1, n1,
		     a_res, ppst, f_str, NULL,
		     NULL, NULL, annot0_p, NULL, annot1_p, CALC_ID, 0,
		     &nc, NULL, NULL, NULL, NULL,
		     NULL, NULL, aln, score_delta, annot_var_dyn,
		     NULL
		     );
}

/* calc_idd looks at domains and variation */

int calc_idd(const unsigned char *aa0, int n0,
	    const unsigned char *aa1, int n1,
	    struct a_struct *aln, 
	    struct a_res_str *a_res,
	    struct pstruct *ppst,
	    const struct annot_str *annot0_p,
	    const struct annot_str *annot1_p,
	    int *score_delta,
	    struct dyn_string_str *annot_var_dyn,
	    struct f_struct *f_str)
{
  int nc;

  return calc_cons_u(
		     aa0, n0, aa1, n1,
		     a_res, ppst, f_str, NULL,
		     NULL, NULL, annot0_p, NULL, annot1_p, CALC_ID_DOM, 0,
		     &nc, NULL, NULL, NULL, NULL,
		     NULL, NULL, aln, score_delta, annot_var_dyn,
		     NULL
		     );
}
