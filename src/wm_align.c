/* $Id: wm_align.c 1166 2013-05-30 01:05:55Z wrp $  */

/* algorithms and code provided by Webb Miller, Penn State
   University */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "defs.h"
#include "rstruct.h"
#include "aln_structs.h"

struct swstr {int H, E;};

int
NW_ALIGN(int IW, const unsigned char *B,
	 int M, int N,
	 int **W, int G, int H, int *S, int *NC);

static int
CHECK_SCORE(int IW, const unsigned char *B,
	    int M, int N,
	    int *S, int **W, int G, int H, int *nres, int *sw);

/* sw_walign() is here, rather than in dropgsw2.c, because it is also
   used by dropnfa.c
*/

struct a_res_str *
merge_ares_chains(struct a_res_str *cur_ares,
		  struct a_res_str *tmp_ares,
		  int score_ix,
		  const char *msg);

int
sw_walign (int **pam2p, int n0,
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
   int qr;
   int     score;
   int cost, I, J, K, L;

   qr = q + r;

   /* initialize 0th row */
   for (ssj=ss; ssj<ss+n0; ssj++) {
     ssj->H = 0;
     ssj->E = -q;
   }

   /* I = saved position in aa1
      J = saved position in aa0
   */
   score = I = J = 0;
   aa1p = aa1;
   i = 0;
   while (*aa1p) {
     h = p = 0;
     f = -q;
     /* pwaa = waa + (*aa1p++ * n0); */
     for (ssj = ss, j=0; j < n0; ssj++, j++) {
       if ((h =   h    - qr) > (f =   f    - r)) f = h;
       if ((h = ssj->H - qr) > (e = ssj->E - r)) e = h;
       /* h = p + *pwaa++; */
       h = p + pam2p[j][*aa1p];
       if (h < 0 ) h = 0;
       if (h < f ) h = f;
       if (h < e ) h = e;
       p = ssj->H;
       ssj->H = h;
       ssj->E = e;
       if (h > score) {
	 score = h;
	 I = i;
	 /* J = (int)(ssj-ss);  */
	 J = j;
       }
     }
     i++;
     aa1p++;
   }				/* done with forward pass */
   if (score <= 0) return 0;

  /* to get the start point, go backwards */

   /* K = begin in aa1
      L = begin in aa0
   */
  cost = K = L = 0;
  for (ssj=ss+J; ssj>=ss; ssj--) {
    ssj->H=ssj->E= -1;
  }

  for (i=I,aa1p=aa1+I; i>=0; i--) {
    h = f = -1;
    p = (i == I) ? 0 : -1;
    for (ssj = ss+J, j=J; ssj>=ss; ssj--,j--) {
      f = max (f,h-q)-r;
      ssj->E=max(ssj->E,ssj->H-q)-r;
      h = max(max(ssj->E,f), p+pam2p[j][aa1[i]]);
      p = ssj->H;
      ssj->H=h;
      if (h > cost) {
	cost = h;
	K = i;
	L = (int)(ssj-ss);
	if (cost >= score) goto found;
      }
    }
  }

found:

  /*  printf(" %d: L: %3d-%3d/%3d; K: %3d-%3d/%3d\n",score,L,J,n0,K,I,n1); */

  a_res->n1 = n1;
  a_res->max0 = J+1; a_res->min0 = L; a_res->max1 = I+1; a_res->min1 = K;

  NW_ALIGN(L,&aa1[K-1],J-L+1,I-K+1,pam2p,q,r,a_res->res,&a_res->nres);

  return score;
}

/* nsw_malign is a recursive interface to nw/sw_walign() that is called
   from do_walign(). nsw_malign() first does an alignment, then checks
   to see if the score is greater than the threshold. If so, it tries
   doing a left and right alignment.

   2009-Mar-22 -- This version generalizes the strategy for
   partitioning the solution by taking the *_walign function as an
   argument

   2009-May-1 -- add code to ensure that returned a_res->chain is
   sorted by score.  One strategy is to simply always insert at the
   appropriate place, which requires re-searching from the top each
   time (which is not a big deal, since the list is short).  We simply
   need another argument to nsw_malign, which is the head of the list.
 */
struct a_res_str *
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
	    )
{
  struct a_res_str *tmpl_ares, *tmpr_ares, *this_ares;
  unsigned char *local_aa1;
  int nc, score_ix;
  int min_alen;
  int max_sub_score = -1;

  min_alen = min(MIN_LOCAL_LEN,n0);

  /* now we need alignment storage - get it */
  if ((cur_ares->res = (int *)calloc((size_t)max_res,sizeof(int)))==NULL) {
    fprintf(stderr," *** cannot allocate alignment results array %d\n",max_res);
    exit(1);
  }

  score_ix = 0;

  cur_ares->next = NULL;

  cur_ares->sw_score = (*fn_walign)(pam2p[0], n0, aa1, n1,
				    gdelval,  ggapval,
				    ss, cur_ares);

  /* The scores in a_res->rst include low-complexity alignment.
     Re-calculate the score of the optimal alignment using the -S matrix.

     This makes sense for secondary HSP's, but not for the initial
     HSP, because it -S score for the non-S alignment could be smaller
     than the original optimal -S score (it cannot be higher).
  */
  CHECK_SCORE(cur_ares->min0, &aa1[cur_ares->min1-1],
	      cur_ares->max0 - cur_ares->min0, cur_ares->max1-cur_ares->min1,
	      cur_ares->res,
	      pam2p[pam_ix], gdelval, ggapval,
	      &nc, &cur_ares->rst.score[0]);

  if (!do_rep || cur_ares->rst.score[score_ix] <= score_thresh) { return cur_ares;}

  if (cur_ares->min1 >= min_alen) { /* try the left  */
    /* allocate a_res */
    tmpl_ares = (struct a_res_str *)calloc(1, sizeof(struct a_res_str));

    local_aa1 = (unsigned char *)calloc(cur_ares->min1+2,sizeof(unsigned char));
    local_aa1++;
    memcpy(local_aa1,aa1,cur_ares->min1);
    /*
    save_res = aa1[cur_ares->min1];
    aa1[cur_ares->min1] = '\0';
    */
    tmpl_ares = nsw_malign(pam2p, pam_ix, n0, local_aa1, cur_ares->min1,
			   score_thresh, max_res,
			   gdelval, ggapval, ss, tmpl_ares,
			   fn_walign, do_rep);
    free(--local_aa1);
    /*
    aa1[cur_ares->min1] = save_res;
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

  if (n1 - cur_ares->max1 >= min_alen) { /* try the right  */
    /* allocate a_res */
    tmpr_ares = (struct a_res_str *)calloc(1, sizeof(struct a_res_str));

    /* find boundaries */

    local_aa1 = (unsigned char *)calloc(n1-cur_ares->max1+2,sizeof(unsigned char));
    local_aa1++;
    memcpy(local_aa1,aa1+cur_ares->max1,n1-cur_ares->max1);
    /*
    save_res = aa1[cur_ares->max1-1];
    aa1[cur_ares->max1-1] = '\0';
    */

    tmpr_ares = nsw_malign(pam2p, pam_ix, n0, local_aa1, n1 - cur_ares->max1,
			   score_thresh, max_res,
			   gdelval, ggapval, ss, tmpr_ares,
			   fn_walign, do_rep);

    free(--local_aa1);
    /*
    aa1[cur_ares->max1-1] = save_res;
    */

    if (tmpr_ares->rst.score[score_ix] > score_thresh) {
      /* adjust the left boundary */
      for (this_ares = tmpr_ares; this_ares; this_ares = this_ares->next) {
	this_ares->min1 += cur_ares->max1;
	this_ares->max1 += cur_ares->max1;
      }

      if (tmpr_ares->rst.score[score_ix] >= max_sub_score) {
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

  /* We have checked both left and right, and better score is in max_sub_score.
     If both scores are <= score_thresh, then forget it */

  if (max_sub_score <= score_thresh) {
    if (tmpl_ares) {
      if (tmpl_ares->res) {free(tmpl_ares->res);}
      free(tmpl_ares);
    }
    if (tmpr_ares) {
      if (tmpr_ares->res) {free(tmpr_ares->res);}
      free(tmpr_ares);
    }
    return cur_ares;
  }

  cur_ares =  merge_ares_chains(cur_ares, tmpl_ares, score_ix, "left");
  cur_ares =  merge_ares_chains(cur_ares, tmpr_ares, score_ix, "right");

  return cur_ares;
}

#define gap(k)  ((k) <= 0 ? 0 : q+r*(k))	/* k-symbol indel cost */

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

/* align(A,B,M,N,tb,te,last) returns the cost of an optimum conversion between
   A[1..M] and B[1..N] that begins(ends) with a delete if tb(te) is zero
   and appends such a conversion to the current script.                   */

static int
nw_align(int iw,	/* beginning of alignment in pam2p profile */
	 const unsigned char *B,	/* second sequence aa1 */
	 int M, int N,			/* length of profile, aa1 */
	 int tb, int te,
	 int **w, int q, int r, 	/* pam2p profile, open, ext */
	 struct swstr *f_ss,		/* forward, reverse row matrix */
	 struct swstr *r_ss,
	 int dir,			/* dir [0..3] is not currently used */
	 int **sapp, int *last)
{

  int midi, midj, type;	/* Midpoint, type, and cost */
  int midc;
  int c1, c2;

  register int   i, j;
  register int c, e, d, s;
  int qr, t, *wa;

/*   print_seq_prof(A,M,B,N,w,iw); */

/*  m = g + h;  */
  qr = q + r;

/* Boundary cases: M <= 1 or N == 0 */

  if (N <= 0) {
    if (M > 0) {DEL(M)}
    return -gap(M);
  }

  if (M <= 1) {
    if (M <= 0) {
      INS(N)
      return -gap(N);
    }

    if (tb < te) tb = te;
    midc = (tb-r) - gap(N);
    midj = 0;
/*  wa = w[A[1]]; */
    wa = w[iw];
    for (j = 1; j <= N; j++) {
      c = -gap(j-1) + wa[B[j]] - gap(N-j);
      if (c > midc) { midc = c; midj = j;}
    }
    if (midj == 0) { DEL(1) INS(N) }
    else  {
      if (midj > 1) { INS(midj-1)}
      *last = (*sapp)[0] = 0;
      (*sapp)++;
      if (midj < N) { INS(N-midj)}
    }
    return midc;
  }

/* Divide: Find optimum midpoint (midi,midj) of cost midc */

  midi = M/2;		/* Forward phase:                          */
  f_ss[0].H = 0;	/*   Compute H(M/2,k) & E(M/2,k) for all k */
  f_ss[0].E = t = -q;
  for (j = 1; j <= N; j++) {
    f_ss[j].H = t = t-r;
    f_ss[j].E = t-q;
  }
  t = tb;
  for (i = 1; i <= midi; i++) {
    s = f_ss[0].H;
    f_ss[0].H = c = t = t-r;
    e = t-q;
/*    wa = w[A[i]]; */
    wa = w[iw+i-1];
    for (j = 1; j <= N; j++) {
      if ((c =   c   - qr) > (e =   e   - r)) e = c;
      if ((c = f_ss[j].H - qr) > (d = f_ss[j].E - r)) d = c;
      c = s + wa[B[j]];
      if (e > c) c = e;
      if (d > c) c = d;
      s = f_ss[j].H;
      f_ss[j].H = c;
      f_ss[j].E = d;
    }
  }
  f_ss[0].E = f_ss[0].H;

  r_ss[N].H = 0;		/* Reverse phase:                  */
  t = -q;			/*   Compute R(M/2,k) & S(M/2,k) for all k */

  for (j = N-1; j >= 0; j--) {
    r_ss[j].H = t = t-r;
    r_ss[j].E = t-q;
  }

  t = te;
  for (i = M-1; i >= midi; i--) {
    s = r_ss[N].H;
    r_ss[N].H = c = t = t-r;
    e = t-q;
/*    wa = w[A[i+1]]; */
    wa = w[iw+i];
    for (j = N-1; j >= 0; j--) {
      if ((c =   c   - qr) > (e =   e   - r)) { e = c; }
      if ((c = r_ss[j].H - qr) > (d = r_ss[j].E - r)) { d = c; }
      c = s + wa[B[j+1]];
      if (e > c) c = e;
      if (d > c) c = d;
      s = r_ss[j].H;
      r_ss[j].H = c;
      r_ss[j].E = d;
    }
  }
  r_ss[N].E = r_ss[N].H;

  midc = f_ss[0].H+r_ss[0].H;		/* Find optimal midpoint */
  midj = 0;
  type = 1;

  for (j = 0; j <= N; j++) {
    if ((c = f_ss[j].H + r_ss[j].H) >= midc) {
      if (c > midc || (f_ss[j].H != f_ss[j].E && r_ss[j].H == r_ss[j].E)) {
	midc = c;
	midj = j;
      }
    }
  }

  for (j = N; j >= 0; j--) {
    if ((c = f_ss[j].E + r_ss[j].E + q) > midc) {
      midc = c;
      midj = j;
      type = 2;
    }
  }

/* Conquer: recursively around midpoint */

  if (type == 1) {
    c1 = nw_align(iw,B,midi,midj,tb,-q,w,q,r,f_ss, r_ss,0,sapp,last);
    c2 = nw_align(iw+midi,B+midj,M-midi,N-midj,-q,te,w,q,r,f_ss, r_ss,1,sapp,last);
  }
  else {
    nw_align(iw,B,midi-1,midj,tb,0,w,q,r,f_ss, r_ss,2,sapp,last);
    DEL(2);
    nw_align(iw+midi+1,B+midj,M-midi-1,N-midj,0,te,w,q,r,f_ss,r_ss,3,sapp,last);
  }
  return midc;
}

/* Interface and top level of comparator */

int
NW_ALIGN(int IW, const unsigned char *B,
	 int M, int N,
	 int **W, int G, int H, int *S, int *NC)
{
  struct swstr *f_ss, *r_ss;
  int *sapp, last;
  int c, ck, sw;

  sapp = S;
  last = 0;

   if ((f_ss = (struct swstr *) calloc (N+2, sizeof (struct swstr)))
       == NULL) {
     fprintf (stderr, " *** cannot allocate f_ss array %3d\n", N+2);
     exit (1);
   }
   f_ss++;

   if ((r_ss = (struct swstr *) calloc (N+2, sizeof (struct swstr)))
       == NULL) {
     fprintf (stderr, " *** cannot allocate r_ss array %3d\n", N+2);
     exit (1);
   }
   r_ss++;

  /*   print_seq_prof(A,M,W,IW); */
  c = nw_align(IW,B,M,N,-G,-G,W,G,H,f_ss, r_ss,0,&sapp,&last);	/* OK, do it */

  ck = CHECK_SCORE(IW,B,M,N,S,W,G,H,NC, &sw);
  if (c != ck) {
    fprintf(stderr," *** Check_score error. %d != %d ***\n",c,ck);
  }

  f_ss--; r_ss--;
  free(r_ss); free(f_ss);

  return c;
}

/* CHECK_SCORE - return the score of the alignment stored in S */

static int
CHECK_SCORE(int iw, const unsigned char *B,
	    int M, int N,
	    int *S, int **w,
	    int g, int h, int *NC, int *sw_score)
{
  register int   i,  j, op, nc;
  int itmp;
  int score;
  int l_score, mx_l_score;

  /*  print_seq_prof(A,M,w,iw); */

  score = i = j = nc = l_score = mx_l_score = 0;
#ifdef SHOW_ALIGN_SCORE
  printf("====start\n");
  printf("#i j pam2 score l_score mx_l_score\n");
#endif
  while (i < M || j < N) {
    op = *S++;
    if (op == 0) {
      itmp = w[iw+i][B[++j]];
      score += itmp;
      i++;
      nc++;
      l_score += itmp;
      if (l_score < 0) l_score = 0;
      if (l_score > mx_l_score) mx_l_score = l_score;
#ifdef SHOW_ALIGN_SCORE
      printf("%d\t%d\t%d\t%d\t%d\t%d\n",i, j, itmp, score, l_score, mx_l_score);
#endif
    }
    else if (op > 0) {
      score = score - (g+op*h);
      j += op;
      nc += op;
      l_score -= (g+op*h);
      if (l_score < 0) l_score = 0;
#ifdef SHOW_ALIGN_SCORE
      printf("%d\t%d\t%d\t%d\t%d\t%d\n",i, j, -(g+op*h) ,score, l_score, mx_l_score);
#endif
    } else {
      score = score - (g-op*h);
      i -= op;
      nc -= op;
      l_score -= (g-op*h);
      if (l_score < 0) l_score = 0;
#ifdef SHOW_ALIGN_SCORE
      printf("%d\t%d\t%d\t%d\t\%d\t%d\n",i, j, -(g-op*h), score, l_score, mx_l_score);
#endif
    }
  }
#ifdef SHOW_ALIGN_SCORE
  printf("%d\t%d\tend\t%d\t%d\n====\n",i, j, score, mx_l_score);
#endif
  *NC = nc;
  /* used to return mx_l_score, which is wrong when CHECK_SCORE is used for global alignments */
#ifndef GGSEARCH
  *sw_score = mx_l_score;
#else
  *sw_score = score;
#endif
  return score;
}

