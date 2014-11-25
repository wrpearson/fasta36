/*
  lsim4.c - calculate non-overlapping local alignments
  
  derived from lsim2.c from Webb Miller
*/

/* $Id: lsim4.c 1070 2012-10-12 16:44:16Z wrp $ */
/* $Revision: 1070 $  */

/*  March 2007 - modified to avoid global references  */

/* October, 2008 - modified following changes from Xiaoqui Huang to
   prevent alignments from crossing the identity diagonal during
   self-comparison */

/* October 27, 2008 - modified to free pair_ptr memory more reliably
   (it is still imperfect) */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "defs.h"
#include "param.h"

#include "lsim4.h"

/* SIM(A,B,M,N,mini_score,Q,R) reports best non-intersecting alignments with
   score >= mini_score of the segments of A and B in order of similarity
   scores, where pam2[a][b] is the score of aligning a and b, and
   -(Q+R*i) is the score of an i-symbol indel.  */

/* SIM uses A[1..M], B[1..N] FORTRAN indexing */

void SIM(const unsigned char *A, /* seq1 indexed A[1..M] */
	 const unsigned char *B, /* seq2 indexed B[1..N] */
	 int M, int N,		 /* len seq1, seq2 */
	 struct pstruct *ppst,	/* parameters */
	 int nseq,		 /* nseq - number of different sequences */
	 int mini_score,	 /* cut-off score */
	 int max_count,		 /* number of alignments */
	 struct a_res_str *a_res)	/* alignment result structure */
{
  int endi, endj, stari, starj;	 /* endpoint and startpoint */ 
  int score;  			 /* the max score in LIST */
  int ck;
  int  i;		 /* row and column indices */

  int flag;
  struct a_res_str *cur_ares, *tmp_ares;

  bool first_pass;			
  int q, r, qr;

  int *sapp, last;
  struct vert_str vLIST;

  struct l_struct *l_ptr;
  pair_ptr z, z_save;
  int count;			 /* maximum size of list */	
  vertex_p v_cur;		 /* temporary pointer */

  /* allocate space for all vectors */

  l_ptr = (struct l_struct *)ckalloc(1, sizeof(struct l_struct));

  l_ptr->CCC = ( space_ptr ) ckalloc(N+1,sizeof(space));

  l_ptr->r_ss = (struct lrr_str *) ckalloc(N + 1,sizeof(struct lrr_str));

  l_ptr->c_ss = (struct lcc_str *) ckalloc(M + 1, sizeof(struct lcc_str));

  l_ptr->row = ( pair_ptr * ) ckalloc( M + 1, sizeof(pair_ptr));

  /* set up list for each row */
  if ( nseq == 2 ) {
    for ( i = 1; i <= M; i++ ) {l_ptr->row[i] = NULL;}
  }
  else {
    for ( i = 1; i <= M; i++ ) {
      l_ptr->row[i] = z = (pair_ptr) ckalloc(1, sizeof(pair));
      z->COL = i;
      z->NEXT = NULL;
    }
  }

#ifdef OLD_FASTA_GAP
  q = -(ppst->gdelval - ppst->ggapval);
#else
  q = -ppst->gdelval;
#endif

  r = -ppst->ggapval;

  qr = q + r;

  vLIST.LIST = vLIST.most = NULL;
  vLIST.numnode = 0;

  /* fill in l_ptr->CCC and vLIST */
  big_pass(A,B,M,N,mini_score,ppst->pam2[0],q,r,nseq, &vLIST, l_ptr);
  first_pass= 1;

  /* Report the K best alignments one by one. After each alignment is
     output, recompute part of the matrix. First determine the size
     of the area to be recomputed, then do the recomputation         */

  count = 0;
  while (count < max_count) {
    if ( vLIST.numnode == 0 ) break;	/* no more alignments */

    if (first_pass) {
      cur_ares = a_res;
    }
    else {	/* need a new a_res */
      tmp_ares = (struct a_res_str *)calloc(1, sizeof(struct a_res_str));
      cur_ares->next = tmp_ares;
      cur_ares = tmp_ares;
    }

    /* get the best next alignment */
    v_cur = findmax(&vLIST);

    score = v_cur->SCORE;
    stari = ++v_cur->STARI;
    starj = ++v_cur->STARJ;
    endi = v_cur->ENDI;
    endj = v_cur->ENDJ;

    l_ptr->m1 = v_cur->TOP;
    l_ptr->mm = v_cur->BOT;
    l_ptr->n1 = v_cur->LEFT;
    l_ptr->nn = v_cur->RIGHT;

    l_ptr->rl = endi - stari + 1;
    l_ptr->cl = endj - starj + 1;

    l_ptr->I = stari - 1;
    l_ptr->J = starj - 1;

    /* minimum allocation for alignment */

    sapp = cur_ares->res =(int *)calloc(2*min(l_ptr->rl,l_ptr->cl), sizeof(int));
    last = 0;

    cur_ares->n1 = N;
    cur_ares->sw_score = cur_ares->rst.score[ppst->score_ix] = score;
    cur_ares->min0 = stari-1;
    cur_ares->min1 = starj-1;
    cur_ares->max0 = stari+l_ptr->rl-1;
    cur_ares->max1 = starj+l_ptr->cl-1;
    cur_ares->next = NULL;

    /* produce an alignment, encoded in sapp - equivalent to "align() in dropgsw.c" */
    (void) diff(&A[stari]-1, &B[starj]-1, l_ptr->rl, l_ptr->cl,
		q, q, (nseq == 2),ppst->pam2[0], q, r, &sapp, &last, l_ptr);

    ck = CHECK_SCORE(&A[stari]-1,&B[starj]-1,l_ptr->rl,l_ptr->cl,
		     cur_ares->res,ppst->pam2[0],q,r,&cur_ares->nres);

#ifdef DEBUG
    /* the same errors are produced by Miller and Huang's sim96 code, so I hope 
       this reflects a mistake in CHECK_SCORE */

    if (score != ck) {
      fprintf(stderr,"*** Check_score error: orig %d != %d recons ***\n aa0[%d-%d] : aa1[%d-%d]\n",
	      score,ck, cur_ares->min0, cur_ares->max0, cur_ares->min1, cur_ares->max1);
    }
#endif

    free(v_cur);

    flag = 0;
    if (first_pass && maxi(l_ptr->rl, l_ptr->cl) > maxi(M,N)/4) {
      /*printf("no locate\n");*/
      flag = 1; l_ptr->n1 = l_ptr->m1 = 0; 
    } 
    else {
      locate(A,B,mini_score,ppst->pam2[0],q,r, nseq, &flag, &vLIST, l_ptr);
    }
    if ( flag ) {
      /*printf("small pass\n");*/
      small_pass(A,B,mini_score,ppst->pam2[0],q,r, nseq, &vLIST, l_ptr);
    }
    first_pass= 0;
    count++;
  }
  /* start cleaning up */

  while (vLIST.numnode > 0) {
    v_cur = findmax(&vLIST);
    if (v_cur) free(v_cur);
  }

  for (i=M; i>0; i--) {
    if ((z=l_ptr->row[i]) != NULL) {
      for (z= l_ptr->row[i], z_save = z->NEXT; z_save != NULL; z = z_save) {
	z_save = z->NEXT;
	free(z);
      }
    }
  }

  free(l_ptr->row);
  free(l_ptr->c_ss);
  free(l_ptr->r_ss);
  free(l_ptr->CCC);
  free(l_ptr);
}

/* A big pass to compute classes scoring over K */
/* fills in the l_ptr->CCC structure for further analysis */
/* adds nodes to v_ptr->LIST when appropriate */
/* does not produce alignments */

static void
big_pass(const unsigned char *A,
	 const unsigned char *B,
	 int M, int N,
	 int mini_score,
	 int **pam2,
	 int Q, int R,
	 int nseq,
	 struct vert_str *v_ptr,
	 struct l_struct *l_ptr)
{
  int  i, j;		/* row and column indices */
  int  c;		/* best score at current point */
  int  f;		/* best score ending with insertion */
  int  d;		/* best score ending with deletion */
  int  p;		/* best score at (i-1, j-1) */
  int  ci, cj;		/* end-point associated with c */
  int  fi, fj;		/* end-point associated with f */
  int  pi, pj;		/* end-point associated with p */
  space_ptr sp;
  pair_ptr z;
  int qr;
  int  *va;				/* pointer to v(A[i], B[j]) */

  qr = Q+R;

  /* Compute the matrix and save the best scores in LIST
     CC : the scores of the current row
     RR and EE : the starting point that leads to score CC
     DD : the scores of the current row, ending with deletion
     SS and FF : the starting point that leads to score DD
  */

  /* Initialize the 0 th row */
  for ( sp=&l_ptr->CCC[1], j = 1; j <= N; j++, sp++ ) {
    sp->CC = sp->RR = 0;
    sp->EE = j;
    sp->DD = - (qr);
    sp->SS = 1;
    sp->FF = j;
  }

  for ( i = 1; i <= M; i++) {
    c = 0;				/* Initialize column 0 */
    f = - (qr);
    va = pam2[A[i]];
    ci = fi = i;
    if ( nseq == 2 ) {
      p = 0;			/* score of current row */
      pi = (i - 1);		/* starting point */
      cj = fj = pj = 0;		/* pj starting point */
    }
    else {
      p = l_ptr->CCC[i].CC;/* score of current row */
      pi = l_ptr->CCC[i].RR;	/* starting point */
      pj = l_ptr->CCC[i].EE;	/* starting point */
      cj = fj = i;
    }
    j = (nseq == 2 ? 1: i+1);
    for ( sp = &l_ptr->CCC[j]; sp <= &l_ptr->CCC[N]; j++, sp++) {
		   
      d = sp->DD;
      c = -1;
      /* assign p+va[B[j]] to c if i, j not part of path */
      DIAG(i, j, c, p+va[B[j]])		/* diagonal */
      if (c < 0) {
	p = sp->CC; pi = sp->RR; pj = sp->EE;
	if (f >= 0) {
	  c = f; ci = fi; cj = fj;
	  /* replace ci, cj with sp->SS, sp->FF if c < d */
	  ORDER1(c, ci, cj,  d, sp->SS, sp->FF)
	  sp->CC = c; sp->RR = ci; sp->EE = cj;
	  sp->DD -= R; f-=R;
	} else if (d >= 0) {
	  sp->CC = d; sp->RR = sp->SS; sp->EE = sp->FF;  
	  sp->DD -= R;
	} else {
	  sp->CC = 0; sp->RR=i; sp->EE = j;
	}
      } else { 
	ci = pi; cj = pj;
	ORDER1(c, ci, cj,  f, fi, fj)
	ORDER1(c, ci, cj,  d, sp->SS, sp->FF)
	p = sp->CC;
	sp->CC = c;
	pi = sp->RR;
	sp->RR = ci;
	pj = sp->EE;
	sp->EE = cj;
	f -= R;
	if (c >= qr) {
	  if ( c > mini_score) {	/* add the score into list */
	    addnode(c, ci, cj, i, j, v_ptr);
	  }
	  d -= R; c -= qr;
	  ORDER1(f, fi, fj, c, ci, cj)
	  ORDER1(d, sp->SS, sp->FF, c, ci, cj)
	  sp->DD = d;
	} else {
	  sp->DD -= R;
	}
      }
    }
  }
}

/* Determine the left and top boundaries of the recomputed area */
/* this function is not recursive */

static void
locate(const unsigned char *A,
       const unsigned char *B,
       int mini_score,
       int **pam2, int Q, int R,
       int nseq,
       int *flag_p,
       struct vert_str *v_ptr,
       struct l_struct *l_ptr) {
  int  i, j;		/* row and column indices */
  int  c;		/* best score at current point */
  int  f;		/* best score ending with insertion */
  int  d;		/* best score ending with deletion */
  int  p;		/* best score at (i-1, j-1) */
  int  ci, cj;		/* end-point associated with c */ 
  int  di, dj;
  int  fi, fj;		/* end-point associated with f */
  int  pi, pj;		/* end-point associated with p */

  space_ptr sp;
  pair_ptr z;
  bool  cflag, rflag;			/* for recomputation */
  int  *va;				/* pointer to v(A[i], B[j]) */
  int  limit;				/* the bound on j */
  int qr; 

  qr = Q + R;

  /* Reverse pass
     rows come from CCC
     CC : the scores on the current row
     RR and EE : the endpoints that lead to CC
     DD : the deletion scores 
     SS and FF : the endpoints that lead to DD

     columns come from c_ss[]
     HH : the scores on the current columns
     II and JJ : the endpoints that lead to HH
     WW : the deletion scores
     XX and YY : the endpoints that lead to WW
  */

  for ( j = l_ptr->nn; j >= l_ptr->n1 ; j-- ) {
    l_ptr->CCC[j].CC = 0;
    l_ptr->CCC[j].EE = j;
    l_ptr->CCC[j].DD = - (Q);
    l_ptr->CCC[j].FF = j;
    if ( nseq == 2 || j > l_ptr->mm )
      l_ptr->CCC[j].RR = l_ptr->CCC[j].SS = l_ptr->mm + 1;
    else
      l_ptr->CCC[j].RR = l_ptr->CCC[j].SS = j;
  }

  for ( i = l_ptr->mm; i >= l_ptr->m1; i-- )  {
    c = p = 0;
    f = - (Q);
    ci = fi = i;
    pi = i + 1;
    cj = fj = pj = l_ptr->nn + 1;
    va = pam2[A[i]];
    if ( nseq == 2 || l_ptr->n1 > i ) limit = l_ptr->n1;
    else limit = i + 1;

    for ( j = l_ptr->nn, sp = &l_ptr->CCC[j]; j >= limit ; j--, sp-- ) {
      f = f - R;
      c = c - qr;
      ORDER(f, fi, fj, c, ci, cj)
      c = sp->CC - qr; 
      d = sp->DD - R;
      ORDER(d, sp->SS, sp->FF, c, sp->RR, sp->EE)
      c = 0;
      DIAG(i, j, c, p+va[B[j]])		/* diagonal */
      if ( c <= 0 ) { c = 0; ci = i; cj = j; }
      else  { ci = pi; cj = pj; }
      ORDER1(c, ci, cj, d, sp->SS, sp->FF)
      ORDER1(c, ci, cj, f, fi, fj)
      p = sp->CC;
      sp->CC = c;
      pi = sp->RR;
      pj = sp->EE;
      sp->RR = ci;
      sp->EE = cj;
      sp->DD = d;
      if ( c > mini_score ) *flag_p = 1;
    }

    if ( nseq == 2 || i < l_ptr->n1 ) {
      l_ptr->c_ss[i].HH = l_ptr->CCC[l_ptr->n1].CC;
      l_ptr->c_ss[i].II = l_ptr->CCC[l_ptr->n1].RR;
      l_ptr->c_ss[i].JJ = l_ptr->CCC[l_ptr->n1].EE;
      l_ptr->c_ss[i].WW = f;
      l_ptr->c_ss[i].XX = fi;
      l_ptr->c_ss[i].YY = fj;
    }
  }
      
  for ( l_ptr->rl = l_ptr->m1, l_ptr->cl = l_ptr->n1; ; ) {
    for ( rflag = cflag = 1; ( rflag && l_ptr->m1 > 1 ) || ( cflag && l_ptr->n1 > 1 ) ;  ) {
      if ( rflag && l_ptr->m1 > 1 ) {	/* Compute one row */
	rflag = 0;
	l_ptr->m1--;
	c = p = 0;
	f = - (Q);
	ci = fi = l_ptr->m1;
	pi = l_ptr->m1 + 1;
	cj = fj = pj = l_ptr->nn + 1;
	va = pam2[A[l_ptr->m1]];
	for ( j = l_ptr->nn, sp = &l_ptr->CCC[j]; j >= l_ptr->n1 ; j--, sp-- ) {
	  f = f - R;
	  c = c - qr;
	  ORDER(f, fi, fj, c, ci, cj)
	  c = sp->CC - qr; 
	  ci = sp->RR;
	  cj = sp->EE;
	  d = sp->DD - R;
	  di = sp->SS;
	  dj = sp->FF;
	  ORDER(d, di, dj, c, ci, cj)
          c = 0;
	  DIAG(l_ptr->m1, j, c, p+va[B[j]])		/* diagonal */
	  if ( c <= 0 ) { c = 0; ci = l_ptr->m1; cj = j; }
	  else { ci = pi; cj = pj; }
	  ORDER1(c, ci, cj, d, di, dj)
	  ORDER1(c, ci, cj, f, fi, fj)
	  sp->SS = di;
	  sp->FF = dj;
	  p = sp->CC;
	  sp->CC = c;
	  pi = sp->RR;
	  pj = sp->EE;
	  sp->RR = ci;
	  sp->EE = cj;
	  sp->DD = d;
	  if ( c > mini_score ) *flag_p = 1;
	  if ( ! rflag && ( (ci > l_ptr->rl && cj > l_ptr->cl) || (di > l_ptr->rl && dj > l_ptr->cl) || (fi > l_ptr->rl && fj > l_ptr->cl) ) ) rflag = 1;
	}

	l_ptr->c_ss[l_ptr->m1].HH = l_ptr->CCC[l_ptr->n1].CC;
	l_ptr->c_ss[l_ptr->m1].II = l_ptr->CCC[l_ptr->n1].RR;
	l_ptr->c_ss[l_ptr->m1].JJ = l_ptr->CCC[l_ptr->n1].EE;
	l_ptr->c_ss[l_ptr->m1].WW = f;
	l_ptr->c_ss[l_ptr->m1].XX = fi;
	l_ptr->c_ss[l_ptr->m1].YY = fj;

	if ( ! cflag && ( (ci > l_ptr->rl && cj > l_ptr->cl) || (di > l_ptr->rl && dj > l_ptr->cl)  || (fi > l_ptr->rl && fj > l_ptr->cl ) )) cflag = 1;
      }

      if ( nseq == 1 && l_ptr->n1 == (l_ptr->m1 + 1) && ! rflag ) cflag = 0;
      if ( cflag && l_ptr->n1 > 1 ) {	/* Compute one column */
	cflag = 0;
	l_ptr->n1--;
	c = 0;
	f = - (Q);
	cj = fj = l_ptr->n1;
	va = pam2[B[l_ptr->n1]];
	if ( nseq == 2 || l_ptr->mm < l_ptr->n1 ) {
	  p = 0;
	  ci = fi = pi = l_ptr->mm + 1;
	  pj = l_ptr->n1 + 1;
	  limit = l_ptr->mm;
	}
	else {
	  p = l_ptr->c_ss[l_ptr->n1].HH;
	  pi = l_ptr->c_ss[l_ptr->n1].II;
	  pj = l_ptr->c_ss[l_ptr->n1].JJ;
	  ci = fi = l_ptr->n1;
	  limit = l_ptr->n1 - 1;
	}

	for ( i = limit; i >= l_ptr->m1 ; i-- ) {
	  f = f - R;
	  c = c - qr;
	  ORDER(f, fi, fj, c, ci, cj)
	  c = l_ptr->c_ss[i].HH - qr; 
	  ci = l_ptr->c_ss[i].II;
	  cj = l_ptr->c_ss[i].JJ;
	  d = l_ptr->c_ss[i].WW - R;
	  di = l_ptr->c_ss[i].XX;
	  dj = l_ptr->c_ss[i].YY;
	  ORDER(d, di, dj, c, ci, cj)
	  c = 0;
	  DIAG(i, l_ptr->n1, c, p+va[A[i]])
	  if ( c <= 0 ) { c = 0; ci = i; cj = l_ptr->n1; }
	  else { ci = pi; cj = pj; }
	  ORDER1(c, ci, cj, d, di, dj)
	  ORDER1(c, ci, cj, f, fi, fj)
	  p = l_ptr->c_ss[i].HH;
	  l_ptr->c_ss[i].HH = c;
	  pi = l_ptr->c_ss[i].II;
	  pj = l_ptr->c_ss[i].JJ;
	  l_ptr->c_ss[i].II = ci;
	  l_ptr->c_ss[i].JJ = cj;
	  l_ptr->c_ss[i].WW = d;
	  l_ptr->c_ss[i].XX = di;
	  l_ptr->c_ss[i].YY = dj;
	  if ( c > mini_score ) *flag_p = 1;
	  if ( ! cflag && ( (ci > l_ptr->rl && cj > l_ptr->cl) || (di > l_ptr->rl && dj > l_ptr->cl)
			    || (fi > l_ptr->rl && fj > l_ptr->cl) ) )  cflag = 1;
	}

	l_ptr->CCC[l_ptr->n1].CC = l_ptr->c_ss[l_ptr->m1].HH;
	l_ptr->CCC[l_ptr->n1].RR = l_ptr->c_ss[l_ptr->m1].II;
	l_ptr->CCC[l_ptr->n1].EE = l_ptr->c_ss[l_ptr->m1].JJ;
	l_ptr->CCC[l_ptr->n1].DD = f;
	l_ptr->CCC[l_ptr->n1].SS = fi;
	l_ptr->CCC[l_ptr->n1].FF = fj;
	if ( ! rflag && ( (ci > l_ptr->rl && cj > l_ptr->cl) || (di > l_ptr->rl && dj > l_ptr->cl)
			  || (fi > l_ptr->rl && fj > l_ptr->cl )) ) rflag = 1;
      }
    }
    if (( l_ptr->m1 == 1 && l_ptr->n1 == 1) || no_cross(flag_p, v_ptr->LIST, l_ptr) ) break;
  }
  l_ptr->m1--;
  l_ptr->n1--;
}

/* recompute the area on forward pass */
static void
small_pass(const unsigned char *A,
	   const unsigned char *B,
	   int mini_score,
	   int **pam2, int Q, int R,
	   int nseq,
	   struct vert_str *v_ptr,
	   struct l_struct *l_ptr) {

  int  i, j;		/* row and column indices */
  int  c;		/* best score at current point */
  int  f;		/* best score ending with insertion */
  int  d;		/* best score ending with deletion */
  int  p;		/* best score at (i-1, j-1) */
  int  ci, cj;		/* end-point associated with c */ 
  int  fi, fj;		/* end-point associated with f */
  int  pi, pj;		/* end-point associated with p */
  space_ptr sp;
  pair_ptr z;
  int q, r, qr;
  int  *va;		/* pointer to pam2(A[i], B[j]) */
  
  int  limit;		/* lower bound on j */

  q = Q; r = R; qr = q + r;

  for ( sp = &l_ptr->CCC[l_ptr->n1 + 1], j = l_ptr->n1+1; sp <= &l_ptr->CCC[l_ptr->nn] ; sp++, j++ ) {
    sp->CC = 0;
    sp->RR = l_ptr->m1;
    sp->EE = j;
    sp->DD = - (qr);
    sp->SS = l_ptr->m1+1;
    sp->FF = j;
  }

  for ( i = l_ptr->m1 + 1; i <= l_ptr->mm; i++) {
    c = 0;				/* Initialize column 0 */
    f = - (qr);
    ci = fi = i;
    va = pam2[A[i]];
    if ( nseq == 2 || i <= l_ptr->n1 ) {
      p = 0;
      pi = i - 1;
      cj = fj = pj = l_ptr->n1;
      limit = l_ptr->n1 + 1;
    }
    else {
      p = l_ptr->CCC[i].CC;
      pi = l_ptr->CCC[i].RR;
      pj = l_ptr->CCC[i].EE;
      cj = fj = i;
      limit = i + 1;
    }

    for ( j = limit, sp = &l_ptr->CCC[j] ; j <= l_ptr->nn ; j++, sp++ )  {  
      d = sp->DD;
      c = -1;
      DIAG(i, j, c, p+va[B[j]])		/* diagonal */
      if (c < 0) {
	p = sp->CC; pi = sp->RR; pj = sp->EE;
	if (f >= 0) {
	  c = f; ci = fi; cj = fj;
	  ORDER1(c, ci, cj,  d, sp->SS, sp->FF)
	  sp->CC = c; sp->RR = ci; sp->EE = cj;
	  sp->DD -= r; f-=r;
	} 
	else if (d >= 0) {
	  sp->CC = d; sp->RR = sp->SS; sp->EE = sp->FF;  
	  sp->DD -= r;
	}
	else {
	  sp->CC = 0;
	  sp->RR=i;
	  sp->EE = j;
	}
      }
      else { 
	ci = pi; cj = pj;
	ORDER1(c, ci, cj,  f, fi, fj)
	  ORDER1(c, ci, cj,  d, sp->SS, sp->FF)
	  p = sp->CC;
	sp->CC = c;
	pi = sp->RR;
	sp->RR = ci;
	pj = sp->EE;
	sp->EE = cj;
	f-=r;
	if (c >= qr) {
	  if ( c > mini_score ) /* add the score into list */
	    addnode(c, ci, cj, i, j, v_ptr);
	  d -= r; c-=qr;
	  ORDER1(f, fi, fj, c, ci, cj)
          ORDER1(d, sp->SS, sp->FF, c, ci, cj)
	  sp->DD = d;
	}
	else {
	  sp->DD -= r;
	}
      }
    }
  }
}

/* Add a new node into list.  */

static void
addnode(int c, int ci, int cj, int i, int j, struct vert_str *v_ptr) {

  bool found;	/* 1 if the node is in LIST */

  found = 0;
  if ( v_ptr->most != NULL && v_ptr->most->STARI == ci && v_ptr->most->STARJ == cj)
    found = 1;
  else {
    for ( v_ptr->most = v_ptr->LIST; v_ptr->most; v_ptr->most = v_ptr->most->next ) {
      if ( v_ptr->most->STARI == ci && v_ptr->most->STARJ == cj) {
	found = 1;
	break;
      }
    }
  }
  if ( found ) {
    if ( v_ptr->most->SCORE < c ) {
      v_ptr->most->SCORE = c;
      v_ptr->most->ENDI = i;
      v_ptr->most->ENDJ = j;
    }
    if ( v_ptr->most->TOP > i ) v_ptr->most->TOP = i;
    if ( v_ptr->most->BOT < i ) v_ptr->most->BOT = i;
    if ( v_ptr->most->LEFT > j ) v_ptr->most->LEFT = j;
    if ( v_ptr->most->RIGHT < j ) v_ptr->most->RIGHT = j;
  }
  else { 
    v_ptr->numnode++;
    v_ptr->most = (vertex_p) ckalloc(1,sizeof(vertex));
    v_ptr->most->SCORE = c;
    v_ptr->most->STARI = ci;
    v_ptr->most->STARJ = cj;
    v_ptr->most->ENDI = i;
    v_ptr->most->ENDJ = j;
    v_ptr->most->TOP = v_ptr->most->BOT = i;
    v_ptr->most->LEFT = v_ptr->most->RIGHT = j;
    v_ptr->most->next = v_ptr->LIST;
    v_ptr->LIST = v_ptr->most;
  }
}

/* Find and remove the largest score in list */

static vertex_p
findmax(struct vert_str *v_ptr) {
  vertex_p  ap, cur;
  register int score;

  for ( score = (v_ptr->LIST)->SCORE, cur = NULL, ap = (v_ptr->LIST); ap->next; ap = ap->next) {
    if ( ap->next->SCORE > score ) {
      cur = ap; score = ap->next->SCORE;
    }
  }
  if (cur) {ap = cur->next; cur->next = ap->next; }
  else { ap = v_ptr->LIST; v_ptr->LIST = (v_ptr->LIST)->next;}
  v_ptr->numnode--;
  v_ptr->most = v_ptr->LIST;
  return ( ap );
}

/* return 1 if no node in LIST share vertices with the area */

static bool
no_cross(int *flag_p, vertex_p LIST, struct l_struct *l_ptr) {

  vertex_p  cur;

  for ( cur = LIST; cur; cur = cur->next ) { 
    if ( cur->STARI <= l_ptr->mm && cur->STARJ <= l_ptr->nn && cur->BOT >= l_ptr->m1-1 && 
	 cur->RIGHT >= l_ptr->n1-1 && (cur->STARI < l_ptr->rl || cur->STARJ < l_ptr->cl)) { 
      if ( cur->STARI < l_ptr->rl ) l_ptr->rl = cur->STARI;
      if ( cur->STARJ < l_ptr->cl ) l_ptr->cl = cur->STARJ;
      *flag_p = 1;
      break;
    }
  }
  return !cur;
}

/* The following definitions are for function diff() */

#define gap(k)  ((k) <= 0 ? 0 : Q+R*(k)) /* k-symbol indel score */

/* Append "Delete k" op */
#define DEL(k)				\
{ l_ptr->I += k;			\
  if (*last < 0)			\
    *last = (*sapp)[-1] -= (k);		\
  else {				\
    *last = (*sapp)[0] = -(k);		\
    (*sapp)++;				\
  }					\
}

/* Append "Insert k" op */
#define INS(k)				\
{ l_ptr->J += k;			\
  if (*last < 0) {			\
   (*sapp)[-1] = (k);			\
   (*sapp)[0] = *last;			\
   (*sapp)++;				\
  }					\
  else {				\
    *last = (*sapp)[0] = (k);		\
    (*sapp)++;				\
  }					\
}

/* diff(A,B,M,N,tb,te) returns the score of an optimum conversion between
   A[1..M] and B[1..N] that begins(ends) with a delete if tb(te) is zero
   and appends such a conversion to the current script.                   */

static int
diff(const unsigned char *A,
     const unsigned char *B,
     int M, int N,
     int tb, int te,
     int two_seq,
     int **pam2, int Q, int R,
     int **sapp, int *last,
     struct l_struct *l_ptr) {

  int   midi, midj, type;	/* Midpoint, type, and cost */
  int limit;
  int midc;

  register int i, j;
  register int c, e, d, s;

  pair_ptr z;

  int t;
  int *va;
  int qr;

  bool tt;

  qr = Q + R;

  /* Boundary cases: M <= 1 or N == 0 */

  if (N <= 0){
    if (M > 0) {DEL(M)}
    return - gap(M);
  }

  if (M <= 1) {
    if (M <= 0) {
      INS(N)
      return - gap(N);
    }

    if (tb > te) tb = te;
    midc = - (tb + R + gap(N) );
    midj = 0;

    va = pam2[A[1]];
    j = 2 + l_ptr->I - l_ptr->J;
    if (two_seq || j < 1) j = 1;
    for ( ; j <= N; j++) {
      for ( tt = 1, z = l_ptr->row[l_ptr->I+1]; z != NULL; z = z->NEXT ) {
	if ( z->COL == j+l_ptr->J ) { tt = 0; break; }		
      }
      if (tt) {
	c = va[B[j]] - ( gap(j-1) + gap(N-j) );
	if (c > midc) { midc = c;  midj = j; }
      }
    }

    if (midj == 0) { INS(N) DEL(1) }
    else  {
      if (midj > 1) INS(midj-1)
      *last = (*sapp)[0] = 0;
      (*sapp)++;

      /* mark (A[I],B[J]) as used: put J into list row[I] */	
      l_ptr->I++; l_ptr->J++;
      z = ( pair_ptr ) ckalloc(1,sizeof(pair));
      z->COL = l_ptr->J;			
      z->NEXT = l_ptr->row[l_ptr->I];				
      l_ptr->row[l_ptr->I] = z;
      if (midj < N) INS(N-midj)
    }
    return midc;
  }

  /* Divide: Find optimum midpoint (midi,midj) of cost midc */

  midi = M/2;		/* Forward phase:                          */
  l_ptr->r_ss[0].CC = 0;		/*   Compute C(M/2,k) & D(M/2,k) for all k */
  t = -Q;
  for (j = 1; j <= N; j++) {
    l_ptr->r_ss[j].CC = t = t-R;
    l_ptr->r_ss[j].DD = t-Q;
  }
  t = -tb;
  for (i = 1; i <= midi; i++) {
    va = pam2[A[i]];
    t = t-R;
    j = i + l_ptr->I - l_ptr->J;
    if (two_seq || j <= 0) {
      j = 0;
      s = l_ptr->r_ss[0].CC;
      l_ptr->r_ss[0].CC = c = t;
    }
    else {
      if ( (c = (s = l_ptr->r_ss[j].CC) - qr) < (d = l_ptr->r_ss[j].DD)) c = d;
      l_ptr->r_ss[j].CC = l_ptr->r_ss[j].DD = c;
    }
    e = c-Q;
    for (j++ ; j <= N; j++) {
      if ((c = c - qr) > (e = e - R)) e = c;
      if ((c = l_ptr->r_ss[j].CC - qr) > (d = l_ptr->r_ss[j].DD - R)) d = c;
      DIAG(i+l_ptr->I, j+l_ptr->J, c, s+va[B[j]])
      if (c < d) c = d;
      if (c < e) c = e;
      s = l_ptr->r_ss[j].CC;
      l_ptr->r_ss[j].CC = c;
      l_ptr->r_ss[j].DD = d;
    }
  }
  l_ptr->r_ss[0].DD = l_ptr->r_ss[0].CC;

  l_ptr->r_ss[N].RR = 0;			/* Reverse phase:                          */
  t = -Q;			/*   Compute R(M/2,k) & S(M/2,k) for all k */
  for (j = N-1; j >= 0; j--) {
    l_ptr->r_ss[j].RR = t = t-R;
    l_ptr->r_ss[j].SS = t-Q;
  }
  t = -te;

  for (i = M-1; i >= midi; i--) {
    s = l_ptr->r_ss[N].RR;
    l_ptr->r_ss[N].RR = c = t = t-R;
    e = t-Q;
    va = pam2[A[i+1]];
    limit = i + l_ptr->I - l_ptr->J + 1;
    if (two_seq || limit < 0) limit = 0;
    for (j = N-1; j >= limit; j--) {
      if ((c = c - qr) > (e = e - R)) e = c;
      if ((c = l_ptr->r_ss[j].RR - qr) > (d = l_ptr->r_ss[j].SS - R)) d = c;
      DIAG(i+1+l_ptr->I, j+1+l_ptr->J, c, s+va[B[j+1]])
      if (c < d) c = d;
      if (c < e) c = e;
      s = l_ptr->r_ss[j].RR;
      l_ptr->r_ss[j].RR = c;
      l_ptr->r_ss[j].SS = d;
    }
  }
  l_ptr->r_ss[N].SS = l_ptr->r_ss[N].RR;

  midc = l_ptr->r_ss[0].CC+l_ptr->r_ss[0].RR;		/* Find optimal midpoint */
  midj = 0;
  type = 1;
  limit = midi + l_ptr->I - l_ptr->J + 1;
  if (two_seq || limit < 0) limit = 0;
  for (j = limit; j <= N; j++) {
    if ((c = l_ptr->r_ss[j].CC + l_ptr->r_ss[j].RR) >= midc) {
      if (c > midc ||
	  (l_ptr->r_ss[j].CC != l_ptr->r_ss[j].DD && l_ptr->r_ss[j].RR == l_ptr->r_ss[j].SS)) {
	midc = c;
	midj = j;
      }
    }
  }
  for (j = N; j >= limit; j--)
    if ((c = l_ptr->r_ss[j].DD + l_ptr->r_ss[j].SS + Q) > midc)
      { midc = c;
      midj = j;
      type = 2;
      }

/* Conquer: recursively around midpoint */

  if (type == 1) {
    (void)diff(A,B,midi,midj,tb,Q,two_seq,pam2,Q,R, sapp, last, l_ptr);
    (void)diff(A+midi,B+midj,M-midi,N-midj,Q,te,two_seq,pam2,Q,R, sapp, last, l_ptr);
  }
  else {
    (void)diff(A,B,midi-1,midj,tb,0,two_seq,pam2,Q,R, sapp, last, l_ptr);
    DEL(2);
    (void)diff(A+midi+1,B+midj,M-midi-1,N-midj,0,te,two_seq,pam2,Q,R, sapp, last, l_ptr);
  }
  return midc;
}

/* CHECK_SCORE - return the score of the alignment stored in S */

static int CHECK_SCORE(const unsigned char *A, const unsigned char *B,
		       int M, int N,
		       int *S, int **w,
		       int qq, int rr, int *NC)
{ 
  register int   i,  j, op, nc;
  int itmp, score;
#ifdef SHOW_ALIGN_SCORE
  int mx_l_score;
#endif

  /*  print_seq_prof(A,M,w,iw); */

  score = i = j = op = nc = 0;
#ifdef SHOW_ALIGN_SCORE
  mx_l_score = 0;
  printf("#===start\n");
  printf("#i j pam2 score mx_l_score\n");
#endif
  while (i < M || j < N) {
    op = *S++;
    if (op == 0) {
      itmp = w[A[++i]][B[++j]];
      score += itmp;
      nc++;
    }
    else if (op > 0) {
      itmp = -(qq + op*rr);
      score += itmp;
      j += op;
      nc += op;
    } else {	/* op < 0 */
      itmp = - (qq - op*rr);
      score += itmp;
      i -= op;	/* i increased */
      nc -= op;	/* nc increased */
    }
#ifdef SHOW_ALIGN_SCORE
    if (score > mx_l_score) mx_l_score = score;
    printf("%d\t%d\t%d\t%d\t%d\n",i, j, itmp, score, mx_l_score);
#endif
  }
#ifdef SHOW_ALIGN_SCORE
  printf("%d\t%d\tend\t%d\t%d\n====\n",i, j, score, mx_l_score);
#endif
  *NC = nc;
  return score;
}

/* ckalloc - allocate space; check for success */
void *ckalloc(size_t amount, size_t size)
{
  void *p;
  static size_t mtotal;

  mtotal += amount * size;

  if ((p = malloc( amount * size )) == NULL) {
    fprintf(stderr,"Ran out of near memory: %ld*%ld/%ld\n",amount,size,mtotal);
    exit(1);
  }
  return(p);
}
