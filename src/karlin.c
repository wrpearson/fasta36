/**************** Statistical Significance Parameter Subroutine ****************

     $Id: karlin.c 625 2011-03-23 17:21:38Z wrp $
     $Revision: 625 $

    Version 1.0	February 2, 1990
    Version 2.0	March 18,   1993

    Program by:	Stephen Altschul

    Address:	National Center for Biotechnology Information
    National Library of Medicine
    National Institutes of Health
    Bethesda, MD  20894

    Internet:	altschul@ncbi.nlm.nih.gov

    See:	Karlin, S. & Altschul, S.F. "Methods for Assessing the Statistical
    Significance of Molecular Sequence Features by Using General Scoring
    Schemes,"  Proc. Natl. Acad. Sci. USA 87 (1990), 2264-2268.

    Computes the parameters lambda and K for use in calculating the
    statistical significance of high-scoring segments or subalignments.

    The scoring scheme must be integer valued.  A positive score must be
    possible, but the expected (mean) score must be negative.

    A program that calls this routine must provide the value of the lowest
    possible score, the value of the greatest possible score, and a pointer
    to an array of probabilities for the occurence of all scores between
    these two extreme scores.  For example, if score -2 occurs with
    probability 0.7, score 0 occurs with probability 0.1, and score 3
    occurs with probability 0.2, then the subroutine must be called with
    low = -2, high = 3, and pr pointing to the array of values
    { 0.7, 0.0, 0.1, 0.0, 0.0, 0.2 }.  The calling program must also provide
    pointers to lambda and K; the subroutine will then calculate the values
    of these two parameters.  In this example, lambda=0.330 and K=0.154.

    The parameters lambda and K can be used as follows.  Suppose we are
    given a length N random sequence of independent letters.  Associated
    with each letter is a score, and the probabilities of the letters
    determine the probability for each score.  Let S be the aggregate score
    of the highest scoring contiguous segment of this sequence.  Then if N
    is sufficiently large (greater than 100), the following bound on the
    probability that S is greater than or equal to x applies:
	
    P( S >= x )   <=   1 - exp [ - KN exp ( - lambda * x ) ].
	
    In other words, the p-value for this segment can be written as
    1-exp[-KN*exp(-lambda*S)].

    This formula can be applied to pairwise sequence comparison by assigning
    scores to pairs of letters (e.g. amino acids), and by replacing N in the
    formula with N*M, where N and M are the lengths of the two sequences
    being compared.

    In addition, letting y = KN*exp(-lambda*S), the p-value for finding m
    distinct segments all with score >= S is given by:

    2             m-1           -y
    1 - [ 1 + y + y /2! + ... + y   /(m-1)! ] e

    Notice that for m=1 this formula reduces to 1-exp(-y), which is the same
    as the previous formula.

*******************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define MAXIT 25  /* Maximum number of iterations used in calculating lambda */
#define NMAP_X 23
#define NMAP 33

#define TINY 1e-6

/* first build a residue map to automatically put residues in score bins */

#include "defs.h"
#include "param.h"

/* initialize the Karlin frequency, probability arrays using
   a specific query sequence */

int karlin(int , int, double *, double *, double *);
static int karlin_k(int , int , double *, double *, double *, double *);

void init_karlin(const unsigned char *aa0, int n0, struct pstruct *ppst,
		 double *aa0_f, double **kp)
{
  int kar_nsq, kar_range, kar_min, kar_max;

  const unsigned char *aa0p;
  int i;
  int r_cnt[NMAP+1];
  double fn0, *kar_p;
  
  kar_range = ppst->pam_h - ppst->pam_l + 1;
  if (*kp == NULL) {
    if ((kar_p=(double *)calloc(kar_range+1,sizeof(double)))==NULL) {
      fprintf(stderr," cannot allocate kar_p array: %d\n",kar_range+1);
      exit(1);
    }
    *kp = kar_p;
  }
  kar_nsq = ppst->nsq;	/* alphabet size */
  kar_min = ppst->pam_l;	/* low pam value */
  kar_max = ppst->pam_h;	/* high pam value */

  /* must have at least 1 residue of each type */
  r_cnt[NMAP]=0;
  for (i=1; i<=kar_nsq; i++) r_cnt[i]=1; 
 
  fn0 = 100.0/(double)(n0+kar_nsq);	/* weight of each residue */

  aa0p = aa0;
  /* increment residue count for each residue in query sequence */
  while (*aa0p) r_cnt[ppst->hsqx[*aa0p++]]++;

  /* map all unmapped residues to 'X' */
  r_cnt[NMAP_X] += r_cnt[NMAP];
 
  for (i=1; i<=kar_nsq; i++) aa0_f[i] = fn0*(double)r_cnt[i];
}

double nt_f[] = {0.0, 0.25, 0.25, 0.25, 0.25 };

/* Robinson and Robinson frequencies */
double aa_f[] = {
/* NULL */ 0.00,
/* A */   0.0780474700897585,
/* R */   0.0512953149316987,
/* N */   0.0448725775979007,
/* D */   0.0536397361638076,
/* C */   0.0192460110427568,
/* Q */   0.0426436013507063,
/* E */   0.0629485981204668,
/* G */   0.0737715654561964,
/* H */   0.0219922696262025,
/* I */   0.0514196403000682,
/* L */   0.090191394464413,
/* K */   0.0574383201866657,
/* M */   0.0224251883196316,
/* F */   0.0385564048655621,
/* P */   0.0520279465667327,
/* S */   0.0711984743501224,
/* T */   0.0584129422708473,
/* W */   0.013298374223799,
/* Y */   0.0321647488738564,
/* V */   0.0644094211988074};

/* initialize the Karlin frequency, probability arrays using
   an "average" composition (average length if n0 <=0) */

void
init_karlin_a(struct pstruct *ppst, double *aa0_f, double **kp)
{
  int kar_nsq, kar_range;

  int i;
  double fn0, *kar_p;

  kar_range = ppst->pam_h - ppst->pam_l + 1;
  if (*kp == NULL) {
    if ((kar_p=(double *)calloc(kar_range+1,sizeof(double)))==NULL) {
      fprintf(stderr," cannot allocate kar_p array: %d\n",kar_range+1);
      exit(1);
    }
  *kp = kar_p;
  }

  if (ppst->nt_align) {
    kar_nsq = 4;
    for (i=1; i<=kar_nsq; i++) aa0_f[i] = nt_f[i];
  }
  else if (ppst->dnaseq==SEQT_PROT || ppst->dnaseq == SEQT_UNK) {
    kar_nsq = 20;
    for (i=1; i<=kar_nsq; i++) aa0_f[i] = aa_f[i];
  }
  else {
    kar_nsq = ppst->nsq;
    fn0 = 1.0/(double)(kar_nsq-1);
    for (i=1; i< kar_nsq; i++) aa0_f[i] = fn0;
    aa0_f[kar_nsq]=0.0;
  }

}

/* calculate set up karlin() to calculate Lambda, K, by calculating
   aa1 frequencies */
int
do_karlin(const unsigned char *aa1, int n1,
	  int **pam2, const struct pstruct *ppst,
	  double *aa0_f, double *kar_p, double *lambda, double *H)
{
  register unsigned const char *aap;
  int kar_range, kar_min, kar_max, kar_nsq;
  int r_cnt[NMAP+1];
  double aa1_f[NMAP];
  double fn1, kar_tot;
  int i, j;

  kar_nsq = ppst->nsq;
  kar_min = ppst->pam_l;
  kar_max = ppst->pam_h;
  kar_range = kar_max - kar_min + 1;

  r_cnt[NMAP]=0;
  for (i=1; i<=kar_nsq; i++) r_cnt[i]=1;

  /* residue counts */

  aap=aa1;
  while (*aap) r_cnt[ppst->hsqx[*aap++]]++;

  r_cnt[NMAP_X] += r_cnt[NMAP];
  
  /* residue frequencies */
  fn1 = 100.0/(double)(n1+kar_nsq);
  for (i=1; i<=kar_nsq; i++) aa1_f[i]= fn1*(double)r_cnt[i];

  for (i=0; i<=kar_range; i++) kar_p[i] = 0.0;

  for (i=1; i<=kar_nsq; i++) {
    for (j=1; j<=kar_nsq; j++)
      kar_p[pam2[i][j]-kar_min] += aa0_f[i]*aa1_f[j];
  }

  kar_tot = 0.0;
  for (i=0; i<=kar_range; i++) kar_tot += kar_p[i];
  if (kar_tot <= 0.00001) return 0;

  for (i=0; i<=kar_range; i++) kar_p[i] /= kar_tot;

  return karlin(kar_min, kar_max, kar_p, lambda, H);
}

int
do_karlin_a(int **pam2, struct pstruct *ppst,
	  double *aa0_f, double *kar_p, double *lambda, double *K, double *H)
{
  double *aa1fp;
  int kar_range, kar_min, kar_max, kar_nsq;
  double aa1_f[NMAP];
  double fn1, kar_tot;
  int i, j;

  kar_min = ppst->pam_l;
  kar_max = ppst->pam_h;
  kar_range = kar_max - kar_min + 1;

  kar_tot = 0.0;
  if (ppst->nt_align ) {
    kar_nsq = 4;
    aa1fp = nt_f;
    for (i=1; i<=kar_nsq; i++) {kar_tot += aa1fp[i];}
    for (i=1; i<=kar_nsq; i++) {aa1_f[i]= aa1fp[i]/kar_tot;}
  }
  else if (!ppst->nt_align) {
    kar_nsq = 20;
    aa1fp = aa_f;
    for (i=1; i<=kar_nsq; i++) {kar_tot += aa1fp[i];}
    for (i=1; i<=kar_nsq; i++) {aa1_f[i]= aa1fp[i]/kar_tot;}
  }
  else {
    kar_nsq = ppst->nsq;
    fn1 = 1.0/(double)(kar_nsq-1);
    for (i=1; i< kar_nsq; i++) aa1_f[i] = fn1;
    aa1_f[kar_nsq]=0.0;
  }

  for (i=0; i<=kar_range; i++) kar_p[i] = 0.0;

  for (i=1; i<=kar_nsq; i++) {
    for (j=1; j<kar_nsq; j++)
      kar_p[pam2[i][j]-kar_min] += aa0_f[i]*aa1_f[j];
  }    

  kar_tot = 0.0;
  for (i=0; i<=kar_range; i++) kar_tot += kar_p[i];
  if (kar_tot <= 0.00001) return 0;

  for (i=0; i<=kar_range; i++) kar_p[i] /= kar_tot;

  return karlin_k(kar_min, kar_max, kar_p, lambda, K, H);
}

/* take a array of letters and pam information and get *lambda, *H */
int
karlin(int low,			/* Lowest score (must be negative)    */
       int high,		/* Highest score (must be positive)   */
       double *pr,		/* Probabilities for various scores   */
       double *lambda_p,	/* Pointer to parameter lambda        */
       double *H_p)		/* Pointer to parameter H              */
{
  int i,range, nit;
  double up,new,sum,av,beta,ftemp;
  double lambda;
  double *ptr1;

  /* Calculate the parameter lambda */

  range = high-low;

  /* check for E() < 0.0 */
  sum = 0;
  ptr1 = pr;
  for (i=low; i <= high ; i++) sum += i* (*ptr1++);
  if (sum >= 0.0) {
#ifdef DEBUG
    fprintf(stderr," (karlin lambda) non-negative expected score: %.4lg\n",
	    sum);
#endif
    return 0;
  }

  /* up is upper bound on lambda */
  up=0.5;
  do {
    up *= 2.0;
    ptr1=pr;

    beta=exp(up);

    ftemp=exp(up*(low-1));
    sum = 0.0;
    for (i=0; i<=range; ++i) sum+= *ptr1++ * (ftemp*=beta);
  }
  while (sum<1.0);

  /* avoid overflow from very large lambda*S */
/*
  do {
    up /= 2.0;
    ptr1=pr;
    beta=exp(up);

    ftemp=exp(up*(low-1));
    sum = 0.0;
    for (i=0; i<=range; ++i) sum+= *ptr1++ * (ftemp*=beta);
  } while (sum > 2.0);

  up *= 2.0;
*/	/* we moved past, now back up */

  /*	for (lambda=j=0;j<25;++j) { */
  lambda = 0.0;
  nit = 0;
  while ( nit++ < MAXIT ) {
    new = (lambda+up)/2.0;
    beta = exp(new);
    ftemp = exp(new*(low-1));
    ptr1=pr;
    sum = 0.0;
    /* multiply by exp(new) for each score */
    for (i=0;i<=range;++i) sum+= *ptr1++ * (ftemp*=beta);

    if (sum > 1.0 + TINY) up=new;
    else {
      if ( fabs(lambda - new) < TINY ) goto done;
      lambda = new;
    }
  }

  if (lambda <= 1e-10) {
    lambda = -1.0;
    return 0;
  }

 done:
  *lambda_p = lambda;

  /* Calculate the parameter K */

  ptr1=pr;
  ftemp=exp(lambda*(low-1));
  for (av=0.0, i=low; i<=high; ++i)
    av+= *ptr1++ *i*(ftemp*=beta);
  *H_p= lambda*av;

  return 1;		/* Parameters calculated successfully */
}

static int a_gcd (int, int);

/* take a array of letters and pam information and get *lambda, *K, *H */
static int
karlin_k(int low,		/* Lowest score (must be negative)    */
	 int high,		/* Highest score (must be positive)   */
	 double *pr,	/* Probabilities for various scores   */
	 double *lambda_p,	/* Pointer to parameter lambda        */
	 double *K_p,
	 double *H_p)		/* Pointer to parameter H              */
{
  int i,j,range,lo,hi,first,last, nit;
  double up,new,sum,Sum,av,beta,oldsum,ratio,ftemp;
  double lambda;
  double *P,*ptrP,*ptr2;
  double *ptr1;

  /* Calculate the parameter lambda */

  range = high-low;

  /* check for E() < 0.0 */
  sum = 0;
  ptr1 = pr;
  for (i=low; i <= high ; i++) sum += i* (*ptr1++);
  if (sum >= 0.0) {
    fprintf(stderr," (karlin lambda) non-negative expected score: %.4lg\n",
	    sum);
    /* perhaps we should return values for BLOSUM50 here to avoid fp
       underflow later */
    return 0;
  }

  /* up is upper bound on lambda */
  up=0.5;
  do {
    up *= 2.0;
    ptr1=pr;

    beta=exp(up);

    ftemp=exp(up*(low-1));
    sum = 0.0;
    for (i=0; i<=range; ++i) sum+= *ptr1++ * (ftemp*=beta);
  }
  while (sum<1.0);

  /* avoid overflow from very large lambda*S */
  /*
  do {
    up /= 2.0;
    ptr1=pr;
    beta=exp(up);

    ftemp=exp(up*(low-1));
    sum = 0.0;
    for (i=0; i<=range; ++i) sum+= *ptr1++ * (ftemp*=beta);
  } while (sum > 2.0);

  up *= 2.0;
  */
  /* we moved past, now back up */

  /*	for (lambda=j=0;j<25;++j) { */
  lambda = 0.0;
  nit = 0;
  while ( nit++ < MAXIT ) {
    new = (lambda+up)/2.0;
    beta = exp(new);
    ftemp = exp(new*(low-1));
    ptr1=pr;
    sum = 0.0;
    /* multiply by exp(new) for each score */
    for (i=0;i<=range;++i) sum+= *ptr1++ * (ftemp*=beta);

    if (sum > 1.0 + TINY) up=new;
    else {
      if ( fabs(lambda - new) < TINY ) goto done;
      lambda = new;
    }
  }

  if (lambda <= 1e-10) {
    lambda = -1.0;
    return 0;
  }

 done:
  *lambda_p = lambda;

  /* Calculate the parameter H */

  ptr1=pr;
  ftemp=exp(lambda*(low-1));
  for (av=0.0, i=low; i<=high; ++i) av+= *ptr1++ *i*(ftemp*=beta);
  *H_p= lambda*av;

  /* Calculate the pamameter K */
  Sum=lo=hi=0;
  P= (double *) calloc(MAXIT*range+1,sizeof(double));
  for (*P=sum=oldsum=j=1;j<=MAXIT && sum>0.001;Sum+=sum/=j++) {
    first=last=range;
    for (ptrP=P+(hi+=high)-(lo+=low); ptrP>=P; *ptrP-- =sum) {
      ptr1=ptrP - first;
      ptr2=pr + first;
      for (sum=0,i=first; i<=last; ++i) sum += *ptr1-- * *ptr2++;
      if (first) --first;
      if (ptrP-P<=range) --last;
    }
    ftemp=exp(lambda*(lo-1));
    for (sum=0,i=lo;i;++i) sum+= *++ptrP * (ftemp*=beta);
    for (;i<=hi;++i) sum+= *++ptrP;
    ratio=sum/oldsum;
    oldsum=sum;
  }
  for (;j<=200;Sum+=oldsum/j++) oldsum*=ratio;
  for (i=low; !pr[i-low]; ++i);
  for (j= -i;i<high && j>1;) if (pr[++i-low]) j=a_gcd(j,i);
  *K_p = (j*exp(-2*Sum))/(av*(1.0-exp(- lambda*j)));
  free(P);

  return 1;		/* Parameters calculated successfully */
}

int
a_gcd(int a, int b)
{
  int c;

  if (b<0) b= -b;
  if (b>a) { c=a; a=b; b=c; }
  for (;b;b=c) { c=a%b; a=b; }
  return a;
}

