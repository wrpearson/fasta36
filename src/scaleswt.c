/* $Id: scaleswt.c 714 2011-05-05 00:33:40Z wrp $ */

/* copyright (c) 1995, 1996, 2000, 2002, 2014 by William R. Pearson and The
   Rectors & Visitors of the University of Virginia */

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

/* as of 24 Sept, 2000 - scaleswn uses no global variables */

/* 
  This version is designed for fasts/f, which used Tatusov
  probabilities for statistical estimates, but still needs a
  quick-and-dirty linear regression fit to rank things

  For comparisons that obey tatusov statistics, we try whenever
  possible to provide accurate e_scores, rather than raw scores.  As a
  result, no lambda/K fitting is required; and process_hist() can be
  called at the very beginning of the search to initialize some of the
  statistics structures and find_zp().

  find_z() must still return a valid z_score surrogate, as
  comp_lib.c/p2_complib.c continue to use z_score's to rank hits, save
  the best, etc.

  If e_score's cannot be calculated, the process_hist() provides
  linear regression fitting for conventional z_score estimates.

*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <limits.h>
#include <float.h>
#include <math.h>

#include <limits.h>

#include "defs.h"
#include "param.h"
#include "structs.h"
#include "best_stats.h"

#define MAXHIST 50
#define MAX_LLEN 200
#define LHISTC 5
#define VHISTC 5
#define MAX_SSCORE 300

#define LENGTH_CUTOFF 10 /* minimum database sequence length allowed, for fitting */

#define LN_FACT 10.0
#ifndef M_LN2
#define M_LN2 0.69314718055994530942
#endif

#define EULER_G 0.57721566490153286060
#define PI_SQRT6 1.28254983016186409554

#ifndef M_SQRT2
#define M_SQRT2 1.41421356237
#endif
#define LN200 5.2983173666
#define ZS_MAX 400.0	/* used to prevent underflow on some machines */
#define TOLERANCE 1.0e-12
#define TINY 1.0e-6

/* used by AVE_STATS, REG_STATS, REGI_STATS, REG2_STATS*/
struct pstat_str {
  int zsflag;
  double ngLambda, ngK, ngH;
  double rho, rho_e, mu, mu_e, mean_var, var_e;  /* ?_e:std. error of ? */
/* used by REG2_STATS */
  double rho2, mu2, var_cutoff;
  int n_trimmed; /* excluded because of high z-score */
  int n1_trimmed, nb_trimmed, nb_tot; /* excluded because of bin */
  double tat_a, tat_b, tat_c, spacefactor;
  int have_tat;
  int tie_j;
  int eval_is_pval;
  long zdb_size;
};

#define AVE_STATS 0	/* no length effect, only mean/variance */
double find_zt(int score, double escore, int len, double comp, struct pstat_str *);

double find_zn(int score, double escore, int len, double comp, struct pstat_str *);

double power(double, int);

void sortbesto(double *, int );
extern void sortbeste(struct beststr **bptr, int nbest);

int proc_hist_n(struct stat_str *sptr, int n, 
		 struct pstruct *ppst, struct hist_str *histp, int do_trim,
		 struct pstat_str *);

#define REG_STATS 1	/* length-regression scaled */
double find_zr(int score, double escore, int len, double comp, struct pstat_str *);

int proc_hist_r(struct stat_str *sptr, int n,
		 struct pstruct *ppst, struct hist_str *histp,
		 int do_trim, struct pstat_str *rs);

static double (*find_zp)(int score, double escore, int len, double comp,
			 struct pstat_str *) = &find_zr;

double find_z(int score, double escore, int len, double comp, struct pstat_str *);

/* print out all pstat_str info for independent calculation */
void
pstat_info(char *info_str, int info_str_n, char *comment, struct pstat_str *pu);

struct llen_str {
  int min, max;
  int max_score, min_score;
  int *hist;
  double *score_sums, *score2_sums;
  double *score_var;
  int max_length, min_length, zero_s;
  int fit_flag;
};

static void inithist(struct llen_str *, struct pstruct *, int);
static void free_hist( struct llen_str *);
static void addhist(struct llen_str *, int, int, int);
static void prune_hist(struct llen_str *, int, int, int, long *);
void inithistz(int, struct hist_str *histp);
void addhistz(double zs, struct hist_str *histp);

static void fit_llen(struct llen_str *, struct pstat_str *);
static void fit_llens(struct llen_str *, struct pstat_str *);

void linreg(double *lny, double *x, double *lnx, int n,
	    double *a, double *b, double *c, int start);

double calc_spacefactor(const unsigned char *, int, int, int);

double det(double a11, double a12, double a13,
	   double a21, double a22, double a23,
	   double a31, double a32, double a33);

double factorial (int a, int b);

/* void set_db_size(int, struct db_str *, struct hist_str *); */

#ifdef DEBUG
FILE *tmpf;
#endif

int
process_hist(struct stat_str *sptr, int nstats, 
	     const struct mngmsg *m_msg,
	     struct pstruct *ppst,
	     struct hist_str *histp,
	     struct pstat_str **rs_sp,
	     int do_hist
	     )
{
  int zsflag, do_trim;
  struct pstat_str *rs_s;

  if (ppst->zsflag < 0) {
    *rs_sp = NULL;
    return ppst->zsflag;
  }

  ppst->zs_off = 0.0;

  if (*rs_sp == NULL) {
    if ((rs_s=(struct pstat_str *)calloc(1,sizeof(struct pstat_str)))==NULL) {
      fprintf(stderr," cannot allocate rs_snion: %ld\n",sizeof(struct pstat_str));
      exit(1);
    }
    else *rs_sp = rs_s;
  }
  else {
    rs_s = *rs_sp;
    memset(rs_s,0,sizeof(struct pstat_str));
  }

  rs_s->zsflag = zsflag = ppst->zsflag;
  rs_s->zdb_size = ppst->zdb_size;

  if (m_msg->escore_flg) {
    find_zp = &find_zt;
    inithistz(MAXHIST,histp);
    rs_s->eval_is_pval = 1;
    return 1;
  }

  if (nstats < 20) {
    fprintf(stderr," too few sequences for sampling: %d\n",nstats);
    free(rs_s);
    *rs_sp = NULL;
    return -1;
  }

  rs_s->ngLambda = m_msg->Lambda;
  rs_s->ngK = m_msg->K;
  rs_s->ngH = m_msg->H;

  if (zsflag >= 20) {
    zsflag = ppst->zsflag2;
    do_trim = 0;
  }
  else if (zsflag >= 10) {
    zsflag -= 10;
    do_trim = 0;
  }
  else do_trim = 1;

  rs_s->eval_is_pval = 0;
  find_zp = &find_zr;

  return rs_s->zsflag = proc_hist_r(sptr, nstats, ppst, histp, do_trim,  rs_s);
}

int
calc_thresh(struct pstruct *ppst, int nstats, 
	    double Lambda, double K, double H, double *zstrim)
{
  int max_hscore;
  double ave_n1, tmp_score, z, l_fact;

  if (ppst->dnaseq == SEQT_DNA || ppst->dnaseq == SEQT_RNA) {
    ave_n1 = 5000.0;
    l_fact = 1.0;
  }
  else {
    ave_n1 = 400.0;
    l_fact = 0.7;
  }

/*  max_hscore = MAX_SSCORE; */
/*  mean expected for ppst->n0 * 400 for protein, 5000 for DNA */
/*  we want a number of offsets that is appropriate for the database size so
    far (nstats)
*/

/*
  the calculation below sets a high-score threshold using an
  ungapped lambda, but errs towards the high-score side by using
  E()=0.001 and calculating with 0.70*lambda, which is the correct for
  going from ungapped to -12/-2 gapped lambda with BLOSUM50
*/

#ifndef NORMAL_DIST
  tmp_score = 0.01/((double)nstats*K*(double)ppst->n0*ave_n1);
  tmp_score = -log(tmp_score)/(Lambda*l_fact);
  max_hscore = (int)(tmp_score+0.5);

  z = 1.0/(double)nstats;
  z = (log(z)+EULER_G)/(-PI_SQRT6);
#else
  max_hscore = 100;
  z = 5.0;
#endif
  *zstrim = 10.0*z+50.0;
  return max_hscore;
}

int
proc_hist_r(struct stat_str *sptr, int nstats,
	    struct pstruct *ppst, struct hist_str *histp,
	    int do_trim, struct pstat_str *rs)
{
  int i, max_hscore;
  double zs, ztrim;
  char s_string[128];
  struct llen_str llen;
  char *f_string;
  llen.fit_flag=1;
  llen.hist=NULL;

  max_hscore = calc_thresh(ppst, nstats, rs->ngLambda,
			   rs->ngK, rs->ngH, &ztrim);

  inithist(&llen, ppst,max_hscore);
  f_string = &(histp->stat_info[0]);

  for (i = 0; i<nstats; i++)
    addhist(&llen,sptr[i].score,sptr[i].n1, max_hscore);
  histp->entries = nstats - llen.zero_s;

  if ((llen.max_score - llen.min_score) < 10) {
    free_hist(&llen);
    llen.fit_flag = 0;
    find_zp = &find_zn;
    return proc_hist_n(sptr, nstats, ppst, histp, do_trim, rs);
  }

  fit_llen(&llen, rs); /* now we have rho, mu, rho2, mu2, mean_var
			  to set the parameters for the histogram */

  if (!llen.fit_flag) {	/* the fit failed, fall back to proc_hist_n */
    free_hist(&llen);
    find_zp = &find_zn;
    return proc_hist_n(sptr,nstats, ppst, histp, do_trim, rs);
  }

  rs->n_trimmed= rs->n1_trimmed = rs->nb_trimmed = 0;

  if (do_trim) {
    if (llen.fit_flag) {
      for (i = 0; i < nstats; i++) {
	zs = find_zr(sptr[i].score,sptr[i].escore,sptr[i].n1,sptr[i].comp, rs);
	if (zs < 20.0 || zs > ztrim) {
	  rs->n_trimmed++;
	  prune_hist(&llen,sptr[i].score,sptr[i].n1, max_hscore,
		     &(histp->entries));
	}
      }
    }

  /*  fprintf(stderr,"Z-trimmed %d entries with z > 5.0\n", rs->n_trimmed); */

    if (llen.fit_flag) fit_llens(&llen, rs);

  /*   fprintf(stderr,"Bin-trimmed %d entries in %d bins\n", rs->n1_trimmed,rs->nb_trimmed); */
  }


  free_hist(&llen);

  /* rst all the scores in the histogram */

  if (ppst->zsflag < 10) s_string[0]='\0';
  else if (ppst->zs_win > 0)
    sprintf(s_string,"(shuffled, win: %d)",ppst->zs_win);
  else strncpy(s_string,"(shuffled)",sizeof(s_string));

  inithistz(MAXHIST, histp);

  sprintf(f_string,"%s Expectation_n fit: rho(ln(x))= %6.4f+/-%6.3g; mu= %6.4f+/-%6.3f\n mean_var=%6.4f+/-%6.3f, 0's: %d Z-trim: %d  B-trim: %d in %d/%d\n Lambda= %6.4f",
	  s_string,
	  rs->rho*LN_FACT,sqrt(rs->rho_e),rs->mu,sqrt(rs->mu_e), rs->mean_var,sqrt(rs->var_e),
	  llen.zero_s, rs->n_trimmed, rs->n1_trimmed, rs->nb_trimmed, rs->nb_tot,
	  PI_SQRT6/sqrt(rs->mean_var));
    return REG_STATS;
}


int
proc_hist_n(struct stat_str *sptr, int nstats,
	    struct pstruct *ppst, struct hist_str *histp,
	    int do_trim, struct pstat_str *rs)
{
  int i, j;
  double s_score, s2_score, ssd;
  double ztrim;
  int nit, max_hscore;
  char s_string[128];
  char *f_string;

  f_string = &(histp->stat_info[0]);
  /*   db->entries = db->length = db->carry = 0; */

  max_hscore = calc_thresh(ppst, nstats, rs->ngLambda,
			   rs->ngK, rs->ngH, &ztrim);

  s_score = s2_score = 0.0;

  histp->entries = 0;

  for ( j = 0, i = 0; i < nstats; i++) {
    if (sptr[i].score > 0 && sptr[i].score <= max_hscore) {
      s_score += (ssd=(double)sptr[i].score);
      s2_score += ssd * ssd;
      histp->entries++;
      /* 
      db->length += sptr[i].n1;
      if (db->length > LONG_MAX) {
	db->carry++;
	db->length -= LONG_MAX;
      }
      */
      j++;
    }
  }

  if (j > 1 ) {
    rs->mu = s_score/(double)j;
    rs->mean_var = s2_score - (double)j * rs->mu * rs->mu;
    rs->mean_var /= (double)(j-1);
  }
  else {
    rs->mu = 50.0;
    rs->mean_var = 10.0;
  }
  
  if (rs->mean_var < 0.01) {
    rs->mean_var = (rs->mu > 1.0) ? rs->mu: 1.0;
  }

  /* now remove some scores */

  nit = 5;
  while (nit-- > 0) {
    rs->n_trimmed = 0;

    for (i=0; i< nstats; i++) {
      if (sptr[i].n1 < 0) continue;
      ssd = find_zn(sptr[i].score,sptr[i].escore,sptr[i].n1,sptr[i].comp, rs);
      if (ssd > ztrim || ssd < 20.0) {
	/*      fprintf(stderr,"removing %3d %3d %4.1f\n",
		sptr[i].score, sptr[i].n1,ssd); */
	ssd = sptr[i].score;
	s_score -= ssd;
	s2_score -= ssd*ssd;
	j--;
	rs->n_trimmed++;
	histp->entries--;
	sptr[i].n1 = -sptr[i].n1;
      }
    }

    if (j > 1 ) {
      rs->mu = s_score/(double)j;
      rs->mean_var = s2_score - (double)j * rs->mu * rs->mu;
      rs->mean_var /= (double)(j-1);
    }
    else {
      rs->mu = 50.0;
      rs->mean_var = 10.0;
    }

    if (rs->mean_var < 0.01) {
      rs->mean_var = (rs->mu > 1.0) ? rs->mu: 1.0;
    }

    if (rs->n_trimmed < LHISTC) {
      /*
	fprintf(stderr,"nprune %d at %d\n",nprune,nit);
	*/
      break;
    }
  }

  if (ppst->zsflag < 10) s_string[0]='\0';
  else if (ppst->zs_win > 0)
    sprintf(s_string,"(shuffled, win: %d)",ppst->zs_win);
  else strncpy(s_string,"(shuffled)",sizeof(s_string));

  sprintf(f_string,"%s unscaled statistics: mu= %6.4f  var=%6.4f; Lambda= %6.4f",
	  s_string, rs->mu,rs->mean_var,PI_SQRT6/sqrt(rs->mean_var));
  return AVE_STATS;
}


/*
This routine calculates the maximum likelihood estimates for the
extreme value distribution exp(-exp(-(-x-a)/b)) using the formula

	<lambda> = x_m - sum{ x[i] * exp (-x[i]<lambda>)}/sum{exp (-x[i]<lambda>)}
	<a> = -<1/lambda> log ( (1/nlib) sum { exp(-x[i]/<lambda> } )

	The <a> parameter can be transformed into and K
	of the formula: 1 - exp ( - K m n exp ( - lambda S ))
	using the transformation: 1 - exp ( -exp -(lambda S + log(K m n) ))
			1 - exp ( -exp( - lambda ( S + log(K m n) / lambda))

			a = log(K m n) / lambda
			a lambda = log (K m n)
			exp(a lambda)  = K m n 
	 but from above: a lambda = log (1/nlib sum{exp( -x[i]*lambda)})
	 so:            K m n = (1/n sum{ exp( -x[i] *lambda)})
			K = sum{}/(nlib m n )

*/

void
alloc_hist(struct llen_str *llen)
{
  int max_llen, i;
  max_llen = llen->max;

  if (llen->hist == NULL) {
    llen->hist = (int *)calloc((size_t)(max_llen+1),sizeof(int));
    llen->score_sums = (double *)calloc((size_t)(max_llen + 1),sizeof(double));
    llen->score2_sums =(double *)calloc((size_t)(max_llen + 1),sizeof(double));
    llen->score_var = (double *)calloc((size_t)(max_llen + 1),sizeof(double));
  }

  for (i=0; i< max_llen+1; i++) {
      llen->hist[i] = 0;
      llen->score_var[i] = llen->score_sums[i] = llen->score2_sums[i] = 0.0;
  }
}
  
void
free_hist(struct llen_str *llen)
{
  if (llen->hist!=NULL) {
    free(llen->score_var);
    free(llen->score2_sums);
    free(llen->score_sums);
    free(llen->hist);
    llen->hist=NULL;
  }
}

void
inithist(struct llen_str *llen, struct pstruct *ppst, int max_hscore)
{
  llen->max = MAX_LLEN;

  llen->max_score = -1;
  llen->min_score=10000;

  alloc_hist(llen);

  llen->zero_s = 0;
  llen->min_length = 10000;
  llen->max_length = 0;
}

void
addhist(struct llen_str *llen, int score, int length, int max_hscore)
{
  int llength; 
  double dscore;

  if ( score<=0 || length < LENGTH_CUTOFF) {
    llen->min_score = 0;
    llen->zero_s++;
    return ;
  }

  if (score < llen->min_score) llen->min_score = score;
  if (score > llen->max_score) llen->max_score = score;

  if (length > llen->max_length) llen->max_length = length;
  if (length < llen->min_length) llen->min_length = length;
  if (score > max_hscore) score = max_hscore;

  llength = (int)(LN_FACT*log((double)length)+0.5);

  if (llength < 0 ) llength = 0;
  if (llength > llen->max) llength = llen->max;
  llen->hist[llength]++;
  dscore = (double)score;
  llen->score_sums[llength] += dscore;
  llen->score2_sums[llength] += dscore * dscore;

  /*
  db->entries++;
  db->length += length;
  if (db->length > LONG_MAX) {db->carry++;db->length -= LONG_MAX;}
  */
}

/* histogram will go from z-scores of 20 .. 100 with mean 50 and z=10 */


void
inithistz(int mh, struct hist_str *histp )
{
  int i;

  histp->min_hist = 20;
  histp->max_hist = 120;

  histp->histint = (int)
    ((double)(histp->max_hist - histp->min_hist + 2)/(double)mh+0.5);
  histp->maxh = (int)
    ((double)(histp->max_hist - histp->min_hist + 2)/(double)histp->histint+0.5);

  if (histp->hist_a==NULL) {
    if ((histp->hist_a=(int *)calloc((size_t)histp->maxh,sizeof(int)))==
	NULL) {
      fprintf(stderr," cannot allocate %d for histogram\n",histp->maxh);
      histp->histflg = 0;
    }
    else histp->histflg = 1;
  }
  else {
    for (i=0; i<histp->maxh; i++) histp->hist_a[i]=0;
  }
}

/* fasts/f will not show any histogram */
void
addhistz(double zs, struct hist_str *histp)
{
}

void
prune_hist(struct llen_str *llen, int score, int length, int max_hscore,
	   long *entries)
{
  int llength;
  double dscore;

  if (score <= 0 || length < LENGTH_CUTOFF) return;

  if (score > max_hscore) score = max_hscore;

  llength = (int)(LN_FACT*log((double)length)+0.5);

  if (llength < 0 ) llength = 0;
  if (llength > llen->max) llength = llen->max;
  llen->hist[llength]--;
  dscore = (double)score;
  llen->score_sums[llength] -= dscore;
  llen->score2_sums[llength] -= dscore * dscore;

  (*entries)--;
  /*
  if (length < db->length) db->length -= length;
  else {db->carry--; db->length += (LONG_MAX - (unsigned long)length);}
  */
}  

/* fit_llen: no trimming
   (1) regress scores vs log(n) using weighted variance
   (2) calculate mean variance after length regression
*/

void
fit_llen(struct llen_str *llen, struct pstat_str *pr)
{
  int j;
  int n;
  int n_size;
  double x, y2, u, z;
  double mean_x, mean_y, var_x, var_y, covar_xy;
  double mean_y2, covar_xy2, var_y2, dllj;

  double sum_x, sum_y, sum_x2, sum_xy, sum_v, delta, n_w;
  
/* now fit scores to best linear function of log(n), using
   simple linear regression */
  
  for (llen->min=0; llen->min < llen->max; llen->min++)
    if (llen->hist[llen->min]) break;
  llen->min--;

  for (n_size=0,j = llen->min; j < llen->max; j++) {
    if (llen->hist[j] > 1) {
      dllj = (double)llen->hist[j];
      llen->score_var[j] = llen->score2_sums[j]/dllj
	- (llen->score_sums[j]/dllj)*(llen->score_sums[j]/dllj);
      llen->score_var[j] /= (double)(llen->hist[j]-1);
      if (llen->score_var[j] <= 0.1 ) llen->score_var[j] = 0.1;
      n_size++;
    }
  }

  pr->nb_tot = n_size;

  n_w = 0.0;
  sum_x = sum_y = sum_x2 = sum_xy = sum_v = 0;
  for (j = llen->min; j < llen->max; j++)
    if (llen->hist[j] > 1) {
      x = j + 0.5;
      dllj = (double)llen->hist[j];
      n_w += dllj/llen->score_var[j];
      sum_x +=   dllj * x / llen->score_var[j] ;
      sum_y += llen->score_sums[j] / llen->score_var[j];
      sum_x2 +=  dllj * x * x /llen->score_var[j];
      sum_xy +=  x * llen->score_sums[j]/llen->score_var[j];
    }

  if (n_size < 5 ) {
    llen->fit_flag=0;
    pr->rho = 0;
    pr->mu = sum_y/n_w;
    return;
  }
  else {
    delta = n_w * sum_x2 - sum_x * sum_x;
    if (delta > 0.001) {
      pr->rho = (n_w * sum_xy  - sum_x * sum_y)/delta;
      pr->rho_e = n_w/delta;
      pr->mu = (sum_x2 * sum_y - sum_x * sum_xy)/delta;
      pr->mu_e = sum_x2/delta;
    }
    else {
      llen->fit_flag = 0;
      pr->rho = 0;
      pr->mu = sum_y/n_w;
      return;
    }
  }

  delta = n_w * sum_x2 - sum_x * sum_x;
  pr->rho = (n_w * sum_xy  - sum_x * sum_y)/delta;
  pr->mu = (sum_x2 * sum_y - sum_x * sum_xy)/delta;

  n = 0;
  mean_x = mean_y = mean_y2 = 0.0;
  var_x = var_y = 0.0;
  covar_xy = covar_xy2 = 0.0;

  for (j = llen->min; j <= llen->max; j++) 
   if (llen->hist[j] > 1 ) {
      n += llen->hist[j];
      x = (double)j + 0.5;
      mean_x += (double)llen->hist[j] * x;
      mean_y += llen->score_sums[j];
      var_x += (double)llen->hist[j] * x * x;
      var_y += llen->score2_sums[j];
      covar_xy += x * llen->score_sums[j];
    }
  mean_x /= n; mean_y /= n;
  var_x = var_x / n - mean_x * mean_x;
  var_y = var_y / n - mean_y * mean_y;
  
  covar_xy = covar_xy / n - mean_x * mean_y;
/*
  pr->rho = covar_xy / var_x;
  pr->mu = mean_y - pr->rho * mean_x;
*/
  mean_y2 = covar_xy2 = var_y2 = 0.0;
  for (j = llen->min; j <= llen->max; j++) 
    if (llen->hist[j] > 1) {
      x = (double)j + 0.5;
      u = pr->rho * x + pr->mu;
      y2 = llen->score2_sums[j] - 2.0 * llen->score_sums[j] * u + llen->hist[j] * u * u;
/*
      dllj = (double)llen->hist[j];
      fprintf(stderr,"%.2f\t%d\t%g\t%g\n",x/LN_FACT,llen->hist[j],
	      llen->score_sums[j]/dllj,y2/dllj);
*/
      mean_y2 += y2;
      var_y2 += y2 * y2;
      covar_xy2 += x * y2;
      /*      fprintf(stderr,"%6.1f %4d %8d %8d %7.2f %8.2f\n",
	      x,llen->hist[j],llen->score_sums[j],llen->score2_sums[j],u,y2); */
    }
  
  pr->mean_var = mean_y2 /= (double)n;
  covar_xy2 = covar_xy2 / (double)n - mean_x * mean_y2;

  if (pr->mean_var <= 0.01) {
    llen->fit_flag = 0;
    pr->mean_var = (pr->mu > 1.0) ? pr->mu: 1.0;
  }

  /*
  fprintf(stderr," rho1/mu1: %.4f/%.4f mean_var %.4f\n",
	  pr->rho*LN_FACT,pr->mu,pr->mean_var);
  */
  if (n > 1) pr->var_e = (var_y2/n - mean_y2 * mean_y2)/(n-1);
  else pr->var_e = 0.0;

  if (llen->fit_flag) {
    pr->rho2 = covar_xy2 / var_x;
    pr->mu2 = pr->mean_var - pr->rho2 * mean_x;
  }
  else {
    pr->rho2 = 0;
    pr->mu2 = pr->mean_var;
  }

  if (pr->rho2 < 0.0 )
    z = (pr->rho2 * LN_FACT*log((double)llen->max_length) + pr->mu2 > 0.0) ? llen->max_length : exp((-1.0 - pr->mu2 / pr->rho2)/LN_FACT);
  else z =  pr->rho2 ? exp((1.0 - pr->mu2 / pr->rho2)/LN_FACT) : LENGTH_CUTOFF;
  if (z < 2*LENGTH_CUTOFF) z = 2*LENGTH_CUTOFF;

  pr->var_cutoff = pr->rho2 * LN_FACT*log(z) + pr->mu2;
}

/* fit_llens: trim high variance bins
   (1) regress scores vs log(n) using weighted variance
   (2) regress residuals vs log(n)
   (3) remove high variance bins
   (4) calculate mean variance after length regression
*/

void
fit_llens(struct llen_str *llen, struct pstat_str *pr)
{
  int j;
  int n, n_u2;
  double x, y, y2, u, u2, v, z;
  double mean_x, mean_y, var_x, var_y, covar_xy;
  double mean_y2, covar_xy2;
  double mean_u2, mean_3u2, dllj;
  double sum_x, sum_y, sum_x2, sum_xy, sum_v, delta, n_w;

/* now fit scores to best linear function of log(n), using
   simple linear regression */
  
  for (llen->min=0; llen->min < llen->max; llen->min++)
    if (llen->hist[llen->min]) break;
  llen->min--;

  for (j = llen->min; j < llen->max; j++) {
    if (llen->hist[j] > 1) {
      dllj = (double)llen->hist[j];
      llen->score_var[j] = (double)llen->score2_sums[j]/dllj
	- (llen->score_sums[j]/dllj)*(llen->score_sums[j]/dllj);
      llen->score_var[j] /= (double)(llen->hist[j]-1);
      if (llen->score_var[j] <= 1.0 ) llen->score_var[j] = 1.0;
    }
  }
	  
  n_w = 0.0;
  sum_x = sum_y = sum_x2 = sum_xy = sum_v = 0;
  for (j = llen->min; j < llen->max; j++)
    if (llen->hist[j] > 1) {
      x = j + 0.5;
      dllj = (double)llen->hist[j];
      n_w += dllj/llen->score_var[j];
      sum_x +=   dllj * x / llen->score_var[j] ;
      sum_y += llen->score_sums[j] / llen->score_var[j];
      sum_x2 +=  dllj * x * x /llen->score_var[j];
      sum_xy +=  x * llen->score_sums[j]/llen->score_var[j];
    }

  delta = n_w * sum_x2 - sum_x * sum_x;
  pr->rho = (n_w * sum_xy  - sum_x * sum_y)/delta;
  pr->mu = (sum_x2 * sum_y - sum_x * sum_xy)/delta;

/* printf(" rho1/mu1: %.2f/%.2f\n",pr->rho*LN_FACT,pr->mu); */

  n = 0;
  mean_x = mean_y = mean_y2 = 0.0;
  var_x = var_y = 0.0;
  covar_xy = covar_xy2 = 0.0;

  for (j = llen->min; j <= llen->max; j++) 
    if (llen->hist[j] > 1 ) {
      n += llen->hist[j];
      x = (double)j + 0.5;
      dllj = (double)llen->hist[j];
      mean_x += dllj * x;
      mean_y += llen->score_sums[j];
      var_x += dllj * x * x;
      var_y += llen->score2_sums[j];
      covar_xy += x * llen->score_sums[j];
    }
  mean_x /= n; mean_y /= n;
  var_x = var_x / n - mean_x * mean_x;
  var_y = var_y / n - mean_y * mean_y;
  
  covar_xy = covar_xy / n - mean_x * mean_y;
/*  pr->rho = covar_xy / var_x;
  pr->mu = mean_y - pr->rho * mean_x;
*/

  mean_y2 = covar_xy2 = 0.0;
  for (j = llen->min; j <= llen->max; j++) 
    if (llen->hist[j] > 1) {
      x = (double)j + 0.5;
      u = pr->rho * x + pr->mu;
      y2 = llen->score2_sums[j] - 2 * llen->score_sums[j] * u + llen->hist[j] * u * u;
      mean_y2 += y2;
      covar_xy2 += x * y2;
    }
  
  mean_y2 /= n;
  covar_xy2 = covar_xy2 / n - mean_x * mean_y2;
  pr->rho2 = covar_xy2 / var_x;
  pr->mu2 = mean_y2 - pr->rho2 * mean_x;

  if (pr->rho2 < 0.0 )
    z = (pr->rho2 * LN_FACT*log((double)llen->max_length) + pr->mu2 > 0.0) ? llen->max_length : exp((-1.0 - pr->mu2 / pr->rho2)/LN_FACT);
  else z =  pr->rho2 ? exp((1.0 - pr->mu2 / pr->rho2)/LN_FACT) : LENGTH_CUTOFF;
  if (z < 2* LENGTH_CUTOFF) z = 2*LENGTH_CUTOFF;

  pr->var_cutoff = pr->rho2*LN_FACT*log(z) + pr->mu2;

/*  fprintf(stderr,"\nminimum allowed predicted variance (%0.2f) at n = %.0f\n",
	 pr->var_cutoff,z);
*/
  mean_u2 = 0.0;
  n_u2 = 0;
  for ( j = llen->min; j < llen->max; j++) {
    y = j+0.5;
    dllj = (double)llen->hist[j];
    x = pr->rho * y + pr->mu;
    v = pr->rho2 * y + pr->mu2;
    if (v < pr->var_cutoff) v = pr->var_cutoff;
    if (llen->hist[j]> 1) {
      u2 =  (llen->score2_sums[j] - 2 * x * llen->score_sums[j] + dllj * x * x) - v*dllj;
      mean_u2 += llen->score_var[j] = u2*u2/(llen->hist[j]-1);
      n_u2++;
      /*      fprintf(stderr," %d (%d) u2: %.2f v*ll: %.2f %.2f\n",
	      j,llen->hist[j],u2,v*dllj,sqrt(llen->score_var[j])); */
    }
    else llen->score_var[j] = -1.0;
  }

  mean_u2 = sqrt(mean_u2/(double)n_u2);
  /* fprintf(stderr," mean s.d.: %.2f\n",mean_u2); */

  mean_3u2 = mean_u2*3.0;

  for (j = llen->min; j < llen->max; j++) {
    if (llen->hist[j] <= 1) continue;
    if (sqrt(llen->score_var[j]) > mean_3u2) {
      /*      fprintf(stderr," removing %d %d %.2f\n",
	     j, (int)(exp((double)j/LN_FACT)-0.5),
	     sqrt(llen->score_var[j]));
	     */
      pr->nb_trimmed++;
      pr->n1_trimmed += llen->hist[j];
      llen->hist[j] = 0;
    }
  }
  fit_llen(llen, pr);
}


double find_z(int score, double escore, int length, double comp, struct pstat_str *pu) {
  return find_zp(score, escore, length, comp, pu);
}

/* REG_STATS - Z() from rho/mu/mean_var */
double find_zr(int score, double escore, int length, double comp, 
	      struct pstat_str *rs)
{
  double log_len, z;
  
  if (score <= 0) return 0.0;
  if ( length < LENGTH_CUTOFF) return 0.0;

  log_len = LN_FACT*log((double)(length));
/*  var = rs->rho2 * log_len + rs->mu2;
  if (var < rs->var_cutoff) var = rs->var_cutoff;
*/

  z = ((double)score - rs->rho * log_len - rs->mu) / sqrt(rs->mean_var);

  return (50.0 + z*10.0);
}

double find_zt(int score, double escore, int length, double comp, 
	      struct pstat_str *rs)
{
  if (!rs->eval_is_pval) {escore /= rs->zdb_size;}

  if (escore > 0.0) return -log(escore)/M_LN2;
  else return 744.440071/M_LN2;
}

double find_zn(int score, double escore, int length, double comp,
	      struct pstat_str *rs)
{
  double z;
  
  z = ((double)score - rs->mu) / sqrt(rs->mean_var);

  return (50.0 + z*10.0);
}

/* computes E value for a given z value, assuming extreme value distribution */
double
z_to_E(double zs, long entries, struct db_str db)
{
  double e, n;

  /*  if (db->entries < 5) return (double)db.entries; */
  if (entries < 1) { n = db.entries;}
  else {n = entries;}

  if (zs > ZS_MAX) return 0.0;

  e = exp(-PI_SQRT6 * zs - EULER_G);
  return n * (e > .01 ? 1.0 - exp(-e) : e);
}

double
zs_to_p(double zs)
{
  return zs;
}

/* this version assumes the probability is in the ->zscore variable,
   which is provided by this file after last_scale()
*/

double
zs_to_bit(double zs, int n0, int n1)
{
  return zs+log((double)(n0*n1))/M_LN2 ;
}

/* computes E-value for a given z value, assuming extreme value distribution */
double
zs_to_E(double zs,int n1, int dnaseq, long entries, struct db_str db)
{
  double e, z, k;

  /*  if (db->entries < 5) return 0.0; */

  if (zs > ZS_MAX ) return 0.0;

  if (entries < 1) entries = db.entries;

  if (dnaseq == SEQT_DNA || dnaseq == SEQT_RNA) {
    k = (double)db.length /(double)n1;
    if (db.carry > 0) { k *= (double)db.carry * (double)LONG_MAX;}
  }
  else k = (double)entries;

  if (k < 1.0) k = 1.0;

  zs *= M_LN2;
  if ( zs > 100.0) e = 0.0;
  else e =  exp(-zs);
  return k * e;
}

/* computes E-value for a given z value, assuming extreme value distribution */
double
E_to_zs(double E, long entries)
{
  double e, z;
  int error;

  e = E/(double)entries;

#ifndef NORMAL_DIST
  z = (log(e)+EULER_G)/(-PI_SQRT6);
  return z*10.0+50.0;
#else
  z = np_to_z(1.0-e,&error);

  if (!error) return z*10.0+50.0;
  else return 0.0;
#endif
}

/* computes 1.0 - E value for a given z value, assuming extreme value
   distribution */
double
zs_to_Ec(double zs, long entries)
{
  double e, z;

  if (entries < 5) return 0.0;

  z = (zs - 50.0)/10.0;

  if (z > ZS_MAX) return 1.0;

  e =  exp(-PI_SQRT6 * z - EULER_G);
  return (double)entries * (e > .01 ?  exp(-e) : 1.0 - e);
}

int
E1_to_s(double e_val, int n0, int n1, int db_size,
	void *pu) {
  double mp, np, a_n0, a_n0f, a_n1;
  double zs, log_len, p_val;
  int score;

  if (n1 < LENGTH_CUTOFF) return 0;

  score = -log(e_val)/log(10.0);

#ifndef NORMAL_DIST
  if (score < 0) score = 0;
#endif
  return score;
}

void
sort_escore(double *v, int n)
{
  int gap, i, j;
  double dtmp;
	
  for (gap=n/2; gap>0; gap/=2) {
    for (i=gap; i<n; i++) {
      for (j=i-gap; j>=0; j -= gap) {
	if (v[j] <= v[j+gap]) break;
	dtmp = v[j];
	v[j] = v[j+gap];
	v[j+gap] = dtmp;
      }
    }
  }
}

/* scale_tat - compute 'a', 'b', 'c' coefficients for scaling fasts/f
   escores 
   5-May-2003 - also calculate index for high ties
*/
void
scale_tat(double *escore, int nstats,
	  long db_entries, int do_trim,
	  struct pstat_str *rs)
{
  int i, j, k, start;
  double *x, *lnx, *lny;

  /*   sort_escore(escore, nstats); */

  while (*escore<0.0) {escore++; nstats--; }

  x = (double *) calloc(nstats, sizeof(double));
  if(x == NULL) {
    fprintf(stderr, "Couldn't calloc tatE/x\n");
    exit(1);
  }

  lnx = (double *) calloc(nstats,sizeof(double));
  if(lnx == NULL) {
    fprintf(stderr, "Couldn't calloc tatE/lnx\n");
    exit(1);
  }
  
  lny = (double *) calloc(nstats,sizeof(double));
  if(lny == NULL) {
    fprintf(stderr, "Couldn't calloc tatE/lny\n");
    exit(1);
  }
  
  for(i = 0 ; i < nstats ; ) {

    lny[i] = log(escore[i]);

    for(j = i+1 ; j < nstats ; j++) {
	if(escore[j] != escore[i]) break;
    }

    x[i] = ((((double)i + (double)(j - i - 1)/2.0)*(double)nstats/(double)db_entries)+1.0)/(double)nstats;
    lnx[i] = log(x[i]);

    for(k = i+1 ; k < j ; k++) {
      lny[k]=lny[i];
      x[k] = x[i];
      lnx[k]=lnx[i];
    }
    i = k;
  }

  if (!do_trim) {
    start = 0;
  } else {
    start = 0.05 * (double) nstats;
    start = start > 500 ? 500 : start;
  }

  linreg(lny, x, lnx, nstats, &rs->tat_a, &rs->tat_b, &rs->tat_c, start);

  /* I have the coefficients I need - a, b, c; free arrays */

  free(lny);
  free(lnx);
  free(x);

  /* calculate tie_j - the index below which all scores are considered
     positional ties */

  rs->tie_j = 0.005 * db_entries;
}

void
linreg(double *lny, double *x, double *lnx, int n,
       double *a, double *b, double *c, int start) {

  double yf1, yf2, yf3;
  double f1f1, f1f2, f1f3;
  double f2f2, f2f3;
  double f3f3, delta;

  int i;

  yf1 = yf2 = yf3 = 0.0;
  f1f1 = f1f2 = f1f3 = f2f2 = f2f3 = f3f3 = 0.0;

  for(i = start; i < n; i++) {
    yf1 += lny[i] * lnx[i];
    yf2 += lny[i] * x[i];
    yf3 += lny[i];

    f1f1 += lnx[i] * lnx[i];
    f1f2 += lnx[i] * x[i];
    f1f3 += lnx[i];

    f2f2 += x[i] * x[i];
    f2f3 += x[i];

    f3f3 += 1.0;
  }

  delta = det(f1f1, f1f2, f1f3, f1f2, f2f2, f2f3, f1f3, f2f3, f3f3);

  *a = det(yf1, f1f2, f1f3, yf2, f2f2, f2f3, yf3, f2f3, f3f3) / delta;
  *b = det(f1f1, yf1, f1f3, f1f2, yf2, f2f3, f1f3, yf3, f3f3) / delta;
  *c = det(f1f1, f1f2, yf1, f1f2, f2f2, yf2, f1f3, f2f3, yf3) / delta;

}

double det(double a11, double a12, double a13,
	   double a21, double a22, double a23,
	   double a31, double a32, double a33)
{
  double result;

  result = a11 * (a22 * a33 - a32 * a23);
  result -= a12 * (a21 * a33 - a31 * a23);
  result += a13 * (a21 * a32 - a31 * a22);
    
  return result;
}

void
last_stats(const unsigned char *aa0, int n0,
	   struct stat_str *sptr, int nstats,
	   struct beststr **bestp_arr, int nbest,
	   const struct mngmsg *m_msg, struct pstruct *ppst,
	   struct hist_str *histp, struct pstat_str **rs_sp)
{
  double *obs_escore;
  int i, nobs, nobs_t, do_trim;
  long db_entries;
  struct pstat_str *rs_s;

  if (*rs_sp == NULL) {
    if ((rs_s=(struct pstat_str *)calloc(1,sizeof(struct pstat_str)))==NULL) {
      fprintf(stderr," cannot allocate rs_s: %ld\n",sizeof(struct pstat_str));
      exit(1);
    }
    else *rs_sp = rs_s;
  }
  else rs_s = *rs_sp;
    
  histp->entries = 0;

  sortbeste(bestp_arr,nbest);

  rs_s->spacefactor = 
    calc_spacefactor(aa0, n0, m_msg->nm0,ppst->nsq);

  if (ppst->zsflag >= 1 && ppst->zsflag <= 4) {
    if (m_msg->escore_flg) {
      nobs = nbest;
      do_trim = 1;
    }
    else {
      nobs = nstats;
      do_trim = 0;
    }

    if ((obs_escore = (double *)calloc(nobs,sizeof(double)))==NULL) {
      fprintf(stderr," cannot allocate obs_escore[%d]\n",nbest);
      exit(1);
    }

    if (m_msg->escore_flg) {
      for (i=nobs=0; i<nbest; i++) {
	if (bestp_arr[i]->rst.escore<= 1.00)
	  obs_escore[nobs++]=bestp_arr[i]->rst.escore;
      }
      /*
      nobs_t = nobs;
      for (i=0; i<nbest; i++) {
	if (bestp_arr[i]->rst.escore >= 0.99 &&
	    bestp_arr[i]->rst.escore <= 1.00)
	  obs_escore[nobs++]=bestp_arr[i]->rst.escore;
      }
      */
      db_entries = m_msg->db.entries;
    }
    else {
      for (i=nobs=0; i<nstats; i++) {
	if (sptr[i].escore <= 1.00 ) obs_escore[nobs++]=sptr[i].escore;
      }
      /*
      nobs_t = nobs;
      for (i=0; i<nstats; i++) {
	if (sptr[i].escore >= 0.99 &&
	    sptr[i].escore <= 1.0) obs_escore[nobs++]=sptr[i].escore;
      }
      */
      db_entries = nobs;
/*    db_entries = m_msg->db.entries;*/
    }

    sortbesto(obs_escore,nobs);
    if (nobs > 100) {
      scale_tat(obs_escore,nobs,db_entries,do_trim,rs_s);
      rs_s->have_tat=1;
      sprintf(histp->stat_info,"scaled Tatusov statistics (%d): tat_a: %6.4f tat_b: %6.4f tat_c: %6.4f",
	      nobs,rs_s->tat_a, rs_s->tat_b, rs_s->tat_c);
    }
    else {
      rs_s->have_tat=0;
      sprintf(histp->stat_info,"Space_factor %.4g scaled statistics",
	    rs_s->spacefactor);
    }
    free(obs_escore);
  }
  else {
    rs_s->have_tat=0;
    histp->stat_info[0] = '\0';
  }
  if (rs_s->have_tat) {
    find_zp = &find_zt;
  }

}

/* scale_scores() takes the best (real) scores and re-scales them;
   beststr bptr[] must be sorted */

void
scale_scores(struct beststr **bptr, int nbest, struct db_str db,
	     struct pstruct *ppst, struct pstat_str *rs)
{
  int i, j, k;
  double obs, r_a, r_b, r_c;

  /* this scale function absolutely requires that the results be sorted
     before it is used */

  sortbeste(bptr,nbest);

  if (!rs->have_tat) {
    for (i=0; i<nbest; i++) {
      bptr[i]->rst.escore *= rs->spacefactor;
    }
  }
  else {

    /* here if more than 1000 scores */

    r_a = rs->tat_a; r_b = rs->tat_b; r_c = rs->tat_c;

  /* the problem with scaletat is that the E() value is related to
     ones position in the list of top scores - thus, knowing the score
     is not enough - one must know the rank */

    for(i = 0 ; i < nbest ; ) {
      /* take the bottom 0.5%, and the ties, and treat them all the same */
      j = i + 1;
      while (j< nbest && 
	     (j <= (0.005 * db.entries) || bptr[j]->rst.escore == bptr[i]->rst.escore)
	     ) {
	j++;
      }

      /* observed frequency */
      obs = ((double)i + ((double)(j - i - 1)/ 2.0) + 1.0)/(double)db.entries;

      /* make certain ties all have the same correction */
      for (k = i ; k < j ; k++) {
	bptr[k]->rst.escore *= obs/exp(r_a*log(obs) + r_b*obs + r_c);
      }
      i = k;
    }
  }

  for (i=0; i<nbest; i++) {
    if(bptr[i]->rst.escore > 0.01)
      bptr[i]->rst.escore = 1.0 - exp(-bptr[i]->rst.escore);
    if (bptr[i]->rst.escore > 0.0)
      bptr[i]->zscore = -log(bptr[i]->rst.escore)/M_LN2;
    else 
      bptr[i]->zscore = 744.440071/M_LN2;
    bptr[i]->rst.escore *= ppst->zdb_size;
  }

  rs->zdb_size = ppst->zdb_size;
  rs->eval_is_pval = 0;
}

double scale_one_score (int ipos, double escore,
                        struct db_str db,
                        struct pstat_str *rs) {
  double obs;
  double a, b, c;

  if (!rs->have_tat)
    return escore * rs->spacefactor;

  if (ipos < rs->tie_j) ipos = rs->tie_j/2;

  a = rs->tat_a; b = rs->tat_b; c = rs->tat_c;

  obs = ((double)ipos + 1.0)/(double)db.entries;

  escore *= obs/exp(a*log(obs) + b*obs + c);

  return escore;
}

double calc_spacefactor(const unsigned char *aa0, int n0,
			int nm0, int nsq) {

#if !defined(FASTF)
  return pow(2.0, (double) nm0) - 1.0;
#else

  int i, j, n, l, nr, bin, k;
  int nmoff;
  int **counts;
  int **factors;
  double tmp, result = 0.0;

  nmoff = (n0 - nm0 + 1)/nm0+1;

  counts = (int **) calloc(nsq, sizeof(int *));
  if(counts == NULL) {
    fprintf(stderr, "couldn't calloc counts array!\n");
    exit(1);
  }

  counts[0] = (int *) calloc(nsq * (nmoff - 1), sizeof(int));
  if(counts[0] == NULL) {
    fprintf(stderr, "couldn't calloc counts array!\n");
    exit(1);
  }

  for(i = 0 ; i < nsq ; i++) {
    counts[i] = counts[0] + (i * (nmoff - 1));
  }

  for(i = 0 ; i < nm0 ; i++) {
    for(j = 0 ; j < (nmoff - 1) ; j++) {
      counts[ aa0[nmoff * i + j] ] [ j ] ++;
    }
  }

  factors = (int **) calloc(nm0 + 1, sizeof(int *));
  if(factors == NULL) {
    fprintf(stderr, "Couldn't calloc factors array!\n");
    exit(1);
  }

  factors[0] = (int *) calloc((nm0 + 1) * (nmoff - 1), sizeof(int));
  if(factors[0] == NULL) {
    fprintf(stderr, "Couldn't calloc factors array!\n");
    exit(1);
  }

  for(i = 0 ; i <= nm0 ; i++) {
    factors[i] = factors[0] + (i * (nmoff - 1));
  }

  /*
    this algorithm was adapted from the GAP4 library's NrArrangement function:
    The GAP Group, GAP --- Groups, Algorithms, and Programming,
    Version 4.1; Aachen, St Andrews, 1999.
    (http://www-gap.dcs.st-and.ac.uk/ gap)
  */

  /* calculate K factors for each column in query: */
  for(j = 0 ; j < (nmoff - 1) ; j++) {

    /* only one way to select 0 elements */
    factors[0][j] = 1;

    /* for each of the possible elements in this column */
    for(n = 0 ; n < nsq ; n++) {

      /* if there aren't any of these, skip it */
      if(counts[n][j] == 0) { continue; }

      /* loop over the possible lengths of the arrangement: K..0 */
      for(l = nm0 ; l >= 0 ; l--) {
	nr = 0;
	bin = 1;

	/*
	  compute the number of arrangements of length <l>
	  using only the first <n> elements of <mset>
	*/
	for(i = 0, k = min(counts[n][j], l); i <= k ; i++) {

	  /* 
	     add the number of arrangements of length <l>
	     that consist of <l>-<i> of the first <n>-1 elements
	     and <i> copies of the <n>th element
	  */
	  nr += bin * factors[l-i][j];
	  bin = (int) ((float) bin * (float) (l - i) / (float) (i + 1));
	}

	factors[l][j] = nr;
      }
    }
  }

  result = 0.0;
  for(i = 1 ; i <= nm0 ; i++) {
    tmp = 1.0;
    for(j = 0 ; j < (nmoff - 1) ; j++) {
      tmp *= (double) factors[i][j];
    }
    tmp /= factorial(i, 1);
    result += tmp;
  }
  
  free(counts[0]);
  free(counts);
  free(factors[0]);
  free(factors);

  return result;
#endif
}

void sortbesto (double *obs, int nobs)
{
  int gap, i, j, k;
  double v;
  int incs[16] = { 1391376, 463792, 198768, 86961, 33936,
		   13776, 4592, 1968, 861, 336, 
		   112, 48, 21, 7, 3, 1 };

  for ( k = 0; k < 16; k++)
    for (gap = incs[k], i=gap; i < nobs; i++) {
      v = obs[i];
      j = i;
      while ( j >= gap && obs[j-gap] > v) {
	obs[j] = obs[j - gap];
	j -= gap;
      }
      obs[j] = v;
    }
}

/* print out all pstat_str info for independent calculation */
void
pstat_info(char *info_str, int info_str_n, char *comment, struct pstat_str *pu) {
  char pstat_buf[MAX_STR];

  sprintf(pstat_buf,"%s zsflag: %d\n",comment,pu->zsflag);
  SAFE_STRNCPY(info_str,pstat_buf,info_str_n);
  sprintf(pstat_buf,"%s ngLambda: %g; ngK: %g; ngH: %g\n",comment,pu->ngLambda,pu->ngK,pu->ngH);
  SAFE_STRNCAT(info_str,pstat_buf,info_str_n);

  sprintf(pstat_buf,"%s rho: %g; rho_e: %g; mu: %g; mu_e: %g;\n",comment,
	  pu->rho,pu->rho_e,pu->mu,pu->mu_e);
  SAFE_STRNCAT(info_str,pstat_buf,info_str_n);

  sprintf(pstat_buf,"%s mean_var: %g; var_e: %gg\n",comment,
	  pu->mean_var, pu->var_e);
  SAFE_STRNCAT(info_str,pstat_buf,info_str_n);

  sprintf(pstat_buf,"%s rho2: %g; mu2: %g; var_cutoff: %g\n",comment,
	  pu->rho2, pu->mu2,pu->var_cutoff);
  SAFE_STRNCAT(info_str,pstat_buf,info_str_n);

  sprintf(pstat_buf,"%s n_trimmed: %d; n1_trimmed: %d; nb_trimmed: %d; nb_tot: %d\n",comment,
	  pu->n_trimmed, pu->n1_trimmed,pu->nb_trimmed, pu->nb_tot);
  SAFE_STRNCAT(info_str,pstat_buf,info_str_n);

  sprintf(pstat_buf,"%s tat_a: %g; tat_b: %g; tat_c: %g; spacefactor: %g\n",comment,
	  pu->tat_a, pu->tat_b,pu->tat_c, pu->spacefactor);
  SAFE_STRNCAT(info_str,pstat_buf,info_str_n);

  sprintf(pstat_buf,"%s have_tat: %d; tie_j: %d; eval_is_pval: %d; zdb_size: %ld\n",comment,
	  pu->have_tat,pu->tie_j,pu->eval_is_pval,pu->zdb_size);
  SAFE_STRNCAT(info_str,pstat_buf,info_str_n);
}
