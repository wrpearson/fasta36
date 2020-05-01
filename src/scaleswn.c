/* scaleswn.c */

/* $Id: scaleswn.c 1245 2013-12-18 18:19:38Z wrp $ */

/* copyright (c) 1995, 1996, 2000, 2014 by William R. Pearson and The
   Rector & Visitors of the University of Virginia */

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

/* Provide statistical estimates using an extreme value distribution 

	This code provides multiple methods for scaling sequence
	similarity scores to correct for length effects.

	Currently, six methods are available:

	ppst->zsflag = 0 - no scaling  (AVE_STATS)
	ppst->zsflag = 1 - regression-scaled scores (REG_STATS)
	ppst->zsflag = 2 - (revised) MLE Lambda/K scaled scores (MLE_STATS)
	ppst->zsflag = 3 - scaling using Altschul's parameters (AG_STATS)
	ppst->zsflag = 4 - regression-scaled with iterative outlier removal (REGI_STATS)
	ppst->zsflag = 5 = like 1, but length scaled variance (REG2_STATS)
	ppst->zsflag = 6 = like 2, but uses lambda composition/scale (MLE2_STATS)
	ppst->zsflag = 11 = 10 + 1 - use random shuffles, method 1
*/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

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
struct rstat_str {
  double rho, rho_e, mu, mu_e, mean_var, var_e;  /* var_e:std. error of var */
  double mean_var_sqrt;		/* sqrt(mean_var) - used frequently */
/* used by REG2_STATS */
  double rho2, mu2, var_cutoff;
  int n_trimmed; /* excluded because of high z-score */
  int n1_trimmed, nb_trimmed, nb_tot; /* excluded because of bin */
};

/* used by AG_STATS, MLE_STATS */
struct ag_stat_str {
  double K, Lambda, H, a_n0f, a_n0;
};

/* used by MLE2_STATS */
struct mle2_stat_str {
  double a_n0;
  double mle2_a0, mle2_a1, mle2_a2, mle2_b1;
  double ave_comp, max_comp, ave_H;
};

struct pstat_str {
  int zsflag;
  double ngLambda, ngK, ngH;
  double ave_n1;
  double sample_fract;
  double zs_off;	/* z-score offset from sampling */
  union {
    struct rstat_str rg;
    struct ag_stat_str ag;
    struct mle2_stat_str m2;
  } r_u;
};

double find_z(int score, double escore, int len, double comp, struct pstat_str *);

#define AVE_STATS 0	/* no length effect, only mean/variance */
double find_zn(int score, double escore, int len, double comp, struct pstat_str *);

static void st_sort (struct stat_str *v, int n);

static int
calc_thresh(struct pstruct *ppst, struct stat_str *sptr, int nstats,
	    struct score_count_s s_info,
	    double Lambda, double K, double H, double *zstrim);

int proc_hist_n(struct stat_str *sptr, int nstats, struct score_count_s s_info,
		struct pstruct *ppst, struct hist_str *histp, int do_trim,
		struct pstat_str *);

#define REG_STATS 1	/* length-regression scaled */
#define REGI_STATS 4	/* length regression, iterative */
double find_zr(int score, double escore, int len, double comp, struct pstat_str *);
int proc_hist_r(struct stat_str *sptr,  int nstats, struct score_count_s s_info,
		struct pstruct *ppst, struct hist_str *histp,
		int do_trim, struct pstat_str *pu);

#define MLE_STATS 2	/* MLE for lambda, K */
double find_ze(int score, double escore, int len, double comp, struct pstat_str *);
int proc_hist_ml(struct stat_str *sptr,  int nstats, struct score_count_s s_info,
		 struct pstruct *ppst, struct hist_str *histp, int do_trim,
		 struct pstat_str *);

#define AG_STATS 3	/* Altschul-Gish parameters */
double find_za(int score, double escore, int len, double comp, struct pstat_str *);
int proc_hist_a(struct stat_str *sptr, int nstats, struct score_count_s s_info,
		struct pstruct *ppst, struct hist_str *histp, int do_trim,
		struct pstat_str *);

#ifdef NORMAL_DIST
int proc_hist_anorm(struct stat_str *sptr, int nstats, struct score_count_s s_info,
		struct pstruct *ppst, struct hist_str *histp, int do_trim,
		struct pstat_str *);
#endif

#define REG2_STATS 5	/* length regression on mean + variance */
double find_zr2(int score, double escore, int len, double comp, struct pstat_str *);
int proc_hist_ri(struct stat_str *sptr, int nstats, struct score_count_s s_info,
		 struct pstruct *ppst, struct hist_str *histp, int do_trim,
		 struct pstat_str *);

#define MLE2_STATS 6	/* MLE stats using comp(lambda) */
double find_ze2(int score, double escore, int length, double comp, struct pstat_str *);
int proc_hist_ml2(struct stat_str *sptr, int nstats, struct score_count_s s_info,
		  struct pstruct *ppst, struct hist_str *histp, int do_trim,
		  struct pstat_str *);

#ifdef USE_LNSTATS
#define LN_STATS 2
double find_zl(int score, double escore, int len, double comp, struct pstat_str *);
int proc_hist_ln(struct stat_str *sptr, int nstats, struct score_count_s s_info,
		 struct pstruct *ppst, struct hist_str *histp, int do_trim,
		 struct pstat_str *);
#endif

double (*find_z_arr[7])(int score, double escore, int len, double comp, struct pstat_str *) = {
  find_zn, 	/* AVE_STATS zsflag ==0 */
  find_zr,	/* REG_STATS zsflag ==1 */
  find_ze,	/* MLE_STATS, zsflag==2 */
#ifndef NORMAL_DIST
  find_za,	/* AG_STATS, zsflag==3 */
#else
  find_zn, 	/* AVE_STATS zsflag ==0 */
#endif  
  find_zr,	/* REGI_STATS, zsflag==4 */
  find_zr2,	/* REG2_STATS, zsflag==5 */
  find_ze2,	/* MLE2_STATS, zsflag==6 */
};

/* print out all pstat_str info for independent calculation */
void
pstat_info(char *info_str, int info_str_n, char *comment, struct pstat_str *pu);

/* scaleswn.c local variables that belong in their own structure */

struct llen_str {
  int min, max;
  int max_score, min_score;
  int *hist;
  double *score_sums, *score2_sums;
  double *score_var;
  int max_length, min_length, zero_s;
  int fit_flag;
};

/* llen bins */
static void inithist(struct llen_str *, struct pstruct *, int);
static void free_hist( struct llen_str *);
static void addhist(struct llen_str *, int, int, int);
static void prune_hist(struct llen_str *, int, int, int, long *);
/* final display histogram */
void inithistz(int, struct hist_str *histp, double zs_off);
void addhistz(double zs, struct hist_str *histp);
void addhistzp(double zs, struct hist_str *histp);

/* calculate rho (slope), mu (intercept) REG_STATS,
   mean_var, rho2, mu2 for variance (REG2_STATS),
   minimal iterative exclusion
*/
static void fit_llen(struct llen_str *, struct rstat_str *);
/* calculate rho (slope), mu (intercept) REG_STATS,
   mean_var, rho2, mu2 for variance (REG2_STATS),
   
*/
static void fit_llen2(struct llen_str *, struct rstat_str *);

/* calculate rho (slope), mu (intercept) REG_STATS,
   mean_var, rho2, mu2 for variance (REG2_STATS) */
static void fit_llens(struct llen_str *, struct rstat_str *);

extern void sortbeste(struct beststr **bptr, int nbest);
extern void sortbestz(struct beststr **bptr, int nbest);

/* void set_db_size(int, struct db_str *, struct hist_str *); */

#ifdef DEBUG
FILE *tmpf;
#endif

int
process_hist(struct stat_str *sptr, int nstats, 
	     const struct mngmsg *m_msp,
	     struct pstruct *ppst,
	     struct hist_str *histp,
	     struct pstat_str **ps_sp,
	     struct score_count_s *s_info,
	     int do_hist)
{
  int zsflag, r_zsflag, do_trim, i;
  struct pstat_str *ps_s;

  if (ppst->zsflag < 0) {	/* no statistics */
    *ps_sp = NULL;
    return ppst->zsflag;
  }

  /* need a pstat_str if its NULL */
  if (*ps_sp == NULL) {
    if ((ps_s=(struct pstat_str *)calloc(1,sizeof(struct pstat_str)))==NULL) {
      fprintf(stderr," cannot allocate pstat_union: %ld\n",sizeof(struct pstat_str));
      exit(1);
    }
    else {
      *ps_sp = ps_s;
    }
  }
  else {
    ps_s = *ps_sp;
    memset(ps_s,0,sizeof(struct pstat_str));
  }

  if (s_info->tot_scores > 10) {
    ps_s->sample_fract = min(1.0, (double)s_info->s_cnt[ppst->score_ix]/(double)s_info->tot_scores);
    if (ps_s->sample_fract > 0.0 && ps_s -> sample_fract < 1.0) {
      ps_s->zs_off = -log(ps_s->sample_fract)/PI_SQRT6;
    }
    else {
      ps_s->zs_off = 0.0;
    }
  }
  else {
    ps_s->sample_fract = 1.0;
    ps_s->zs_off = 0.0;
  }
  ppst->zs_off = ps_s->zs_off;

  ps_s->ngLambda = m_msp->Lambda;
  ps_s->ngK = m_msp->K;
  ps_s->ngH = m_msp->H;

  zsflag = ppst->zsflag;
  if (nstats < 10) zsflag = AG_STATS;

  if (m_msp->n0 <= LENGTH_CUTOFF && (zsflag == REG_STATS || zsflag == REGI_STATS || zsflag == REG2_STATS)) {
    zsflag = MLE_STATS;
  }

/*
#ifdef DEBUG
  if (ppst->debug_lib) {
    tmpf=fopen("tmp_stats.res","w+");
    for (i=0; i<nstats; i++) fprintf(tmpf,"%d\t%d\n",sptr[i].score,sptr[i].n1);
    fclose(tmpf);
  }
#endif
*/

#ifdef DEBUG
  st_sort(sptr, nstats);
#endif

  if (zsflag >= 20) {	/* working with shuffled top sequences */
    if (ppst->zsflag2 == AG_STATS) ppst->zsflag2 = MLE_STATS;
    zsflag = ppst->zsflag2;
    do_trim = 0;
  }
  else if (zsflag >= 10) {
    zsflag -= 10;
    do_trim = 0;
  }
  else do_trim = 1;

  if (!ppst->zdb_size_set) ppst->zdb_size = m_msp->db.entries;

#ifdef USE_LNSTATS
  if (zsflag==LN_STATS) {	/* -z 2 */
    find_z_arr[LN_STATS] = &find_zl;
    ps_s->zsflag=r_zsflag = proc_hist_ln(sptr, nstats, m_msp->s_info, histp, do_trim, ps_s);
  }
#else				/* -z 2 */
  if (zsflag==MLE_STATS) {
    find_z_arr[MLE_STATS] = &find_ze;
    ps_s->zsflag=r_zsflag = proc_hist_ml(sptr, nstats, m_msp->s_info, ppst, histp, do_trim, ps_s);
  }
#endif
  else if (zsflag==REG_STATS) {	/* -z 1 */
    ps_s->zsflag=r_zsflag = proc_hist_r(sptr, nstats, m_msp->s_info, ppst, histp, do_trim,  ps_s);
  }
  else if (zsflag==AG_STATS) {	/* -z 3 */
#ifndef NORMAL_DIST
    ps_s->zsflag=r_zsflag = proc_hist_a(sptr, nstats, m_msp->s_info, ppst, histp, do_trim,  ps_s);
#else
    ps_s->zsflag=r_zsflag = proc_hist_anorm(sptr, nstats, m_msp->s_info, ppst, histp, do_trim,  ps_s);
#endif

  }
  else if (zsflag==REGI_STATS) {	/* -z 4, iterated outlier removal */
    ps_s->zsflag=r_zsflag = proc_hist_ri(sptr,nstats, m_msp->s_info, ppst, histp, do_trim,  ps_s);
  }
  else if (zsflag==REG2_STATS) {	/* -z 5, length-scaled variance */
    ps_s->zsflag=r_zsflag = proc_hist_r(sptr,nstats, m_msp->s_info, ppst, histp, do_trim,  ps_s);
  }
#if !defined(TFAST) && !defined(FASTX)
  else if (zsflag == MLE2_STATS) {	/* -z 6, MLE + composition */
    ps_s->zsflag=r_zsflag = proc_hist_ml2(sptr, nstats, m_msp->s_info, ppst, histp, do_trim,  ps_s);
  }
#endif
  else {	/* AVE_STATS */
    ps_s->zsflag=r_zsflag = proc_hist_n(sptr,nstats, m_msp->s_info, ppst, histp, do_trim,  ps_s);
  }

  if (!do_hist && histp != NULL) {
    histp->entries = nstats; /* db->entries = 0; */
    inithistz(MAXHIST, histp, ps_s->zs_off);
    for (i = 0; i < nstats; i++) {
      if (sptr[i].n1 < 0) sptr[i].n1 = -sptr[i].n1;
      addhistz(find_z(sptr[i].score,sptr[i].escore,sptr[i].n1,sptr[i].comp,ps_s),
	       histp);
    }
  }
  return r_zsflag;
}

/* calc_thresh() sets the threshold used to exclude high scores from
   regression calculations
*/
static int
calc_thresh(struct pstruct *ppst, struct stat_str *sptr, int nstats,
	    struct score_count_s s_info,
	    double Lambda, double K, double H, double *zstrim)
{
  int max_hscore;
  int i;
  double ave_n1, tmp_score, z, l_fact;

  ave_n1 = 0.0;
  for (i=0; i<nstats; i++) {
    ave_n1 += (double)sptr[i].n1;
  }
  if (nstats >= 1) ave_n1 /= nstats;

  if (ppst->dnaseq == SEQT_DNA || ppst->dnaseq == SEQT_RNA) {
    l_fact = 1.0;
  }
  else {
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
  E()=0.01 and calculating with 0.70*lambda, which is the correct for
  going from ungapped to -12/-2 gapped lambda with BLOSUM50
*/

#ifndef NORMAL_DIST
  if (K < 1.0e-50) K = 0.02;	/* prevent floating pt underflow with K=0 */
  tmp_score = 0.01/((double)nstats*K*(double)ppst->n0*ave_n1);
  tmp_score = -log(tmp_score)/(Lambda*l_fact);
  max_hscore = (int)(tmp_score+0.5);

  /*   z = 1.0/(double)nstats; */
  z = (double)s_info.tot_scores/((double)nstats*(double)s_info.s_cnt[ppst->score_ix]);
  z = (log(z)+EULER_G)/(- PI_SQRT6);
  if (z < 3.0) { *zstrim = -1.0; }
  else { *zstrim = 10.0*z+50.0;}
#else
  max_hscore = 100;
  z = 5.0;
  *zstrim = 10.0*z+50.0;
#endif
  return max_hscore;
}

int
proc_hist_r(struct stat_str *sptr, int nstats, struct score_count_s s_info,
	    struct pstruct *ppst, struct hist_str *histp,
	    int do_trim, struct pstat_str *pu)
{
  int i, max_hscore;
  double zs, zstrim;
  char s_string[128];
  char rho_str[32];	/* used for simplifying output label */
  double rho_val;
  struct llen_str llen;
  char *f_string;
  llen.fit_flag=1;
  llen.hist=NULL;

  max_hscore = calc_thresh(ppst, sptr,  nstats, s_info, pu->ngLambda,
			   pu->ngK, pu->ngH, &zstrim);


  if (do_trim && zstrim < 0.0) {
    /* too few scores to do a z-trim, shift to Karlin-Altscul */
    return proc_hist_a(sptr, nstats, s_info, ppst, histp, do_trim, pu);
  }

  inithist(&llen, ppst,max_hscore);
  f_string = histp->stat_info;

  for (i = 0; i<nstats; i++)
    addhist(&llen,sptr[i].score,sptr[i].n1, max_hscore);

  if ((llen.max_score - llen.min_score) < 10) {
    free_hist(&llen);
    llen.fit_flag = 0;
    return proc_hist_n(sptr, nstats, s_info, ppst, histp, do_trim, pu);
  }

  fit_llen(&llen, &(pu->r_u.rg)); /* get parameters for REG_STATS */

  if (!llen.fit_flag) {	/* the fit failed, fall back to proc_hist_ml */
    free_hist(&llen);
    return proc_hist_ml(sptr,nstats, s_info, ppst, histp, do_trim, pu);
  }

  pu->r_u.rg.n_trimmed= pu->r_u.rg.n1_trimmed = pu->r_u.rg.nb_trimmed = 0;

  if (do_trim) {
    for (i = 0; i < nstats; i++) {
      zs = find_zr(sptr[i].score,sptr[i].escore,sptr[i].n1,sptr[i].comp, pu);
      if (zs < 20.0 || zs > zstrim) {
	pu->r_u.rg.n_trimmed++;
	prune_hist(&llen,sptr[i].score,sptr[i].n1, max_hscore,
		   &(histp->entries));
      }
    }

    if ((nstats - pu->r_u.rg.n_trimmed) < 50) {
      llen.fit_flag = 0;
      free_hist(&llen);
      return proc_hist_ml(sptr,nstats, s_info, ppst, histp, do_trim, pu);
    }

    /*  fprintf(stderr,"Z-trimmed %d entries with z > 5.0\n",
	pu->r_u.rg.n_trimmed); */

    /* calculate parameters for REG2_STATS; fit_llens() is always
       called, but the parameters are only used for REG2_STATS
       (find_zr2) */

    fit_llens(&llen, &(pu->r_u.rg));

    /*   fprintf(stderr,"Bin-trimmed %d entries in %d bins\n",
	 pu->r_u.rg.n1_trimmed,pu->r_u.rg.nb_trimmed); */
  }

  free_hist(&llen);

  /* put all the scores in the histogram */

  if (ppst->zsflag < 10) s_string[0]='\0';
  else if (ppst->zs_win > 0) {sprintf(s_string,"(shuffled [%d], win: %d)",nstats,ppst->zs_win);}
  else if (ppst->shuffle_dna3) { sprintf(s_string,"(shuffled3 [%d])",nstats);}
  else { sprintf(s_string,"(shuffled [%d])",nstats);}

  /* #ifdef LOCAL_SCORE */
  strncpy(rho_str,"ln(x)",sizeof(rho_str));
  rho_val = pu->r_u.rg.rho*LN_FACT;
  /*
#else
  strncpy(rho_str,"x",sizeof(rho_str));
  rho_val = pu->r_u.rg.rho;
#endif
  */

  if (ppst->zsflag == REG2_STATS || ppst->zsflag == 10+REG2_STATS ||
      ppst->zsflag == 20+REG2_STATS) {
    sprintf(f_string,"%s Expectation_v fit: rho(%s)= %6.4f+/-%6.3g; mu= %6.4f+/-%6.3f;\n rho2=%6.2f; mu2= %6.2f, 0's: %d Z-trim: %d  B-trim: %d in %d/%d",
	    s_string, rho_str, rho_val ,sqrt(pu->r_u.rg.rho_e),pu->r_u.rg.mu,sqrt(pu->r_u.rg.mu_e),
	    pu->r_u.rg.rho2,pu->r_u.rg.mu2,llen.zero_s,
	    pu->r_u.rg.n_trimmed, pu->r_u.rg.n1_trimmed, pu->r_u.rg.nb_trimmed, pu->r_u.rg.nb_tot);
  }
  else {
    sprintf(f_string,"%s Expectation_n fit: rho(%s)= %6.4f+/-%6.3g; mu= %6.4f+/-%6.3f\n mean_var=%6.4f+/-%6.3f, 0's: %d Z-trim(%.1f): %d  B-trim: %d in %d/%d\n Lambda= %8.6f",
	    s_string, rho_str, rho_val,sqrt(pu->r_u.rg.rho_e),pu->r_u.rg.mu,sqrt(pu->r_u.rg.mu_e), pu->r_u.rg.mean_var,sqrt(pu->r_u.rg.var_e),
	    llen.zero_s, zstrim, pu->r_u.rg.n_trimmed, pu->r_u.rg.n1_trimmed, pu->r_u.rg.nb_trimmed, pu->r_u.rg.nb_tot,
	    PI_SQRT6/sqrt(pu->r_u.rg.mean_var));
  }
  return REG_STATS;
}

/* proc_hist_ri() -- iterative removal of outliers before regression of mean vs log(len) */
int
proc_hist_ri(struct stat_str *sptr, int nstats, struct score_count_s s_info,
	     struct pstruct *ppst, struct hist_str *histp,
	     int do_trim, struct pstat_str *pu)
{
  int i, nit, nprune, max_hscore;
  double zs, zstrim;
  char s_string[128];
  char *f_string;
  struct llen_str llen;

  llen.fit_flag=1;
  llen.hist=NULL;

  max_hscore = calc_thresh(ppst, sptr, nstats, s_info, pu->ngLambda,
			   pu->ngK, pu->ngH, &zstrim);

  if (do_trim && zstrim < 0.0) {
    /* too few scores to do a z-trim, shift to Karlin-Altscul */
    return proc_hist_a(sptr, nstats, s_info, ppst, histp, do_trim, pu);
  }

  inithist(&llen, ppst,max_hscore);
  f_string = histp->stat_info;

  for (i = 0; i<nstats; i++)
    addhist(&llen,sptr[i].score,sptr[i].n1,max_hscore);

  pu->r_u.rg.n_trimmed= pu->r_u.rg.n1_trimmed = pu->r_u.rg.nb_trimmed = 0;
  if (do_trim) nit = 5;
  else nit = 0;

  while (nit-- > 0) {
    nprune = 0;
    fit_llen2(&llen, &(pu->r_u.rg));
    if (llen.fit_flag == 0) break;

    for (i = 0; i < nstats; i++) {
      if (sptr[i].n1 < 0) continue;
      zs = find_zr(sptr[i].score,sptr[i].escore,sptr[i].n1,sptr[i].comp,pu);
      if (zs < 20.0 || zs > zstrim ) {
	nprune++;
	pu->r_u.rg.n_trimmed++;
	prune_hist(&llen,sptr[i].score,sptr[i].n1,max_hscore,
		   &(histp->entries));
	sptr[i].n1 = -sptr[i].n1;
      }
    }
    /*    fprintf(stderr," %d Z-trimmed at %d\n",nprune,nit); */
    if (nprune < LHISTC) { break; }
  }

  fit_llen(&llen, &(pu->r_u.rg));
  free_hist(&llen);

  if (!llen.fit_flag) {	/* the fit failed, fall back to proc_hist_ml */
    return proc_hist_ml(sptr,nstats, s_info, ppst, histp, do_trim, pu);
  }

  if (ppst->zsflag < 10) s_string[0]='\0';
  else if (ppst->zs_win > 0)
    sprintf(s_string,"(shuffled [%d], win: %d)",nstats,ppst->zs_win);
  else 
    sprintf(s_string,"(shuffled [%d])",nstats);

  sprintf(f_string,"%s Expectation_i fit: rho(ln(x))= %6.4f+/-%6.3g; mu= %6.4f+/-%6.3f;\n mean_var=%6.4f+/-%6.3f 0's: %d Z-trim: %d N-it: %d\n Lambda= %8.6f",
	  s_string,
	  pu->r_u.rg.rho*LN_FACT,sqrt(pu->r_u.rg.rho_e),pu->r_u.rg.mu,sqrt(pu->r_u.rg.mu_e),
	  pu->r_u.rg.mean_var,sqrt(pu->r_u.rg.var_e),llen.zero_s,pu->r_u.rg.n_trimmed, nit,
	  PI_SQRT6/sqrt(pu->r_u.rg.mean_var));
  return REGI_STATS;
}

/* this procedure implements Altschul's pre-calculated values for lambda, K */

#include "alt_parms.h"

int
look_p(struct alt_p parm[], int gap, int ext, 
       double *K, double *Lambda, double *H);

int
proc_hist_a(struct stat_str *sptr, int nstats, struct score_count_s s_info,
	    struct pstruct *ppst, struct hist_str *histp,
	    int do_trim, struct pstat_str *pu)
{
  double Lambda, K, H;
  char *f_string;
  int r_v;
  int t_gdelval, t_ggapval;

#ifdef OLD_FASTA_GAP
  t_gdelval = ppst->gdelval;
  t_ggapval = ppst->ggapval;
#else
  t_gdelval = ppst->gdelval+ppst->ggapval;
  t_ggapval = ppst->ggapval;
#endif

  f_string = histp->stat_info;

  if (ppst->dnaseq==0) {
    if (strcmp(ppst->pam_name,"BL50")==0 || strcmp(ppst->pam_name,"BLOSUM50")==0) {
      r_v = look_p(bl50_p,t_gdelval,t_ggapval,&K,&Lambda,&H);
    }
    else if (strcmp(ppst->pam_name,"BL62")==0 || strcmp(ppst->pam_name,"BLOSUM62")==0) {
      r_v = look_p(bl62_p,t_gdelval,t_ggapval,&K,&Lambda,&H);
    }
    else if (strcmp(ppst->pam_name,"BL80")==0 || strcmp(ppst->pam_name,"BLOSUM80")==0) {
      r_v = look_p(bl80_p,t_gdelval,t_ggapval,&K,&Lambda,&H);
    }
    else if (strcmp(ppst->pam_name,"PAM250")==0) {
      r_v = look_p(p250_p,t_gdelval,t_ggapval,&K,&Lambda,&H);
    }
    else if ((strcmp(ppst->pam_name,"PAM120")==0) || (strcmp(ppst->pam_name,"VT120")==0)) {
      r_v = look_p(p120_p,t_gdelval,t_ggapval,&K,&Lambda,&H);
    }
    else if ((strcmp(ppst->pam_name,"MD10")==0) || (strcmp(ppst->pam_name,"VT10")==0)) {
      r_v = look_p(md10_p,t_gdelval,t_ggapval,&K,&Lambda,&H);
    }
    else if ((strcmp(ppst->pam_name,"MD20")==0) || (strcmp(ppst->pam_name,"VT20")==0)) {
      r_v = look_p(md20_p,t_gdelval,t_ggapval,&K,&Lambda,&H);
    }
    else if ((strcmp(ppst->pam_name,"MD40")==0) || (strcmp(ppst->pam_name,"VT40")==0)) {
      r_v = look_p(md40_p,t_gdelval,t_ggapval,&K,&Lambda,&H);
    }
    else if (strcmp(ppst->pam_name,"OPTIMA5")==0) {
      r_v = look_p(opt5_p,t_gdelval,t_ggapval,&K,&Lambda,&H);
    }
    else r_v = 0;
  }
  else {
    r_v = look_p(nt32_p,t_gdelval,t_ggapval,&K,&Lambda,&H);
    if (strcmp(ppst->pam_name,"DNA")==0 || (ppst->pam_h==5 && ppst->pam_l == -4)) {
      r_v = look_p(nt54_p,t_gdelval,t_ggapval,&K,&Lambda,&H);
    }
    else if (strcmp(ppst->pam_name,"+3/-2")==0 || (ppst->pam_h==3 && ppst->pam_l == -2)) {
      r_v = look_p(nt32_p,t_gdelval,t_ggapval,&K,&Lambda,&H);
    }
    else if (strcmp(ppst->pam_name,"+1/-3")==0 || (ppst->pam_h==1 && ppst->pam_l == -3)) {
      r_v = look_p(nt13_p,t_gdelval,t_ggapval,&K,&Lambda,&H);
    }
  }

  if (r_v != 1) {
      r_v = look_p(bl62_p,t_gdelval,t_ggapval,&K,&Lambda,&H);
#ifdef DEBUG
    fprintf(stderr,"*** Warning : [%s:%d] Parameters not available for: %s: %d/%d -- using BL62\n",
	    __FILE__, __LINE__, ppst->pam_name,t_gdelval-t_ggapval,t_ggapval);
#endif
  }
  pu->r_u.ag.Lambda = Lambda;
  pu->r_u.ag.K = K;
  pu->r_u.ag.H = H;

  /*
    fprintf(stderr," the parameters are: Lambda: %5.3f K: %5.3f H: %5.3f\n",
	    Lambda, K, H);
	    */

    pu->r_u.ag.a_n0 = (double)ppst->n0;
    pu->r_u.ag.a_n0f = log (K * pu->r_u.ag.a_n0)/H;

    sprintf(f_string,"Altschul/Gish params: n0: %d Lambda: %5.3f K: %5.3f H: %5.3f",
	    ppst->n0,Lambda, K, H);
    return AG_STATS;
}

int 
ag_parm(char *pam_name, int gdelval, int ggapval, struct pstat_str *pu)
{
  double Lambda, K, H;
  int r_v;

  if (strcmp(pam_name,"BL50")==0)
    r_v = look_p(bl50_p,gdelval,ggapval,&K,&Lambda,&H);
  else if (strcmp(pam_name,"BL62")==0)
      r_v = look_p(bl62_p,gdelval,ggapval,&K,&Lambda,&H);
  else if (strcmp(pam_name,"P250")==0)
      r_v = look_p(p250_p,gdelval,ggapval,&K,&Lambda,&H);
  else if (strcmp(pam_name,"P120")==0)
      r_v = look_p(p120_p,gdelval,ggapval,&K,&Lambda,&H);
  else if (strcmp(pam_name,"MD10")==0 || strcmp(pam_name,"VT10")==0)
      r_v = look_p(md10_p,gdelval,ggapval,&K,&Lambda,&H);
  else if (strcmp(pam_name,"MD20")==0 || strcmp(pam_name,"VT20")==0)
      r_v = look_p(md20_p,gdelval,ggapval,&K,&Lambda,&H);
  else if (strcmp(pam_name,"MD40")==0 || strcmp(pam_name,"VT40")==0)
      r_v = look_p(md40_p,gdelval,ggapval,&K,&Lambda,&H);
  else if (strcmp(pam_name,"DNA")==0 || strcmp(pam_name,"+5/-4")==0)
      r_v = look_p(nt54_p,gdelval,ggapval, &K,&Lambda,&H);
  else if (strcmp(pam_name,"+3/-2")==0)
      r_v = look_p(nt32_p,gdelval,ggapval, &K,&Lambda,&H);
  else if (strcmp(pam_name,"+1/-3")==0)
      r_v = look_p(nt13_p,gdelval,ggapval, &K,&Lambda,&H);
  else r_v = 0;

  pu->r_u.ag.K = K;
  pu->r_u.ag.Lambda = Lambda;
  pu->r_u.ag.H = H;

  /*
  if (r_v == 0) {
    fprintf(stderr,"Parameters not available for: %s: %d/%d\n",
	    pam_name,gdelval,ggapval);
  }
  */
  return r_v;
}

int
look_p(struct alt_p parm[], int gap, int ext,
       double *K, double *Lambda, double *H)
{
  int i;

  gap = -gap;
  ext = -ext;

  if (gap > parm[1].gap) {
    *K = parm[0].K;
    *Lambda = parm[0].Lambda;
    *H = parm[0].H;
    return 1;
  }

  for (i=1; parm[i].gap > 0; i++) {
    if (parm[i].gap > gap) continue;
    else if (parm[i].gap <= gap && parm[i].ext > ext ) continue;
    else if (parm[i].gap <= gap && parm[i].ext <= ext) {
      *K = parm[i].K;
      *Lambda = parm[i].Lambda;
      *H = parm[i].H;
      return 1;
    }
    else break;
  }
  return 0;
}

#ifdef NORMAL_DIST
int
proc_hist_anorm(struct stat_str *sptr, int nstats, struct score_count_s s_info,
		struct pstruct *ppst, struct hist_str *histp,
		int do_trim, struct pstat_str *pu)
{
  char s_string[128];
  char *f_string;
  int med_ix, q1_ix, q3_ix;
  double iqr_sd;

  f_string = histp->stat_info;

  if (ppst->zsflag < 10) s_string[0]='\0';
  else if (ppst->zs_win > 0) {
    sprintf(s_string,"(shuffled [%d], win: %d)",nstats,ppst->zs_win);
  }
  else {
    sprintf(s_string,"(shuffled [%d])",nstats);
  }

  /* find median, 1st/3rd quartile */
  
  med_ix = nstats/2;
  q3_ix = med_ix - med_ix/2;
  q1_ix = med_ix + med_ix/2;

  st_sort(sptr,nstats);

  pu->r_u.rg.mu = sptr[med_ix].score;

  iqr_sd = (sptr[q1_ix].score - sptr[q3_ix].score)/1.35;

  pu->r_u.rg.mean_var_sqrt = iqr_sd;
  pu->r_u.rg.mean_var = iqr_sd * iqr_sd;
  pu->r_u.rg.n_trimmed = 0;

  sprintf(f_string,
	  "%s Unscaled normal (median, IQR) statistics: mu= %6.4f  var=%6.4f Ztrim: %d",
	  s_string, pu->r_u.rg.mu,pu->r_u.rg.mean_var, pu->r_u.rg.n_trimmed
	  );
  return AG_STATS;
}
#endif

/* uncensored and censored maximum likelihood estimates developed
   by Aaron Mackey based on a preprint from Sean Eddy */

int mle_cen  (struct stat_str *, int, int, double, double *, double *);

int
proc_hist_ml(struct stat_str *sptr, int nstats, struct score_count_s s_info,
	     struct pstruct *ppst, struct hist_str *histp,
	     int do_trim, struct pstat_str *pu)
{
  double f_cen;
  char s_string[128];
  char *f_string;
  double opt_correct;

  f_string = histp->stat_info;
  pu->r_u.ag.a_n0 = (double)ppst->n0;

  if (ppst->zsflag < 10) s_string[0]='\0';
  else {	/* shuffled, do_trim = 0 */
    if (ppst->zs_win > 0) {
      sprintf(s_string,"(shuffled [%d], win: %d)",nstats,ppst->zs_win);
    }
    else {
      sprintf(s_string,"(shuffled [%d])",nstats);
    }
  }

  if (!do_trim) {
    if (mle_cen(sptr, nstats, ppst->n0, 0.0, &pu->r_u.ag.Lambda, &pu->r_u.ag.K) == -1)
      goto bad_mle;
    sprintf(f_string,"%s MLE statistics: Lambda= %6.4f;  K=%6.4g",
	    s_string,pu->r_u.ag.Lambda,pu->r_u.ag.K);
  }
  else {
    /* this section attempts to censor the appropriate fraction of
       high-scores.  In general, we would like 5% censored.
       Unfortunately, for fasta36 with statistical thresholds for
       optimization, the thresholds mean that only a small fraction of
       the highest scoring library sequences are optimized (and used
       for scores).  As a result, there may be no low scores, and a
       large excess of high-scores.  Thus, when only the best scores
       are calculated initially, more censoring must be done.

       The pre-selection information is available in struct
       score_count_s m_msp->s_info, by comparing
       s_info.s_cnt[ppst->score_ix] to sinfo.tot_scores.
     */

    opt_correct = (double)s_info.s_cnt[ppst->score_ix]/(double)s_info.tot_scores;

    if ((double)nstats/20.0 > 1000.0) f_cen = 1000.0/((double)nstats*opt_correct);
    else f_cen = 0.05/opt_correct;

    if (f_cen > 0.4) f_cen = 0.2;

    if (mle_cen(sptr, nstats, ppst->n0, f_cen, &pu->r_u.ag.Lambda, &pu->r_u.ag.K) == -1)
      goto bad_mle;
    sprintf(f_string, "MLE_cen statistics: Lambda= %6.4f;  K=%6.4g (cen=%d)",
	    pu->r_u.ag.Lambda,pu->r_u.ag.K,(int)((double)nstats*f_cen));
  }    

  return MLE_STATS;
 bad_mle:
  return proc_hist_n(sptr, nstats, s_info, ppst, histp, do_trim, pu);
}

int
mle_cen2  (struct stat_str *, int, int, double, double *, double *, double *, double *);


int
proc_hist_ml2(struct stat_str *sptr, int nstats, struct score_count_s s_info,
	      struct pstruct *ppst, struct hist_str *histp,
	      int do_trim, struct pstat_str *pu)
{
  int i, ns=0, nneg=0;
  double f_cen, ave_lambda;
  char s_string[128], ex_string[64];
  char *f_string;

  f_string = histp->stat_info;
  pu->r_u.m2.a_n0 = (double)ppst->n0;

  if (ppst->zsflag < 10) s_string[0]='\0';
  else if (ppst->zs_win > 0)
    sprintf(s_string,"(shuffled [%d], win: %d)",nstats,ppst->zs_win);
  else 
    sprintf(s_string,"(shuffled [%d])",nstats);


  pu->r_u.m2.ave_comp = 0.0;
  pu->r_u.m2.max_comp = -1.0;

  ns = nneg = 0;
  for (i=0; i<nstats; i++) {
    if (sptr[i].comp > pu->r_u.m2.max_comp) pu->r_u.m2.max_comp = sptr[i].comp;
    if (sptr[i].comp > 0.0) {
      pu->r_u.m2.ave_comp += log(sptr[i].comp);
      ns++;
    }
    else nneg++;
  }
  pu->r_u.m2.ave_comp /= (double)ns;
  pu->r_u.m2.ave_comp = exp(pu->r_u.m2.ave_comp);
  for (i=0; i<nstats; i++) if (sptr[i].comp < 0.0) {
    sptr[i].comp = pu->r_u.m2.ave_comp;
  }

  if (nneg > 0)
    sprintf(ex_string,"composition = -1 for %d sequences",nneg);
  else ex_string[0]='\0';

  if (!do_trim) {
    if (mle_cen2(sptr, nstats, ppst->n0, 0.0,
	     &pu->r_u.m2.mle2_a0, &pu->r_u.m2.mle2_a1,
	     &pu->r_u.m2.mle2_a2, &pu->r_u.m2.mle2_b1) == -1) goto bad_mle2;
    ave_lambda = 1.0/(pu->r_u.m2.ave_comp*pu->r_u.m2.mle2_b1);

    sprintf(f_string,"%s MLE-2 statistics: a0= %6.4f;  a1=%6.4f; a2=%6.4f; b1=%6.4f\n  ave Lamdba: %6.4f",
	    s_string, pu->r_u.m2.mle2_a0, pu->r_u.m2.mle2_a1, pu->r_u.m2.mle2_a2, pu->r_u.m2.mle2_b1,ave_lambda);
  }
  else {
    if (nstats/20 > 500) f_cen = 500.0/(double)nstats;
    else f_cen = 0.05;
    if (mle_cen2(sptr, nstats, ppst->n0, f_cen, &pu->r_u.m2.mle2_a0, &pu->r_u.m2.mle2_a1, &pu->r_u.m2.mle2_a2, &pu->r_u.m2.mle2_b1)== -1) goto bad_mle2;

    ave_lambda = 1.0/(pu->r_u.m2.ave_comp*pu->r_u.m2.mle2_b1);

    sprintf(f_string,"%s MLE-2-cen statistics: a0= %6.4f;  a1=%6.4f; a2=%6.4f; b1=%6.4f (cen=%d)\n  ave Lambda:%6.4f",
	    s_string, pu->r_u.m2.mle2_a0, pu->r_u.m2.mle2_a1, pu->r_u.m2.mle2_a2, pu->r_u.m2.mle2_b1, (int)((double)nstats*f_cen),ave_lambda);
  }    

  return MLE2_STATS;
 bad_mle2:
  return proc_hist_n(sptr, nstats, s_info, ppst, histp, do_trim, pu);
}

double first_deriv_cen(double lambda, struct stat_str *sptr, 
		       int start, int stop,
		       double sumlenL, double cenL,
		       double sumlenH, double cenH);

double second_deriv_cen(double lambda, struct stat_str *sptr,
			int start, int stop,
			double sumlenL, double cenL,
			double sumlenH, double cenH);

static void
st_sort (struct stat_str *v, int n) {
  int gap, i, j;
  int tmp;
  double dtmp;

  for (gap = 1; gap < n/3; gap = 3*gap +1) ;

  for (; gap > 0; gap = (gap-1)/3) {
    for (i = gap; i < n; i++) {
      for (j = i - gap; j >= 0; j -= gap) {
	if (v[j].score <= v[j + gap].score) break;

	tmp = v[j].score;
	v[j].score = v[j + gap].score;
	v[j + gap].score = tmp;

	tmp = v[j].n1;
	v[j].n1 = v[j + gap].n1;
	v[j + gap].n1 = tmp;

	dtmp = v[j].comp;
	v[j].comp = v[j + gap].comp;
	v[j + gap].comp = dtmp;

	dtmp = v[j].H;
	v[j].H = v[j + gap].H;
	v[j + gap].H = dtmp;
      }
    }
  }
}


/* sptr[].score, sptr[].n1; sptr[] must be sorted
   int n = total number of samples
   int M = length of query
   double fn = fraction of scores to be censored fn/2.0 from top, bottom
   double *Lambda = Lambda estimate
   double *K = K estimate
*/

#define MAX_NIT 100

int
mle_cen(struct stat_str *sptr, int n, int M, double fc,
	double *Lambda, double *K) {

  double sumlenL, sumlenH, cenL, cenH;
  double sum_s, sum2_s, mean_s, var_s, dtmp;
  int start, stop;
  int i, nf;
  int nit = 0;
  double deriv, deriv2, lambda, old_lambda, sum = 0.0;

  /*
   int sumlenL, int sumlenghtsR = sum of low (Left), right (High) seqs.
   int cenL, cenH = censoring score low, high 
  */

  nf = (fc/2.0) * n;
  start = nf;
  stop = n - nf;

  st_sort(sptr,n);

  sum_s = sum2_s = 0.0;
  for (i=start; i<stop; i++) {
    sum_s += sptr[i].score;
  }
  dtmp = (double)(stop-start);
  mean_s = sum_s/dtmp;

  for (i=start; i<stop; i++) {
    sum2_s += sptr[i].score * sptr[i].score;
  }
  var_s = sum2_s/(dtmp-1.0);

  sumlenL = sumlenH = 0.0;
  for (i=0; i<start; i++) sumlenL += (double)sptr[i].n1;
  for (i=stop; i<n; i++) sumlenH += (double)sptr[i].n1;

  if (nf > 0) {
    cenL = (double)sptr[start].score;
    cenH = (double)sptr[stop].score;
  }
  else {
    cenL = (double)sptr[start].score/2.0;
    cenH = (double)sptr[start].score*2.0;
  }

  if (cenL >= cenH) return -1;

  /* initial guess for lambda is 0.2 - this does not work for matrices
     with very different scales */
  /*  lambda = 0.2; */
  lambda = PI_SQRT6/sqrt(var_s);
  if (lambda > 1.0) {
    fprintf(stderr,"*** Warning [%s:%d] - Lambda initial estimate error: lambda: %6.4g; var_s: %6.4g\n",
	    __FILE__, __LINE__, lambda,var_s);
    lambda = 0.2;
  }

  do {
    deriv = first_deriv_cen(lambda, sptr, start, stop,
			    sumlenL, cenL, sumlenH, cenH);
    /*   (uncensored version)
	 first_deriv(lambda, &sptr[start], stop - start))
    */

    /*  (uncensored version)
    deriv2 = second_deriv(lambda, &sptr[start], stop-start);
    */
    deriv2 = second_deriv_cen(lambda, sptr, start, stop,
			     sumlenL, cenL, sumlenH, cenH); 

    old_lambda = lambda;
    if (lambda - deriv/deriv2 > 0.0) lambda = lambda - deriv/deriv2;
    else lambda = lambda/2.0;
    nit++;
  } while (fabs((lambda - old_lambda)/lambda) > TINY && nit < MAX_NIT);

  /*  fprintf(stderr," mle_cen nit: %d\n",nit); */

  if (nit >= MAX_NIT) return -1;
  
  for(i = start; i < stop ; i++) {
    sum += (double) sptr[i].n1 * exp(- lambda * (double)sptr[i].score);
  }

  *Lambda = lambda;
  /* 
  *K = (double)(stop-start)/((double)M*sum);
  */
  *K = (double)n/((double)M*
		  (sum+sumlenL*exp(-lambda*cenL)-sumlenH*exp(-lambda*cenH)));
  return 0;
}

/*
double
first_deriv(double lambda, struct stat_str *sptr, int n) {

  int i;
  double sum = 0.0, sum1 = 0.0, sum2 = 0.0;
  double s, l, es;

  for(i = 0 ; i < n ; i++) {
    s = (double)sptr[i].score;
    l = (double)sptr[i].n1;
    es = exp(-lambda * s );
    sum += s;
    sum2 += l * es;
    sum1 += s * l * es;
  }

  return (1.0/lambda) - (sum/(double)n) + (sum1/sum2);
}
*/

/*
double
second_deriv(double lambda, struct stat_str *sptr, int n) {
  double sum1 = 0.0, sum2 = 0.0, sum3 = 0.0;
  double s, l, es;
  int i;

  for(i = 0 ; i < n ; i++) {
    l = (double)sptr[i].n1;
    s = (double)sptr[i].score;
    es = exp(-lambda * s);
    sum2 += l * es;
    sum1 += l * s * es;
    sum3 += l * s * s * es;
  }

  return ((sum1*sum1)/(sum2*sum2)) - (sum3/sum2) -  (1.0/(lambda*lambda));
}
*/

double
first_deriv_cen(double lambda, struct stat_str *sptr, int start, int stop,
		double sumlenL, double cenL, double sumlenH, double cenH) {
  int i;
  double sum = 0.0, sum1 = 0.0, sum2 = 0.0;
  double s, l, es;

  for(i = start ; i < stop ; i++) {
    s = (double)sptr[i].score;
    l = (double)sptr[i].n1;
    es = exp(-lambda * s );
    sum += s;
    sum2 += l * es;
    sum1 += s * l * es;
  }

  sum1 += sumlenL*cenL*exp(-lambda*cenL) - sumlenH*cenH*exp(-lambda*cenH);
  sum2 += sumlenL*exp(-lambda*cenL) - sumlenH*exp(-lambda*cenH);

  return (1.0 / lambda) - (sum /(double)(stop-start)) + (sum1 / sum2);
}

double
second_deriv_cen(double lambda, struct stat_str *sptr, int start, int stop,
		 double sumlenL, double cenL, double sumlenH, double cenH) {

  double sum1 = 0.0, sum2 = 0.0, sum3 = 0.0;
  double s, l, es;
  int i;

  for(i = start ; i < stop ; i++) {
    s = (double)sptr[i].score;
    l = (double)sptr[i].n1;
    es = exp(-lambda * s);
    sum2 += l * es;
    sum1 += l * s * es;
    sum3 += l * s * s * es;
  }

  sum1 += sumlenL*cenL*exp(-lambda*cenL) - sumlenH*cenH*exp(-lambda*cenH);
  sum2 += sumlenL*exp(-lambda * cenL) -  sumlenH*exp(-lambda * cenH);
  sum3 += sumlenL*cenL*cenL * exp(-lambda * cenL) -
    sumlenH*cenH*cenH * exp(-lambda * cenH);
  return ((sum1 * sum1) / (sum2 * sum2)) - (sum3 / sum2)
    - (1.0 / (lambda * lambda));
}

double mle2_func(double *params,
		 double *consts,
		 struct stat_str *values,
		 int n, int start, int stop);

void simplex(double *fitparams,
	     double *lambda,
	     int nparam,
	     double (*minfunc) (double *tryparams, double *consts, 
				struct stat_str *data, int ndata,
				int start, int stop),
	     double *consts,
	     void *data,
	     int ndata, int start, int stop
	     );

int
mle_cen2(struct stat_str *sptr, int n, int M, double fc,
	double *a0, double *a1, double *a2, double *b1) {

  double params[4], lambdas[4], consts[9];
  double avglenL, avglenH, avgcompL, avgcompH, cenL, cenH;
  int start, stop;
  int i, nf;

  nf = (fc/2.0) * n;
  start = nf;
  stop = n - nf;

  st_sort(sptr,n);

  /* choose arithmetic or geometic mean for compositions by appropriate commenting */

  if (nf > 0) {
    avglenL = avglenH = 0.0;
    avgcompL = avgcompH = 0.0;
    /* avgcompL = avgcompH = 1.0 */
    for (i=0; i<start; i++) {
      avglenL += (double)sptr[i].n1;
      avgcompL += (double)sptr[i].comp;
      /* avgcompL *= (double) sptr[i].comp; */
    }
    avglenL /= (double) start;
    avgcompL /= (double) start;
    /* avgcompL = pow(avgcompL, 1.0/(double) start); */
  
    for (i=stop; i<n; i++) {
      avglenH += (double)sptr[i].n1;
      avgcompH += (double)sptr[i].comp;
      /* avgcompH *= (double) sptr[i].comp; */
    }
    avglenH /= (double) (n - stop);
    avgcompH /= (double) (n - stop);
    /* avgcompL = pow(avgcompL, 1.0/(double) (n - stop)); */

    cenL = (double)sptr[start].score;
    cenH = (double)sptr[stop].score;
    if (cenL >= cenH) return -1;
  }
  else {
    avglenL = avglenH = cenL = cenH = 0.0;
    avgcompL = avgcompH = 1.0;
  }

  params[0] = 10.0;
  params[1] = -10.0;
  params[2] = 1.0;
  params[3] = 1.0;

  lambdas[0] = 1.0;
  lambdas[1] = 0.5;
  lambdas[2] = 0.1;
  lambdas[3] = 0.01;

  consts[0] = M;
  consts[1] = (double) start;
  consts[2] = (double) stop;
  consts[3] = cenL;
  consts[4] = cenH;
  consts[5] = avglenL;
  consts[6] = avglenH;
  consts[7] = avgcompL;
  consts[8] = avgcompH;

  simplex(params, lambdas, 4,
	  (double (*) (double *, double *, struct stat_str *, int, int, int) )mle2_func,
	  consts, sptr, n, start, stop);

  *a0 = params[0];
  *a1 = params[1];
  *a2 = params[2];
  *b1 = params[3];

  return 0;
}

double mle2_func(double *params,
		 double *consts,
		 struct stat_str *values,
		 int n, int start, int stop
		 ) {

  double a0, a1, a2, b1, M;
  double score, length, comp;
  double cenL, cenH, avglenL, avglenH, avgcompL, avgcompH;
  double L, y;

  int i;

  a0 = params[0];
  a1 = params[1];
  a2 = params[2];
  b1 = params[3];

  M = consts[0];
  /*
  start = (int) consts[1];
  stop = (int) consts[2];
  */
  cenL = consts[3];
  cenH = consts[4];
  avglenL = consts[5];
  avglenH = consts[6];
  avgcompL = consts[7];
  avgcompH = consts[8];

  L = 0;
  y = 0;

  if (start > 0) {
    y = -(cenL - (a0 + a1*avgcompL +a2*avgcompL*log(M*avglenL)))/(b1*avgcompL);
    L += (double) start * exp(y);
  }

  for(i = start ; i < stop ; i++) {
    score = (double) values[i].score;
    length = (double) values[i].n1;
    comp = (double) values[i].comp;

    y = - (score - (a0 + a1*comp + a2 * comp * log(M*length))) / (b1*comp);

    L += -y + exp(y) + log(b1 * comp);
  }

  if (stop < n) {
    y = -(cenH -(a0 + a1*avgcompH + a2*avgcompH*log(M*avglenH)))/(b1*avgcompH);
    L -= (double) (n - stop) * exp(y);
  }
  return L;
}

/* Begin Nelder-Mead simplex code: */

double evalfunc(double **param,
		double *vals,
		double *psums,
		double *ptry,
		int nparam,
		double (*minfunc) (double *params, double *consts,
				   struct stat_str *data, int ndata,
				   int start, int stop),
		double *consts,
		void *data,
		int ndata, int start, int stop,
		int ihi,
		double factor);

void simplex(double *fitparams,
	     double *lambda,
	     int nparam,
	     double (*minfunc) (double *tryparams, double *consts,
				struct stat_str *data, int ndata, 
				int start, int stop),
	     double *consts,
	     void *data,
	     int ndata,
	     int start,
	     int stop
	     )
{

  int i, j, ilo, ihi, inhi;
  double rtol, sum, tmp, ysave, ytry;
  double *psum, *vals, *ptry, **param;


  psum = (double *) calloc(nparam, sizeof(double));
  ptry = (double *) calloc(nparam, sizeof(double));

  vals = (double *) calloc(nparam + 1, sizeof(double));

  param = (double **) calloc(nparam + 1, sizeof(double *));
  param[0] = (double *) calloc((nparam + 1) * nparam, sizeof(double));
  for( i = 1 ; i < (nparam + 1) ; i++ ) {
    param[i] = param[0] + i * nparam;
  }

  /* Get our N+1 initial parameter values for the simplex */

  for( i = 0 ; i < nparam ; i++ ) {
    param[0][i] = fitparams[i];
  }

  for( i = 1 ; i < (nparam + 1) ; i++ ) {
    for( j = 0 ; j < nparam ; j++ ) {
      param[i][j] = fitparams[j] + lambda[j] * ( (i - 1) == j ? 1 : 0 );
    }
  }

  /* calculate initial values at the simplex nodes */

  for( i = 0 ; i < (nparam + 1) ; i++ ) {
    vals[i] = minfunc(param[i], consts, data, ndata, start, stop);
  }

  /* Begin Nelder-Mead simplex algorithm from Numerical Recipes in C */

  for( j = 0 ; j < nparam ; j++ ) {
    for( sum = 0.0, i = 0 ; i < nparam + 1 ; i++ ) {
      sum += param[i][j];
    }
    psum[j] = sum;
  }


  while( 1 ) {
/*
      determine which point is highest (ihi), next highest (inhi) and
      lowest (ilo) by looping over the points in the simplex
*/
    ilo = 0;

/*  ihi = vals[0] > vals[1] ? (inhi = 1, 0) : (inhi = 0, 1); */
    if(vals[0] > vals[1]) { ihi = 0; inhi = 1; }
    else { ihi = 1; inhi = 0; }

    for( i = 0 ; i < nparam + 1 ; i++) {
      if( vals[i] <= vals[ilo] ) ilo = i;
      if( vals[i] > vals[ihi] ) {
	inhi = ihi;
	ihi = i;
      } else if ( vals[i] > vals[inhi] && i != ihi ) inhi = i;
    }

    /* Are we finished? */

    rtol = 2.0 * fabs(vals[ihi] - vals[ilo]) / 
      (fabs(vals[ihi]) + fabs(vals[ilo]) + TINY);

    if( rtol < TOLERANCE ) {

/* put the best value and best parameters into the first index */

      tmp = vals[0];
      vals[0] = vals[ilo];
      vals[ilo] = tmp;

      for( i = 0 ; i < nparam ; i++ ) {
	tmp = param[0][i];
	param[0][i] = param[ilo][i];
	param[ilo][i] = tmp;
      }

      /* et voila, c'est finis */
      break;
    }

    /* Begin a new iteration */

    /* first, extrapolate by -1 through the face of the simplex across from ihi */

    ytry = evalfunc(param, vals, psum, ptry, nparam, minfunc, consts,
		    data, ndata, start, stop, ihi, -1.0);

    if( ytry <= vals[ilo] ) {

      /* Good result, try additional extrapolation by 2 */

      ytry = evalfunc(param, vals, psum, ptry, nparam, minfunc, consts, 
		      data, ndata, start, stop, ihi, 2.0);

    } else if ( ytry >= vals[inhi] ) {

      /* no good, look for an intermediate lower point by contracting */

      ysave = vals[ihi];
      ytry = evalfunc(param, vals, psum, ptry, nparam, minfunc, consts,
		      data, ndata, start, stop, ihi, 0.5);

      if( ytry >= ysave ) {

	/* Still no good.  Contract around lowest (best) point. */

	for( i = 0 ; i < nparam + 1 ; i++ ) {
	  if( i != ilo ) {
	    for ( j = 0 ; j < nparam ; j++ ) {
	      param[i][j] = psum[j] = 0.5 * (param[i][j] + param[ilo][j]);
	    }
	    vals[i] = minfunc(psum, consts, data, ndata, start, stop);
	  }
	}


	for( j = 0 ; j < nparam ; j++ ) {
	  for( sum = 0.0, i = 0 ; i < nparam + 1 ; i++ ) {
	    sum += param[i][j];
	  }
	  psum[j] = sum;
	}

      }
    }
  }
			   
  for( i = 0 ; i < nparam ; i++ ) {
    fitparams[i] = param[0][i];
  }

  if (ptry!=NULL) {
    free(ptry);
    ptry=NULL;
  }
  free(param[0]);
  free(param);
  free(vals);
  free(psum);
}


double evalfunc(double **param,
		double *vals,
		double *psum,
		double *ptry,
		int nparam,
		double (*minfunc)(double *tryparam, double *consts,
				  struct stat_str *data, int ndata,
				  int start, int stop),
		double *consts,
		void *data,
		int ndata, int start, int stop,
		int ihi,
		double factor) {

  int j;
  double fac1, fac2, ytry;


  fac1 = (1.0 - factor) / nparam;
  fac2 = fac1 - factor;

  for( j = 0 ; j < nparam ; j++ ) {
    ptry[j] = psum[j] * fac1 - param[ihi][j] * fac2;
  }

  ytry = minfunc(ptry, consts, data, ndata, start, stop);

  if( ytry < vals[ihi] ) {
    vals[ihi] = ytry;
    for( j = 0 ; j < nparam ; j++ ) {
      psum[j] += ptry[j] - param[ihi][j];
      param[ihi][j] = ptry[j];
    }
  }

  return ytry;
}

/* end of Nelder-Mead simplex code */

int
proc_hist_n(struct stat_str *sptr, int nstats, struct score_count_s s_info,
	    struct pstruct *ppst, struct hist_str *histp,
	    int do_trim, struct pstat_str *pu)
{
  int i, j, t_j;
  double s_score, s2_score, ssd, zstrim, mean, var;
  double t_score;
  int max_trim, min_trim;
  int nit, max_hscore;
  char s_string[128];
  char *f_string;
  int no_trim = 0;

  f_string = histp->stat_info;

  max_hscore = calc_thresh(ppst, sptr, nstats, s_info, pu->ngLambda,
			   pu->ngK, pu->ngH, &zstrim);

  if (do_trim && zstrim < 0.0) {
    /* too few scores to do a z-trim, shift to Karlin-Altscul */
    /* cannot use proc_hist_a(), because proc_hist_a() calls proc_hist_n() */
    /* return proc_hist_a(sptr, nstats, s_info, ppst, histp, do_trim, pu); */
    no_trim = 1;
  }

  s_score = s2_score = 0.0;

  /* calculate mean */
  for (j=0, i=0; i < nstats; i++) {
    if ( sptr[i].score <= max_hscore) {
      s_score += (double)sptr[i].score;
      j++;
    }
    else { sptr[i].n1 = -sptr[i].n1;}
  }

  if (j > 1) {
    pu->r_u.rg.mu = mean = s_score/(double)j;

    /* calculate variance */
    for (i=0; i < nstats; i++) {
      if (sptr[i].n1 < 0) continue;
      ssd = (double)sptr[i].score - mean;
      s2_score += ssd * ssd;
    }
    
    var = pu->r_u.rg.mean_var = s2_score/(double)(j-1);
    pu->r_u.rg.mean_var_sqrt = sqrt(var);
  }
  else {
    pu->r_u.rg.mu = 50.0;
    pu->r_u.rg.mean_var = 10.0;
    pu->r_u.rg.mean_var_sqrt = sqrt(10.0);
  }
  
  if (pu->r_u.rg.mean_var < 0.01) {
    /*     pu->r_u.rg.mean_var = (pu->r_u.rg.mu > 1.0) ? pu->r_u.rg.mu: 1.0; */
    return proc_hist_a(sptr, nstats, s_info, ppst, histp, do_trim, pu);
  }

  /* now remove some scores */
  if (no_trim==0) {
    nit = 5;
    pu->r_u.rg.n_trimmed = 0;
    max_trim = -BIGNUM;
    min_trim = BIGNUM;
    while (nit-- > 0) {
      t_score = 0.0;
      t_j = 0;
      for (i=0; i< nstats; i++) {
	if (sptr[i].n1 < 0) continue;
	ssd = find_zn(sptr[i].score,sptr[i].escore,sptr[i].n1,sptr[i].comp, pu);
	if ((ssd > zstrim)
#ifndef NORMAL_DIST
	    || ssd < 20.0
#else
	    || ssd < 0.0
#endif
	    )
	  {
	    /*      fprintf(stderr,"removing %3d %3d %4.1f\n",
		    sptr[i].score, sptr[i].n1,ssd); */

	    ssd = sptr[i].score;
	    if (ssd > max_trim) max_trim = ssd;
	    if (ssd < min_trim) min_trim = ssd;

	    t_score += (double)ssd;
	    t_j++;
	    pu->r_u.rg.n_trimmed++;
	    histp->entries--;
	    sptr[i].n1 = -sptr[i].n1;
	  }
      }

      if (j - t_j > 1 ) {
	mean = pu->r_u.rg.mu = (s_score - t_score)/(double)(j - t_j);

	/* calculate variance */
	s2_score = 0.0;
	for (i=0; i < nstats; i++) {
	  if ( sptr[i].n1 > 0 && sptr[i].score <= max_hscore) {
	    ssd = (double)sptr[i].score - mean;
	    s2_score += ssd * ssd;
	  }
	}
	var = pu->r_u.rg.mean_var = s2_score/(double)(j-1);
	pu->r_u.rg.mean_var_sqrt = sqrt(var);
      }
      else {
	pu->r_u.rg.mu = 50.0;
	pu->r_u.rg.mean_var = 10.0;
	pu->r_u.rg.mean_var_sqrt = sqrt(10.0);
      }

      if (pu->r_u.rg.mean_var < 0.01) {
	pu->r_u.rg.mean_var = (pu->r_u.rg.mu > 1.0) ? pu->r_u.rg.mu: 1.0;
      }

      if (pu->r_u.rg.n_trimmed < LHISTC) {
	/*
	  fprintf(stderr,"nprune %d at %d\n",nprune,nit);
	*/
	break;
      }
    }
  }

  if (ppst->zsflag < 10) s_string[0]='\0';
  else if (ppst->zs_win > 0) {
    sprintf(s_string,"(shuffled [%d], win: %d)",nstats,ppst->zs_win);
  }
  else {
    sprintf(s_string,"(shuffled [%d])",nstats);
  }

  sprintf(f_string,
#ifndef NORMAL_DIST
	  "%s Unscaled statistics: mu= %6.4f  var=%6.4f; Lambda= %6.4f",
	  s_string, pu->r_u.rg.mu,pu->r_u.rg.mean_var,PI_SQRT6/sqrt(pu->r_u.rg.mean_var_sqrt)
#else
	  "%s Unscaled normal statistics: mu= %6.4f  var=%6.4f Ztrim: %d",
	  s_string, pu->r_u.rg.mu,pu->r_u.rg.mean_var, pu->r_u.rg.n_trimmed
#endif
	  );
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

  if ( length < LENGTH_CUTOFF) {
    llen->min_score = 0;
    llen->zero_s++;
    return;
  }

  if (score < llen->min_score) llen->min_score = score;
  if (score > llen->max_score) llen->max_score = score;

  if (length > llen->max_length) llen->max_length = length;
  if (length < llen->min_length) llen->min_length = length;
  if (score > max_hscore) score = max_hscore;

#ifdef LOCAL_SCORE
  llength = (int)(LN_FACT*log((double)length)+0.5);
#else
  llength = (int)(LN_FACT*log((double)length)+0.5);
  /*  llength = length; */
#endif

  if (llength < 0 ) llength = 0;
  if (llength > llen->max) llength = llen->max;
  llen->hist[llength]++;
  dscore = (double)score;
  llen->score_sums[llength] += dscore;
  llen->score2_sums[llength] += dscore * dscore;
}

/* histogram will go from z-scores of 20 .. 100 with mean 50 and z=10 */

void
inithistz(int mh, struct hist_str *histp, double zs_off)
{
  int i, izs_off;

  izs_off = (int)(zs_off + 0.5);

  histp->z_calls = 0;

  histp->min_hist = 20 + izs_off*10;
  histp->max_hist = 120 + izs_off * 10;

  histp->histint = (int)
    ((double)(histp->max_hist - histp->min_hist + 2)/(double)mh+0.5);
  histp->maxh = (int)
    ((double)(histp->max_hist - histp->min_hist + 2)/(double)histp->histint+0.5);

  if (histp->hist_a==NULL) {
    if ((histp->hist_a=(int *)calloc((size_t)histp->maxh,sizeof(int)))==
	NULL) {
      fprintf(stderr,"*** Warning [%s:%d] cannot allocate %d for histogram\n",__FILE__, __LINE__, histp->maxh);
      histp->histflg = 0;
    }
    else histp->histflg = 1;
  }
  else {
    for (i=0; i<histp->maxh; i++) histp->hist_a[i]=0;
  }
  histp->entries = 0;
}

static double nrv[100]={
  0.3098900570,-0.0313400923, 0.1131975903,-0.2832547606, 0.0073672659,
  0.2914489107, 0.4209306311,-0.4630181404, 0.3326537896, 0.0050140359,
 -0.1117435426,-0.2835630301, 0.2302997065,-0.3102716394, 0.0819894916,
 -0.1676455701,-0.3782225018,-0.3204509938,-0.3594969187,-0.0308950398,
  0.2922813812, 0.1337170751, 0.4666577031,-0.2917784349,-0.2438179916,
  0.3002301394, 0.0231147123, 0.5687927366,-0.2318208709,-0.1476839273,
 -0.0385043851,-0.1213476523, 0.1486341995, 0.1027917167, 0.1409192644,
 -0.3280652579, 0.4232041455, 0.0775993309, 0.1159071787, 0.2769424442,
  0.3197284751, 0.1507346903, 0.0028580909, 0.4825103412,-0.0496843610,
 -0.2754357656, 0.6021881753,-0.0816123956,-0.0899148991, 0.4847183201,
  0.2151621865,-0.4542246220, 0.0690709102, 0.2461894193, 0.2126042295,
 -0.0767060668, 0.4819746149, 0.3323031326, 0.0177600676, 0.1143185210,
  0.2653977455, 0.0921872958,-0.1330986718, 0.0412287716,-0.1691604748,
 -0.0529679078,-0.0194157955,-0.6117493924, 0.1199067932, 0.0210243193,
 -0.5832259838,-0.1685528664, 0.0008591271,-0.1120347822, 0.0839125069,
 -0.2787486831,-0.1937017962,-0.1915733940,-0.7888453635,-0.3316745163,
  0.1180885226,-0.3347001067,-0.2477492636,-0.2445697600, 0.0001342482,
 -0.0015759812,-0.1516473992,-0.5202267615, 0.2136975210, 0.2500423188,
 -0.2402926401,-0.1094186280,-0.0618869933,-0.0815221188, 0.2623337275,
  0.0219427302 -0.1774469919, 0.0828245026,-0.3271952808,-0.0632898028};

void
addhistz(double zs, struct hist_str *histp)
{
  int ih, zi;
  double rv;

  if (histp == NULL) return;

  rv = nrv[histp->z_calls++ % 100];
  zi = (int)(zs + 0.5+rv );

  if ((zi >= 0) && (zi <= 120)) histp->entries++;

  if (zi < histp->min_hist) zi = histp->min_hist;
  if (zi > histp->max_hist) zi = histp->max_hist;
  
  ih = (zi - histp->min_hist)/histp->histint;

  histp->hist_a[ih]++;
}

/* addhistzp() does not increase histp->entries since addhist did it already */
/*
void
addhistzp(double zs, struct hist_str *histp)
{
  int ih, zi;
  double rv;

  rv = nrv[histp->z_calls++ %100];
  zi = (int)(zs + 0.5 + rv);

  if (zi < histp->min_hist) zi = histp->min_hist;
  if (zi > histp->max_hist) zi = histp->max_hist;
  
  ih = (zi - histp->min_hist)/histp->histint;

  histp->hist_a[ih]++;
}
*/

void
prune_hist(struct llen_str *llen, int score, int length, int max_hscore,
	   long *entries)
{
  int llength;
  double dscore;

#ifdef LOCAL_SCORE
  if (score <= 0 || length < LENGTH_CUTOFF) return;
#endif

  if (score > max_hscore) score = max_hscore;

#ifdef LOCAL_SCORE
  llength = (int)(LN_FACT*log((double)length)+0.5);
#else
  llength = (int)(LN_FACT*log((double)length)+0.5);
  /*  llength = length; */
#endif


  if (llength < 0 ) llength = 0;
  if (llength > llen->max) llength = llen->max;
  llen->hist[llength]--;
  dscore = (double)score;
  llen->score_sums[llength] -= dscore;
  llen->score2_sums[llength] -= dscore * dscore;

/*  (*entries)--; histp->entries is not yet initialized */
}  

/* find_zr uses rg.rho, rg.mu, and rg.mean_var_sqrt */
/* find_zr2 uses rho,  mu, but calculates var from regression using rho2, mu2 */

/* fit_llen: no trimming
   (1) regress scores vs log(n) using weighted variance
   (2) calculate mean variance after length regression
   (3) regress residual variance vs log(n), calculate REG2_STATS parameters
   (4) set variance cutoff for bin exclusion
*/

void
fit_llen(struct llen_str *llen, struct rstat_str *pr)
{
  int j;
  int n;
  int n_size;
  double x, y2, u, z;
  double mean_x, mean_y, var_x, var_y, covar_xy;
  double mean_y2, covar_xy2, var_y2, dllj;

  double sum_x, sum_y, sum_x2, sum_xy, sum_v, det, n_w;
  
  /* now fit scores to best linear function of log(n), using
     simple linear regression */
  
  /* find the shortest bin with data */
  for (llen->min=0; llen->min < llen->max; llen->min++)
    if (llen->hist[llen->min]) break;

  /* calculate the mean/variance for each length interval */
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
  for (j = llen->min; j < llen->max; j++) {
    if (llen->hist[j] > 1) {
      x = j + 0.5;
      dllj = (double)llen->hist[j];
      n_w += dllj/llen->score_var[j];
      sum_x +=   dllj * x / llen->score_var[j] ;
      sum_y += llen->score_sums[j] / llen->score_var[j];
      sum_x2 +=  dllj * x * x /llen->score_var[j];
      sum_xy +=  x * llen->score_sums[j]/llen->score_var[j];
    }
  }

  if (n_size < 5 ) {
    llen->fit_flag=0;
    pr->rho = 0;
    pr->mu = sum_y/n_w;
    return;
  }
  else {
    det = n_w * sum_x2 - sum_x * sum_x;
    if (det > 0.001) {
      llen->fit_flag = 1;
      pr->rho = (n_w * sum_xy  - sum_x * sum_y)/det;
      pr->rho_e = n_w/det;
      pr->mu = (sum_x2 * sum_y - sum_x * sum_xy)/det;
      pr->mu_e = sum_x2/det;
    }
    else {
      llen->fit_flag = 0;
      pr->rho = 0;
      pr->mu = sum_y/n_w;
      return;
    }
  }

  /* this code duplicates the calculation above
  det = n_w * sum_x2 - sum_x * sum_x;
  pr->rho = (n_w * sum_xy  - sum_x * sum_y)/det;
  pr->mu = (sum_x2 * sum_y - sum_x * sum_xy)/det;
  */

  n = 0;
  mean_x = mean_y = mean_y2 = 0.0;
  var_x = var_y = 0.0;
  covar_xy = covar_xy2 = 0.0;

  for (j = llen->min; j <= llen->max; j++) {
   if (llen->hist[j] > 1 ) {
      n += llen->hist[j];
      x = (double)j + 0.5;
      mean_x += (double)llen->hist[j] * x;
      mean_y += llen->score_sums[j];
      var_x += (double)llen->hist[j] * x * x;
      var_y += llen->score2_sums[j];
      covar_xy += x * llen->score_sums[j];
    }
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
  for (j = llen->min; j <= llen->max; j++) {
    if (llen->hist[j] > 1) {
      x = (double)j + 0.5;
      u = pr->rho * x + pr->mu;
      y2 = llen->score2_sums[j] - 2.0 * llen->score_sums[j] * u + llen->hist[j] * u * u;
      mean_y2 += y2;
      var_y2 += y2 * y2;
      covar_xy2 += x * y2;
    }
  }
  
  pr->mean_var = mean_y2 /= (double)n;
  pr->mean_var_sqrt = sqrt(pr->mean_var);
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

/* used to get parameters for REG2_STATS

   fit_llens: trim high variance bins
   (1) regress scores vs log(n) using weighted variance
   (2) regress residuals vs log(n)
   (3) remove high variance bins
   (4) calculate mean variance after length regression
*/

void
fit_llens(struct llen_str *llen, struct rstat_str *pr)
{
  int j;
  int n, n_u2, n_size;
  double x, y, y2, u, u2, v, z;
  double mean_x, mean_y, var_x, var_y, covar_xy;
  double mean_y2, covar_xy2;
  double mean_u2, mean_3u2, dllj;
  double sum_x, sum_y, sum_x2, sum_xy, sum_v, det, n_w;

  /* now fit scores to best linear function of log(n), using
     simple linear regression */
  
  /* find the shortest bin with data */
  for (llen->min=0; llen->min < llen->max; llen->min++)
    if (llen->hist[llen->min] > 1) break;

  /* calculate the mean/variance for each length interval */
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
  for (n_size = 0, j = llen->min; j < llen->max; j++) {
    if (llen->hist[j] > 1) {
      x = j + 0.5;
      dllj = (double)llen->hist[j];
      n_w += dllj/llen->score_var[j];
      sum_x +=   dllj * x / llen->score_var[j] ;
      sum_y += llen->score_sums[j] / llen->score_var[j];
      sum_x2 +=  dllj * x * x /llen->score_var[j];
      sum_xy +=  x * llen->score_sums[j]/llen->score_var[j];
      n_size++;
    }
  }

  pr->nb_tot = n_size;

  if (n_size < 5 ) {
    llen->fit_flag=0;
    pr->rho = 0;
    pr->mu = sum_y/n_w;
    return;
  }
  else {
    det = n_w * sum_x2 - sum_x * sum_x;
    if (det > 0.001) {
      llen->fit_flag = 1;
      pr->rho = (n_w * sum_xy  - sum_x * sum_y)/det;
      pr->rho_e = n_w/det;
      pr->mu = (sum_x2 * sum_y - sum_x * sum_xy)/det;
      pr->mu_e = sum_x2/det;
    }
    else {
      llen->fit_flag = 0;
      pr->rho = 0;
      pr->mu = sum_y/n_w;
      return;
    }
  }

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

  /* to this point, fit_llen(), fit_llens(), and fit_llen2() -- this
     function -- are the same */

  /* now regress on the residual variance */
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
  /* calculate parameters for variance regression */
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
    }
    else llen->score_var[j] = -1.0;
  }

  /* mean residual variance after both length and variance with length
     considered */
  mean_u2 = sqrt(mean_u2/(double)n_u2);
  mean_3u2 = mean_u2*3.0;

  /* trim all bins with variance > mean_3u2 */
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
  /* fit what is left */
  fit_llen(llen, pr);
}

/* s2str is used to regress against residual variance */
struct s2str {double s; int n;};
void s2_sort ( struct s2str *sptr, int n);

/* fit_llen2  - used by REGI_STATS

   (1) does the normal fit_llen() regression fitting
   (2) sorts the residual variance
   (3) excludes 5% lowest, highest bins by residual variance
   (4) recalculates mean_var without excluded bins
*/
void
fit_llen2(struct llen_str *llen, struct rstat_str *pr)
{
  int j;
  int n, n_y2, llen_delta, llen_del05;
  int n_size;
  double x, y2, u;
  double mean_x, mean_y, var_x, var_y, covar_xy;
  double mean_y2, covar_xy2, dllj;
  struct s2str *ss2;

  double sum_x, sum_y, sum_x2, sum_xy, sum_v, det, n_w;
  
  /* now fit scores to best linear function of log(n), using
     simple linear regression */
  
  /* find the shortest bin with data */
  for (llen->min=0; llen->min < llen->max; llen->min++)
    if (llen->hist[llen->min]) break;

  /* find the longest bin with data */
  for ( ; llen->max > llen->min; llen->max--)
    if (llen->hist[llen->max]) break;

  /* calculate the mean/variance for each length interval */
  for (n_size=0,j = llen->min; j < llen->max; j++) {
    if (llen->hist[j] > 1) {
      dllj = (double)llen->hist[j];
      llen->score_var[j] = llen->score2_sums[j]/dllj
	- (llen->score_sums[j]/dllj) * (llen->score_sums[j]/dllj);
      llen->score_var[j] /= (double)(llen->hist[j]-1);
      if (llen->score_var[j] <= 1.0 ) llen->score_var[j] = 1.0;
      n_size++;
    }
  }
	  
  pr->nb_tot = n_size;

  n_w = 0.0;
  sum_x = sum_y = sum_x2 = sum_xy = sum_v = 0;
  for (j = llen->min; j < llen->max; j++) {
    if (llen->hist[j] > 1) {
      x = j + 0.5;
      dllj = (double)llen->hist[j];
      n_w += (double)llen->hist[j]/dllj;
      sum_x +=   (double)llen->hist[j] * x / llen->score_var[j] ;
      sum_y += llen->score_sums[j] / llen->score_var[j];
      sum_x2 +=  dllj * x * x /llen->score_var[j];
      sum_xy +=  x * llen->score_sums[j]/llen->score_var[j];
    }
  }

  if (n_size < 5 ) {
    llen->fit_flag=0;
    pr->rho = 0;
    pr->mu = sum_y/n_w;
    return;
  }
  else {
    det = n_w * sum_x2 - sum_x * sum_x;
    if (det > 0.001) {
      llen->fit_flag = 1;
      pr->rho = (n_w * sum_xy  - sum_x * sum_y)/det;
      pr->rho_e = n_w/det;
      pr->mu = (sum_x2 * sum_y - sum_x * sum_xy)/det;
      pr->mu_e = sum_x2/det;
    }
    else {
      llen->fit_flag = 0;
      pr->rho = 0;
      pr->mu = sum_y/n_w;
      return;
    }
  }

  /*
  det = n_w * sum_x2 - sum_x * sum_x;
  pr->rho = (n_w * sum_xy  - sum_x * sum_y)/det;
  pr->mu = (sum_x2 * sum_y - sum_x * sum_xy)/det;
  */
/*   fprintf(stderr," rho1/mu1: %.2f/%.2f\n",pr->rho*LN_FACT,pr->mu); */

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

  /* to this point, fit_llen(), fit_llens(), and fit_llen2() -- this
     function -- are the same */

  /* now prepare to exclude some bins because of excess variance */

  /* allocate space for s2str */
  if ((ss2=(struct s2str *)calloc(llen->max+1,sizeof(struct s2str)))==NULL) {
    llen->fit_flag = 0;
    fprintf(stderr,"*** Warning [%s:%d] cannot allocate ss2\n",__FILE__, __LINE__);
    return;
  }

  /* calculate residual variance after fit, bin into ss2[] */
  mean_y2 = 0.0;
  n_y2 = n = 0;
  for (j = llen->min; j <= llen->max; j++) {
    if (llen->hist[j] > VHISTC) {
      n++;
      n_y2 += ss2[j].n = llen->hist[j];
      x = (double)j + 0.5;
      u = pr->rho * x + pr->mu;
      ss2[j].s = y2 = llen->score2_sums[j] - 2*llen->score_sums[j]*u + llen->hist[j]*u*u;
      mean_y2 += y2;
    }
  }

  pr->mean_var = mean_y2/(double)n_y2;
  pr->mean_var_sqrt = sqrt(pr->mean_var);

  /* sort the ss2[.squared residuals, n] by decreasing squared residuals */
  s2_sort(ss2+llen->min,llen->max-llen->min+1);
  
  /*  fprintf(stderr,"llen->min: %d, max: %d\n",llen->min,llen->max); */
  /* llen_delta has the number of bins with data */
  llen_delta = 0;
  for (j=llen->min; j<=llen->max; j++) {
    if (ss2[j].n > 1) { llen_delta++;}
  }

  llen_del05 = llen_delta/20;	/* top 5% of bins */

  /* exclude bottom 5% */
  for (j = llen->min; j<llen->min+llen_del05; j++) {
    pr->n1_trimmed += ss2[j].n;
    pr->nb_trimmed++;
  }

  /* calculate mean_y2 using middle 90% */
  mean_y2 = 0.0;
  n_y2 = 0;
  for (j = llen->min+llen_del05; j <= llen->min+llen_delta-llen_del05; j++) 
    if (ss2[j].n > 1) {
      mean_y2 += ss2[j].s;
      n_y2 += ss2[j].n;
    }

  /* exclude top 5% */
  for (j = llen->min+llen_delta-llen_del05+1; j< llen->max; j++) {
    pr->n1_trimmed += ss2[j].n;
    pr->nb_trimmed++;
  }
  
  free(ss2);

  /* return mean_var from middle 90% */
  if (n_y2 > 1) {
    pr->mean_var = mean_y2/(double)n_y2;
    pr->mean_var_sqrt = sqrt(pr->mean_var);
  }
  else {
    llen->fit_flag = 0;
  }

  /*    fprintf(stderr," rho1/mu1: %.4f/%.4f mean_var: %.4f/%d\n",
	  pr->rho*LN_FACT,pr->mu,pr->mean_var,n); */
  pr->var_e = 0.0;
}

double find_z(int score, double escore, int length, double comp, struct pstat_str *pu) {
  if (pu == NULL) {return 0.0;}
  return find_z_arr[pu->zsflag](score, escore, length, comp, pu);
}


/* REG_STATS - Z() from rho/mu/mean_var */
double find_zr(int score, double escore, int length, double comp, struct pstat_str *pu)
{
  double log_len, z;
  
  if (pu == NULL) return score;

#ifdef LOCAL_SCORE
  if (score <= 0) return 0;
#endif
  if ( length < LENGTH_CUTOFF) return 0;

#ifdef LOCAL_SCORE
  log_len = LN_FACT*log((double)(length));
#else
  /*   log_len = length; */
  log_len = LN_FACT*log((double)(length));
#endif

/*  var = pu->r_u.rg.rho2 * log_len + pu->r_u.rg.mu2;
  if (var < pu->r_u.rg.var_cutoff) var = pu->r_u.rg.var_cutoff;
*/

  z = pu->zs_off + ((double)score - pu->r_u.rg.rho * log_len - pu->r_u.rg.mu) / pu->r_u.rg.mean_var_sqrt;

  return (50.0 + z*10.0);
}

/* REG2_STATS Z() from rho/mu, rho2/mu2 */
double find_zr2(int score, double escore, int length, double comp, struct pstat_str *pu)
{
  double log_len, var;
  double z;
  
  if ( length < LENGTH_CUTOFF) return 0;

#ifdef LOCAL_SCORE
  log_len = LN_FACT*log((double)(length));
#else
  /*   log_len = length; */
  log_len = LN_FACT*log((double)(length));
#endif

  var = pu->r_u.rg.rho2 * log_len + pu->r_u.rg.mu2;
  if (var < pu->r_u.rg.var_cutoff) var = pu->r_u.rg.mean_var;

  z = pu->zs_off + ((double)score - pu->r_u.rg.rho * log_len - pu->r_u.rg.mu) / sqrt(var);

  return (50.0 + z*10.0);
}

#ifdef USE_LNSTATS
/* LN_STATS - ln()-scaled mu, mean_var */
double find_zl(int score, int length, double comp, struct pstat_str *pu)
{
  double ls, z;
  
  ls = (double)score*LN200/log((double)length);

  z = (ls - pu->r_u.rg.mu) / pu->r_u.rg.mean_var_sqrt;

  return (50.0 + z*10.0);
}
#endif

/* MLE_STATS - Z() from MLE for lambda, K */
double
find_ze(int score, double escore, int length, double comp, struct pstat_str *pu)
{
  double z, mp, np, a_n1;
  
  a_n1 = (double)length; 

  mp = pu->r_u.ag.a_n0;
  np = a_n1;

  if (np < 1.0) np = 1.0;
  if (mp < 1.0) mp = 1.0;

  z = pu->r_u.ag.Lambda * score - log(pu->r_u.ag.K * np * mp);

  z = -z + EULER_G;
  z /= - PI_SQRT6;
  z += pu->zs_off;

  return (50.0 + z*10.0);
}

/* MLE2_STATS - Z() from MLE for mle_a0..2, mle_b1, length, comp */
double
find_ze2(int score, double escore, int length, double comp, struct pstat_str *pu)
{
  double z, mp, np, a_n1;
  
  a_n1 = (double)length; 

  if (comp <= 0.0) comp = pu->r_u.m2.ave_comp;

  /* avoid very biased comp estimates */
  /* comp = exp((4.0*log(comp)+log(pu->r_u.m2.ave_comp))/5.0); */

  mp = pu->r_u.m2.a_n0;
  np = a_n1;

  if (np < 1.0) np = 1.0;
  if (mp < 1.0) mp = 1.0;

  z = (-(pu->r_u.m2.mle2_a0 + pu->r_u.m2.mle2_a1 * comp + pu->r_u.m2.mle2_a2 * comp * log(np * mp)) + score) / (pu->r_u.m2.mle2_b1 * comp);


  z = -z + EULER_G;
  z /= - PI_SQRT6;
  z += pu->zs_off;

  return (50.0 + z*10.0);
}

/* AG_STATS - Altschul-Gish Lamdba, K */
double
find_za(int score, double escore, int length, double comp, struct pstat_str *pu)
{
  double z, mp, np, a_n1, a_n1f;
  
  a_n1 = (double)length; 
  a_n1f = log(a_n1)/pu->r_u.ag.H;

  mp = pu->r_u.ag.a_n0 - pu->r_u.ag.a_n0f - a_n1f;
  np = a_n1 - pu->r_u.ag.a_n0f - a_n1f;

  if (np < 1.0) np = 1.0;
  if (mp < 1.0) mp = 1.0;

  z = pu->r_u.ag.Lambda * score - log(pu->r_u.ag.K * np * mp);

  z = -z + EULER_G;
  z /= - PI_SQRT6;

  return (50.0 + z*10.0);
}

double find_zn(int score, double escore, int length, double comp, struct pstat_str *pu)
{
  double z;
  
  z = pu->zs_off + ((double)score - pu->r_u.rg.mu) / pu->r_u.rg.mean_var_sqrt;

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

#ifndef NORMAL_DIST
  e = exp(- PI_SQRT6 * zs - EULER_G);
  return n * (e > .01 ? 1.0 - exp(-e) : e);
#else
  return n * erfc(zs/M_SQRT2)/2.0; 
#endif
}

double
zs_to_p(double zs)
{
  double e, z;

  /*  if (db.entries < 5) return 0.0; */

  z = (zs - 50.0)/10.0;

  if (z > ZS_MAX) return 0.0;

#ifndef NORMAL_DIST
  e = exp(- PI_SQRT6 * z - EULER_G);
  return (e > .01 ? 1.0 - exp(-e) : e);
#else
  return erfc(z/M_SQRT2)/2.0; 
#endif
}

double
zs_to_bit(double zs, int n0, int n1)
{
  double z, a_n0, a_n1;

  z = (zs - 50.0)/10.0;
  a_n0 = (double)n0;
  a_n1 = (double)n1;

  return (PI_SQRT6 * z + EULER_G + log(a_n0*a_n1))/M_LN2 ;
}

/* computes E-value for a given z value, assuming extreme value distribution */
double
zs_to_E(double zs,int n1, int dnaseq, long entries, struct db_str db)
{
  double e, z, k;

  /*  if (db->entries < 5) return 0.0; */

  z = (zs - 50.0)/10.0;

  if (z > ZS_MAX ) return 0.0;

  if (entries < 1) entries = db.entries;

  if (dnaseq == SEQT_DNA || dnaseq == SEQT_RNA) {
    k = (double)db.length /(double)n1;
    if (db.carry > 0) {
      k += ((double)db.carry * (double)LONG_MAX)/(double)n1;
    }
  }
  else k = (double)entries;

  if (k < 1.0) k = 1.0;

#ifndef NORMAL_DIST
  z *= PI_SQRT6;
  z += EULER_G;
  e = exp(-z);
  return k * (e > .01 ? 1.0 - exp(-e) : e);
#else
  return k * erfc(z/M_SQRT2)/2.0; 
#endif
}

#ifdef NORMAL_DIST
double np_to_z(double, int *);
double np1_to_z(double, int *);
#endif

/* compute z-score for given E()-value, assuming normal or
   extreme-value dist */
double
E_to_zs(double E, long entries)
{
  double e, z;
  int error;

  e = E/(double)entries;

#ifndef NORMAL_DIST
  z = (log(e)+EULER_G)/(- PI_SQRT6);
  return z*10.0 + 50.0;
#else
  /* this formula does not work for E() >= 1 */
  if (e >= 1.0) e = 0.99;
  z = np1_to_z(e,&error);
  if (!error) return z*10.0 + 50.0;
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

#ifndef NORMAL_DIST
  e =  exp(- PI_SQRT6 * z - EULER_G);
  return (double)entries * (e > .01 ?  exp(-e) : 1.0 - e);
#else
  return (double)entries*erf(z/M_SQRT2)/2.0; 
#endif
}

int
ELK_to_s(double e_val, int n0, int n1,
	 double Lambda, double K, double H) {

  double a_n0, a_n1;
  double a_n0f, mp, np;

  a_n0 = (double)n0;
  a_n1 = (double)n1;

  a_n0f = log(K * a_n0 * a_n1)/H;
  mp = a_n0 - a_n0f;
  np = a_n1 - a_n0f;

  if (np < 1.0) np = 1.0;
  if (mp < 1.0) mp = 1.0;

  return (int)((log(K * mp * np/e_val))/Lambda);
}

/* calculate a threshold score, given an E() value and Lambda,K,H */

int
E1_to_s(double e_val, int n0, int n1, int db_size,
	struct pstat_str *pu) {
  double mp, np, a_n0, a_n0f, a_n1;
  double zs, log_len, p_val;
  int score;

  if (pu->zsflag < 0 || n0 < LENGTH_CUTOFF || n1 < LENGTH_CUTOFF) return BIGNUM;

  a_n0 = (double)n0;
  a_n1 = (double)n1;
  zs = E_to_zs(e_val, db_size);
  /* convert from "zscore" to "zvalue" */
  zs = (zs - 50.0)/10.0;
  p_val = e_val / db_size;

  switch (pu->zsflag) {

  case AVE_STATS:
    score = zs * pu->r_u.rg.mean_var_sqrt + pu->r_u.rg.mu;
    break;

  case REG_STATS:
  case REGI_STATS:
#ifdef LOCAL_SCORE
    log_len = LN_FACT * log(a_n1);
#else
    /*     log_len = a_n1; */
    log_len = LN_FACT * log(a_n1);
#endif
    score = zs * pu->r_u.rg.mean_var_sqrt + pu->r_u.rg.rho * log_len + pu->r_u.rg.mu;
    break;

  case MLE_STATS:
    score = (int)((log( pu->r_u.ag.K * a_n0 * a_n1) - log(p_val))/pu->r_u.ag.Lambda +0.5);
    break;

  case AG_STATS:
    a_n0f = log(pu->r_u.ag.K * a_n0 * a_n1)/pu->r_u.ag.H;
    mp = a_n0 - a_n0f;
    np = a_n1 - a_n0f;

    if (np < 1.0) np = 1.0;
    if (mp < 1.0) mp = 1.0;

    score = (int)((log( pu->r_u.ag.K * mp * np) - log(p_val))/pu->r_u.ag.Lambda +0.5);
    break;

  default: 
    fprintf(stderr,"*** Warning [%s:%d] statistics method: %d not yet supported ***\n", __FILE__, __LINE__, pu->zsflag);
    score = 999;
  }

#ifndef NORMAL_DIST
  if (score < 0) score = 0;
#endif
  return score;
}

/* calculate E()-value from score, lengths, database size */
double
s_to_bit(int score, int n0, int  n1, struct pstat_str *pu) { 

  double bit;
  double a_n0, a_n1;
  double z, log_len;

  if (n0 < LENGTH_CUTOFF || n1 < LENGTH_CUTOFF) return -1.0;

  a_n0 = (double)n0;
  a_n1 = (double)n1;

  switch (pu->zsflag) {

  case AVE_STATS:
    z = ((double)score - pu->r_u.rg.mu)/pu->r_u.rg.mean_var_sqrt;
    bit = (PI_SQRT6 * z + EULER_G + log(a_n0*a_n1))/M_LN2 ;
    break;

  case REG_STATS:
  case REGI_STATS:
    log_len = LN_FACT * log(a_n1);
    z = ((double)score - pu->r_u.rg.rho * log_len - pu->r_u.rg.mu);
    z /= pu->r_u.rg.mean_var_sqrt ;
    bit = (PI_SQRT6 * z + EULER_G + log(a_n0*a_n1))/M_LN2 ;
    break;

  case MLE_STATS:
  case AG_STATS:
    bit = ((double)score * pu->r_u.ag.Lambda - log(pu->r_u.ag.K))/M_LN2;
    break;

  default: 
    fprintf(stderr,"*** Warning [%s:%d] s_to_bit -- statistics method: %d not yet supported ***\n", 
	    __FILE__, __LINE__, pu->zsflag);
    bit = -1.0;
  }

  return bit;
}

double
bit_to_E (double bit, int n0, int n1, long db_size,
	  struct pstat_str *pu)
{
  double a_n0, a_n1, a_n0f, p_val;

  a_n0 = (double)n0;
  a_n1 = (double)n1;

  if (pu->zsflag == AG_STATS) {
    a_n0f = log(pu->r_u.ag.K * a_n0 * a_n1)/pu->r_u.ag.H;

    a_n0 -= a_n0f;
    a_n1 -= a_n0f;

    if (a_n0 < 1.0) a_n0 = 1.0;
    if (a_n1 < 1.0) a_n1 = 1.0;
  }

  p_val = a_n0 * a_n1 / pow(2.0, bit);
  if (p_val > 0.01) p_val = 1.0 - exp(-p_val);

  return (double)db_size * p_val;
}


void s2_sort (struct s2str *ptr, int n)
{
  int gap, i, j;
  struct s2str tmp;

  /* shell sort using sequence from Knuth */
  /* find a gap larger than 1/3 n */
  for (gap = 1; gap < n/3; gap = 3*gap +1) ;

  /* do a shell() sort, shrinking the gap */
  for ( ; gap > 0; gap = (gap-1)/3) {
    for (i = gap; i < n; i++) {
      for (j = i - gap; j >= 0; j-= gap) {
	if (ptr[j].s >= ptr[j + gap].s) break;
	tmp.s = ptr[j].s;
	tmp.n = ptr[j].n;
	ptr[j].s = ptr[j + gap].s;
	ptr[j].n = ptr[j + gap].n;
	ptr[j + gap].s = tmp.s;
	ptr[j + gap].n = tmp.n;
      }
    }
  }
}

void last_stats() {}

void
scale_scores(struct beststr **bptr, int nbest, struct db_str db,
	     struct pstruct *ppst, struct pstat_str *rs)
{
  int i, ix;
  double zscore;

  ix = ppst->score_ix;

  if (ppst->zsflag < 0 || ppst->zsflag_f < 0) {
    /* even though no statistics, we need to sort the scores */
    for (i=0; i<nbest; i++) {
      bptr[i]->zscore = bptr[i]->rst.score[ix];
    }
    sortbestz(bptr,nbest);
  }
  else {
    for (i=0; i<nbest; i++) {
      zscore = find_z(bptr[i]->rst.score[ppst->score_ix], bptr[i]->rst.escore,
		       bptr[i]->seq->n1,bptr[i]->rst.comp,rs);
      bptr[i]->zscore = zscore;
      bptr[i]->rst.escore
	=zs_to_E(zscore,bptr[i]->seq->n1,ppst->dnaseq, ppst->zdb_size,db);
    }
    sortbeste(bptr,nbest);
  }
}

#ifdef NORMAL_DIST
/*     ALGORITHM AS241  APPL. STATIST. (1988) VOL. 37, NO. 3

       Produces the normal deviate Z corresponding to a given lower
       tail area of P; Z is accurate to about 1 part in 10**16.

       The hash sums below are the sums of the mantissas of the
       coefficients.   They are included for use in checking
       transcription.
*/

/* added 30-Dec-2017 to deal with very low E()-value thresholds */

double np1_to_z(double p1, int *fault) {
  double zs1;

  if (p1 < 0.1) {
    zs1 = np_to_z(p1, fault);
    return -zs1;
  }
  else {
    return np_to_z(1.0 - p1, fault);
  }
}

double np_to_z(double p, int *fault) {

  double q, r, ppnd16;

  double zero = 0.0, one = 1.0, half = 0.5;
  double split1 = 0.425, split2 = 5.0;
  double const1 = 0.180625, const2 = 1.6;

/*   Coefficients for P close to 0.5 */

  double a0 = 3.3871328727963666080e0;
  double a1 = 1.3314166789178437745e+2;
  double a2 = 1.9715909503065514427e+3;
  double a3 = 1.3731693765509461125e+4;
  double a4 = 4.5921953931549871457e+4;
  double a5 = 6.7265770927008700853e+4;
  double a6 = 3.3430575583588128105e+4;
  double a7 = 2.5090809287301226727e+3;
  double b1 = 4.2313330701600911252e+1;
  double b2 = 6.8718700749205790830e+2;
  double b3 = 5.3941960214247511077e+3;
  double b4 = 2.1213794301586595867e+4;
  double b5 = 3.9307895800092710610e+4;
  double b6 = 2.8729085735721942674e+4;
  double b7 = 5.2264952788528545610e+3;

  double sum_ab= 55.8831928806149014439;
/* 
  Coefficients for P not close to 0, 0.5 or 1.
*/

  double c0 = 1.42343711074968357734;
  double c1 = 4.63033784615654529590;
  double c2 = 5.76949722146069140550;
  double c3 = 3.64784832476320460504;
  double c4 = 1.27045825245236838258;
  double c5 = 2.41780725177450611770e-1;
  double c6 = 2.27238449892691845833e-2;
  double c7 = 7.74545014278341407640e-4;
  double d1 = 2.05319162663775882187;
  double d2 = 1.67638483018380384940;
  double d3 = 6.89767334985100004550e-1;
  double d4 = 1.48103976427480074590e-1;
  double d5 = 1.51986665636164571966e-2;
  double d6 = 5.47593808499534494600e-4;
  double d7 = 1.05075007164441684324e-9;

  double sum_cd=49.33206503301610289036;
/*
       Coefficients for P near 0 or 1.
*/
  double e0 = 6.65790464350110377720e0;
  double e1 = 5.46378491116411436990e0;
  double e2 = 1.78482653991729133580e0;
  double e3 = 2.96560571828504891230e-1;
  double e4 = 2.65321895265761230930e-2;
  double e5 = 1.24266094738807843860e-3;
  double e6 = 2.71155556874348757815e-5;
  double e7 = 2.01033439929228813265e-7;
  double f1 = 5.99832206555887937690e-1;
  double f2 = 1.36929880922735805310e-1;
  double f3 = 1.48753612908506148525e-2;
  double f4 = 7.86869131145613259100e-4;
  double f5 = 1.84631831751005468180e-5;
  double f6 = 1.42151175831644588870e-7;
  double f7 = 2.04426310338993978564e-15;

  double sum_ef=47.52583317549289671629;

  double sum_tmp = 0.0;

  /*
  sum_tmp = a0+a1+a2+a3+a4+a5+a6+a7+b1+b2+b3+b4+b5+b6+b7;
  if (fabs(sum_tmp - sum_ab) > 1e-12) {
    fprintf (stderr," sum_ab error: %lg %lg\n",sum_tmp,sum_ab);
    *fault = 1;
    return zero;
  }

  sum_tmp = c0+c1+c2+c3+c4+c5+c6+c7+d1+d2+d3+d4+d5+d6+d7;
  if (fabs(sum_tmp - sum_cd) > 1e-12) {
    fprintf (stderr," sum_cd error: %lg %lg\n",sum_tmp,sum_cd);
    *fault = 1;
    return zero;
  }
  sum_tmp = e0+e1+e2+e3+e4+e5+e6+e7+f1+f2+f3+f4+f5+f6+f7;
  if (fabs(sum_tmp - sum_ef) > 1e-12) {
    fprintf (stderr," sum_ef error: %lg %lg\n",sum_tmp,sum_ef);
    *fault = 1;
    return zero;
  }
  */

  *fault = 0;
  q = p - half;
  if (fabs(q) <= split1) {
    r = const1 - q * q;
    return q * (((((((a7 * r + a6) * r + a5) * r + a4) * r + a3)
		    * r + a2) * r + a1) * r + a0) /
      (((((((b7 * r + b6) * r + b5) * r + b4) * r + b3)
	 * r + b2) * r + b1) * r + one);
  }
  else {
    r = (q < zero) ?  p : one - p;
    if (r <= zero) {
      *fault = 1;
      return zero;
    }
    r = sqrt(-log(r));
    if (r <= split2) {
      r -= const2;
      ppnd16 = (((((((c7 * r + c6) * r + c5) * r + c4) * r + c3)
		  * r + c2) * r + c1) * r + c0) /
	(((((((d7 * r + d6) * r + d5) * r + d4) * r + d3)
	   * r + d2) * r + d1) * r + one);
    }
    else {
      r -= split2;
      ppnd16 = (((((((e7 * r + e6) * r + e5) * r + e4) * r + e3)
		  * r + e2) * r + e1) * r + e0) /
	(((((((f7 * r + f6) * r + f5) * r + f4) * r + f3)
	   * r + f2) * r + f1) * r + one);
    }
    if (q < zero) return -ppnd16;
    else return ppnd16;
  }
}
#endif

/* print out all pstat_str info for independent calculation */
void
pstat_info(char *info_str, int info_str_n, char *comment, struct pstat_str *pu) {
  char pstat_buf[MAX_STR];

  sprintf(pstat_buf,"%s zsflag: %d\n",comment,pu->zsflag);
  SAFE_STRNCPY(info_str,pstat_buf,info_str_n);
  sprintf(pstat_buf,"%s ngLambda: %g; ngK: %g; ngH: %g\n",comment,pu->ngLambda,pu->ngK,pu->ngH);
  SAFE_STRNCAT(info_str,pstat_buf,info_str_n);
  sprintf(pstat_buf,"%s ave_n1: %g; sample_fract: %g; zs_off: %g\n",comment,
	  pu->ave_n1,pu->sample_fract,pu->zs_off);
  if (pu->zsflag == MLE2_STATS) {	/* print r_u.m2 */
    sprintf(pstat_buf,"%s mle2_stat_str: {\n",comment);
    SAFE_STRNCAT(info_str,pstat_buf,info_str_n);
    sprintf(pstat_buf,"%s a_n0 %g;\n",comment,pu->r_u.m2.a_n0);
    SAFE_STRNCAT(info_str,pstat_buf,info_str_n);
    sprintf(pstat_buf,"%s mle2_a0: %g; mle2_a1: %g; mle2_a2: %g; mle2_b1: %g\n",comment,
	    pu->r_u.m2.mle2_a0,pu->r_u.m2.mle2_a1,pu->r_u.m2.mle2_a2,pu->r_u.m2.mle2_b1);
    SAFE_STRNCAT(info_str,pstat_buf,info_str_n);
    sprintf(pstat_buf,"%s ave_comp: %g; max_comp: %g; ave_H: %g }\n",comment,
	    pu->r_u.m2.ave_comp,pu->r_u.m2.max_comp,pu->r_u.m2.ave_H);
    SAFE_STRNCAT(info_str,pstat_buf,info_str_n);
  }
  else if (pu->zsflag  == AG_STATS || pu->zsflag == MLE_STATS) {
    sprintf(pstat_buf,"%s ag_stat_str: {\n",comment);
    SAFE_STRNCAT(info_str,pstat_buf,info_str_n);
    sprintf(pstat_buf,"%s K: %g; Lambda: %g; a_n0f: %g; a_n0: %g }\n",comment,
	    pu->r_u.ag.K,pu->r_u.ag.Lambda,pu->r_u.ag.a_n0f,pu->r_u.ag.a_n0);
    SAFE_STRNCAT(info_str,pstat_buf,info_str_n);
  }
  else {
    sprintf(pstat_buf,"%s rstat_str (LN_FACT: %.1f): {\n",comment,LN_FACT);
    SAFE_STRNCAT(info_str,pstat_buf,info_str_n);
    sprintf(pstat_buf,"%s  rho: %g; rho_e: %g; mu: %g; mu_e: %g;\n",comment,
	    pu->r_u.rg.rho,pu->r_u.rg.rho_e,pu->r_u.rg.mu,pu->r_u.rg.mu_e);
    SAFE_STRNCAT(info_str,pstat_buf,info_str_n);
    sprintf(pstat_buf,"%s  mean_var: %g; var_e: %g; mean_var_sqrt: %g\n",comment,
	    pu->r_u.rg.mean_var, pu->r_u.rg.var_e,pu->r_u.rg.mean_var_sqrt);
    SAFE_STRNCAT(info_str,pstat_buf,info_str_n);
    sprintf(pstat_buf,"%s  rho2: %g; mu2: %g; var_cutoff: %g\n",comment,
	    pu->r_u.rg.rho2, pu->r_u.rg.mu2,pu->r_u.rg.var_cutoff);
    SAFE_STRNCAT(info_str,pstat_buf,info_str_n);
    sprintf(pstat_buf,"%s  n_trimmed: %d; n1_trimmed: %d; nb_trimmed: %d; nb_tot: %d }\n",comment,
	    pu->r_u.rg.n_trimmed, pu->r_u.rg.n1_trimmed,pu->r_u.rg.nb_trimmed, pu->r_u.rg.nb_tot);
    SAFE_STRNCAT(info_str,pstat_buf,info_str_n);
  }
}
