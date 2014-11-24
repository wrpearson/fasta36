/* $Id: ag_stats.c $ */

/* this procedure implements Altschul's pre-calculated values for lambda, K */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "alt_parms.h"

static double K, Lambda, H;
extern int
look_p(struct alt_p parm[], int gap, int ext, 
       double *K, double *Lambda, double *H);

int
ag_parm(char *pam_type, int gdelval, int ggapval)
{
  int r_v, t_gdelval, t_ggapval;

#ifdef OLD_FASTA_GAP
  t_gdelval = gdelval;
  t_ggapval = ggapval;
#else
  t_gdelval = gdelval+ggapval;
  t_ggapval = ggapval;
#endif

  if (strcmp(pam_type,"BL50")==0 || strcmp(pam_type,"BLOSUM50")==0)
      r_v = look_p(bl50_p,t_gdelval,t_ggapval,&K,&Lambda,&H);
  else if (strcmp(pam_type,"BL62")==0 || strcmp(pam_type,"BLOSUM62")==0)
      r_v = look_p(bl62_p,t_gdelval,t_ggapval,&K,&Lambda,&H);
  else if (strcmp(pam_type,"P250")==0)
      r_v = look_p(p250_p,t_gdelval,t_ggapval,&K,&Lambda,&H);
  else if (strcmp(pam_type,"P120")==0)
      r_v = look_p(p120_p,t_gdelval,t_ggapval,&K,&Lambda,&H);
  else if (strcmp(pam_type,"MD10")==0)
      r_v = look_p(md10_p,t_gdelval,t_ggapval,&K,&Lambda,&H);
  else if (strcmp(pam_type,"MD20")==0)
      r_v = look_p(md20_p,t_gdelval,t_ggapval,&K,&Lambda,&H);
  else if (strcmp(pam_type,"MD40")==0)
      r_v = look_p(md40_p,t_gdelval,t_ggapval,&K,&Lambda,&H);
  else if (strcmp(pam_type,"DNA")==0 || strcmp(pam_type,"+5/-4")==0)
      r_v = look_p(nt54_p,t_gdelval,t_ggapval,&K,&Lambda,&H);
  else r_v = 0;

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
    else if (parm[i].gap == gap && parm[i].ext > ext ) continue;
    else if (parm[i].gap == gap && parm[i].ext == ext) {
      *K = parm[i].K;
      *Lambda = parm[i].Lambda;
      *H = parm[i].H;
      return 1;
    }
    else break;
  }
  return 0;
}

int E1_to_s(double e_val, int n0, int n1) {
  double mp, np, a_n0, a_n0f, a_n1, a_n1f, u;
  int score;

  a_n0 = (double)n0;
  a_n0f = log(a_n0)/H;

  a_n1 = (double)n1;
  a_n1f = log(a_n1)/H;

  mp = a_n0 - a_n0f - a_n1f;
  np = a_n1 - a_n0f - a_n1f;

  if (np < 1.0) np = 1.0;
  if (mp < 1.0) mp = 1.0;

  /*
  e_val = K * np * mp * exp ( - Lambda * score);
  log(e_val) = log(K np mp) - Lambda * score;
  (log(K np mp)-log(e_val)) / Lambda = score;
  */
  score = (int)((log( K * mp * np) - log(e_val))/Lambda +0.5);
  if (score < 0) score = 0;
  return score;
}

double s_to_E4(int score, int n0, int  n1)
{
  double p_val;
  double mp, np, a_n0, a_n0f, a_n1, a_n1f, u;
  
  a_n0 = (double)n0;
  a_n0f = log(a_n0)/H;

  a_n1 = (double)n1;
  a_n1f = log(a_n1)/H;

  mp = a_n0 - a_n0f - a_n1f;
  np = a_n1 - a_n0f - a_n1f;

  if (np < 1.0) np = 1.0;
  if (mp < 1.0) mp = 1.0;

  p_val = K * np * mp * exp ( - Lambda * score);

  if (p_val > 0.01) p_val = 1.0 - exp(-p_val);

  return p_val * 10000.0;
}

