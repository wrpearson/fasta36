
/* copyright (c) 1996, 1997, 1998, 1999 William R. Pearson and the
   U. of Virginia */

/*  $Id: showrss.c 625 2011-03-23 17:21:38Z wrp $ */
/* $Revision: 625 $  */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "defs.h"
#include "structs.h"
#include "param.h"
#include "best_stats.h"

extern double
zs_to_E(double zs, int n1, int isdna, long entries,struct db_str db);
extern double zs_to_bit(double zs, int n0, int n1);
extern double zs_to_p(double zs);

extern char *prog_func;

void showbest (FILE *fp, unsigned char **aa0, unsigned char *aa1, int maxn,
	       struct beststr **bptr, int nbest, int qlib, struct mngmsg *m_msg,
	       struct pstruct pst, struct db_str db,
	       char *info_gstring2, void **f_str)
{
  double zs;
  int score;
  char *rlabel;
  struct beststr *bbp;

  if ((rlabel=strrchr(m_msg->label,' '))==NULL) rlabel = m_msg->label;

  fprintf(fp,"\n %s - %d shuffles; ",prog_func,m_msg->shuff_max);
  if (m_msg->shuff_wid > 0)
    fprintf(fp," window shuffle, window size: %d\n",m_msg->shuff_wid);
  else
    fprintf(fp," uniform shuffle\n");

  bbp = bptr[0];

  fprintf(fp," unshuffled %s score: %d;  bits(s=%d|n_l=%d): %4.1f p(%d) < %g\n",
	  rlabel,bbp->score[0],bbp->score[0], bbp->n1,
	  zs_to_bit(bbp->zscore,m_msg->n0,bbp->n1),bbp->score[0],zs_to_p(bbp->zscore));

  fprintf(fp,"For %ld sequences, a score >= %d is expected %4.4g times\n\n", 
	  pst.zdb_size,bbp->score[0],zs_to_E(bbp->zscore,bbp->n1,0l,pst.zdb_size,db)); 
}

void showalign (FILE *fp, unsigned char *aa0, unsigned char *aa1, int maxn,
		struct beststr **bptr, int nbest,int qlib,
		const struct mngmsg *m_msg, struct pstruct *ppst,
		void *f_str, char *info_gstring2)
{
}

void
aancpy(char *to, char *from, int count,
       struct pstruct pst)
{
  char *tp;

  tp=to;
  while (count-- && *from) {
    if (*from <= pst.nsq) *tp++ = pst.sq[*(from++)];
    else *tp++ = *from++;
  }
  *tp='\0';
}
