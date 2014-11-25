/* $Id: last_tat.c 938 2012-06-04 16:15:06Z wrp $ */

/* copyright (c) 1996, 1997, 1998, 1999, 2014 by William R. Pearson and
   The Rector & Visitors of the University of Virginia */

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

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "defs.h"

#include "structs.h"
#include "param.h"
#include "mm_file.h"
#include "best_stats.h"


extern int (*ranlib) (char *str, int cnt,
	       fseek_t libpos, char *libstr,
	       struct lmf_str *lm_fd);

#define RANLIB (m_fptr->ranlib)

#define MAX_BLINE 200

int
re_getlib(unsigned char *, struct annot_str **,int, int, int, int, int, long *, long *, 
	  struct lmf_str *m_fptr);

void
do_work(unsigned char *aa0, int n0, unsigned char *aa1, int n1, int frame, 
	struct pstruct *ppst, void *f_str, int qr_flg, struct rstruct *rst);

extern void
do_opt (unsigned char *aa0, int n0, unsigned char *aa1, int n1,
	int frame, struct pstruct *pst, void *f_str,
	struct rstruct *rst);

struct lmf_str *re_openlib(struct lmf_str *, int outtty);

void sortbestz (struct beststr **bptr, int nbest);

double zs_to_E(double zs,int n1, int isdna, long entries, struct db_str db);

double scale_one_score(int ipos, double escore, struct db_str db, void *rs_str);

void sortbests (struct beststr **bptr, int nbest)
{
    int gap, i, j;
    struct beststr *tmp;

    for (gap = nbest/2; gap > 0; gap /= 2)
	for (i = gap; i < nbest; i++)
	    for (j = i - gap; j >= 0; j-= gap) {
	      if (bptr[j]->rst.score[0] >= bptr[j + gap]->rst.score[0]) break;
	      tmp = bptr[j];
	      bptr[j] = bptr[j + gap];
	      bptr[j + gap] = tmp;
	    }
}

int
last_calc(
	  unsigned char **aa0, unsigned char *aa1save, int maxn,
	  struct beststr **bptr, int nbest, 
	  const struct mngmsg *m_msg, struct pstruct *ppst
	  , void **f_str
	  , void *rstat_str)
{
  unsigned char *aa1;
  int nopt, ib;
  struct beststr *bbp;
  long loffset, l_off;
  int n0, n1;
  struct rstruct rst;
  struct lmf_str *m_fptr;
  char bline[60];
  int tat_samp, tat_inc, loop_cnt, i;
  double min_escore, ess;

  n0 = m_msg->n0;

  sortbestz(bptr,nbest);

  tat_inc = 500;
/*
  if (zs_to_E(bptr[0]->zscore,bptr[0]->n1,0,ppst->zdb_size,m_msg->db)/ 
      zs_to_E(bptr[nbest-1]->zscore,bptr[nbest-1]->n1,0,ppst->zdb_size,m_msg->db)
	< 1e-20) { tat_inc /= 4 ;}
*/

/* || (zs_to_E(bptr[0]->zscore,bptr[0]->n1,0,ppst->zdb_size,m_msg->db)< 1e-5); */

  ib = tat_samp = 0;
  for (loop_cnt = 0; loop_cnt < 5; loop_cnt++) {
    tat_samp += tat_inc;
    nopt = min(nbest,tat_samp);
    min_escore = 1000000.0;
    for ( ; ib<nopt; ib++) {
      bbp = bptr[ib];

      if (bbp->rst.score[0] < 0) break;

      if (bbp->seq->aa1b == NULL) {

	if ((m_fptr=re_openlib(bbp->mseq->m_file_p,!m_msg->quiet))==NULL) {
	  fprintf(stderr,"*** cannot re-open %s\n",bbp->mseq->m_file_p->lb_name);
	  exit(1);
	}

	RANLIB(bline,sizeof(bline),bbp->mseq->lseek,bbp->mseq->libstr,m_fptr);

	n1 = re_getlib(aa1save,NULL,maxn,m_msg->ldb_info.maxt3,
		       m_msg->ldb_info.l_overlap,bbp->mseq->cont,
		       m_msg->ldb_info.term_code,
		       &loffset,&l_off,bbp->mseq->m_file_p);
	aa1 = aa1save;
      }
      else {
	n1 = bbp->seq->n1;
	aa1 = bbp->seq->aa1b;
	loffset = bbp->seq->l_offset;
	l_off = bbp->seq->l_off;
      }

      do_opt(aa0[bbp->frame],m_msg->n0,aa1,n1,bbp->frame,ppst,
	     f_str[bbp->frame],&rst);
      memcpy(&(bbp->rst),&rst,sizeof(struct rstruct));

      if ((ess=scale_one_score(ib, bbp->rst.escore, m_msg->db, rstat_str)) < 
	min_escore) { min_escore = ess;}
      /*
      fprintf(stderr,"%d: %4d %2d %3d %.4g %.4g\n",
	      ib, bbp->rst.score[0], bbp->segnum,bbp->seglen,bbp->escore, ess);
      */
    }

    if (min_escore > m_msg->e_cut) goto done;
  }
 done:
  return ib;
}
