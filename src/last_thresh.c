/* $Id: last_thresh.c 808 2011-07-19 20:05:24Z wrp $ */

/* copyright (c) 2011, 2014 by William R. Pearson and The
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

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "defs.h"

#include "structs.h"
#include "param.h"
#include "mm_file.h"
#include "best_stats.h"

#ifdef PCOMPLIB
#include "msg.h"
void do_stage2(struct beststr **bptr, int nbest, const struct mngmsg *m_msp0,
	       int s_func, struct qmng_str *qm_msp);
#endif

static char thresh_str[MAX_STR];

int E1_to_s(double e_val, int n0, int n1, int db_size, 	void *pu);

int
last_calc(
#ifndef PCOMPLIB
	  unsigned char **aa0, unsigned char *aa1, int maxn,
#endif	  
	  struct beststr **bptr, int nbest, 
	  const struct mngmsg *m_msg, struct pstruct *ppst
#ifdef PCOMPLIB
	  , struct qmng_str *qm_msp
#else
	  , void **f_str
#endif
	  , void *pstat_str)
{

  if (ppst->zdb_size < 0 ) ppst->zdb_size = m_msg->db.entries;
  ppst->repeat_thresh = E1_to_s(ppst->e_cut, m_msg->n0, bptr[0]->seq->n1, ppst->zdb_size, pstat_str);

  ppst->other_info = thresh_str;
  sprintf(thresh_str,"Threshold: E() < %.2g score: %d\n",ppst->e_cut, ppst->repeat_thresh);

  return nbest;
}
