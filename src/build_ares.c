/* $Id: build_ares.c $ */

/* copyright (c) 2010, 2014 by William R. Pearson and The Rector &
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

/* build_ares_code is called by showbest() (in threaded/serial code) or
   p2_workcomp in PCOMPLIB code to produce the cur_ares-> chain that
   is displayed in showbest().

   For PCOMPLIB, the cur_ares->chain is passed to bbp->a_res by
   do_stage2(), where it is available to showbest();

   By using this code, the a_res chain used in either mode will be the
   same, so the code required to display an a_res should be the same.
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "defs.h"
#include "structs.h"
#include "param.h"

/* #include "mm_file.h" */
#include "best_stats.h"
#include "drop_func.h"

extern void calc_coord(int n0, int n1, long qoffset, long loffset,
		      struct a_struct *aln);

extern void calc_astruct(struct a_struct *aln_p, struct a_res_str *a_res_p, void *f_str);

/* in build_ares_code, *aa1 is separate from *seq because *seq has
   permanent information about aa1, but aa1 may be temporary

   build_ares_code() calculates various annotation strings, depending
   on what kinds of annotations are requested with -m "F# file" and -m #

   There are three fundamentally different annotation formats:
   (1) annot_var_s -- the long version shown with alignments
   (2) annot_code -- a very compact version, shown with -m 9[cC] and -m 8CC
   (3) annot_id -- a different compact version

   These need to be saved separately, and used at the right time.
   They are NOT exclusive (which older code assumed).

*/

struct a_res_str *
build_ares_code(unsigned char *aa0, int n0, 
		unsigned char *aa1, struct seq_record *seq,
		int frame, int *have_ares, int repeat_thresh,
		const struct mngmsg *m_msp, struct pstruct *ppst,
		void *f_str
		)
{
  unsigned char *aa1_ann;
  struct rstruct rst;
  struct a_res_str *my_ares_p, *cur_ares_p;
  struct a_struct *aln_p;
  struct dyn_string_str *annot_str_dyn, *align_code_dyn;
  long loffset;		/* loffset is offset from beginning of real sequence */
  long l_off;		/* l_off is the the virtual coordinate of residue 1 */
  int seqc_max, annc_max;
  char *seq_code;
  int seq_code_len, annot_str_len;
  int score_delta;
  int variant_calc_done = 0;

  align_code_dyn = init_dyn_string(2048, 2048);
  annot_str_dyn = init_dyn_string(2048, 2048);

  if (seq->annot_p) {aa1_ann = seq->annot_p->aa1_ann;}
  else aa1_ann = NULL;
  loffset = seq->l_offset;
  l_off = seq->l_off;

  if (! (*have_ares & 0x1)) {	/* we don't have an a_res, and we need one */

      my_ares_p = do_walign(aa0, n0, aa1, seq->n1,
                            frame, 
                            repeat_thresh, ppst, f_str,
                            have_ares);
  }
  else {	/* we already have the a_res */
      pre_cons(aa1,seq->n1,frame,f_str);
      my_ares_p = NULL;
  }

  /* here, we need to loop through all the alignments, and produce
     the statistics/codes for each */

  for (cur_ares_p = my_ares_p; cur_ares_p != NULL; cur_ares_p = cur_ares_p->next) {

    seqc_max = my_ares_p->nres + 4*m_msp->aln.llen+4;
    cur_ares_p->aln_code = seq_code = NULL;
    cur_ares_p->aln_code_n = seq_code_len = 0;
    cur_ares_p->annot_code = NULL;
    cur_ares_p->annot_code_n = 0;
    cur_ares_p->annot_var_s = NULL;
    cur_ares_p->annot_var_id = NULL;
    cur_ares_p->annot_var_idd = NULL;

    aln_p = &cur_ares_p->aln;

    /* this sets a number of constants, from the alignment function
       and frame, and only needs to be called once */
    aln_func_vals(frame, aln_p);

    if (m_msp->tot_show_code & (SHOW_CODE_ALIGN+SHOW_CODE_CIGAR+SHOW_CODE_EXT)) {
      cur_ares_p->aln_code = seq_code=(char *)calloc(seqc_max,sizeof(char));
      /* if we have an annotation string, allocate space for the
	 encoded annotation */

      if (seq_code != NULL) {

	calc_astruct(aln_p, cur_ares_p, f_str);

	/* we need this for offset information for calc_code, but it is
	 incomplete so we must do it again */
	
	calc_coord(m_msp->n0,seq->n1,
		 m_msp->q_offset + (m_msp->q_off-1) + (m_msp->sq0off-1),
		 loffset + (l_off-1) + (m_msp->sq1off-1),
		 aln_p);

	aln_p->lc=calc_code(aa0, m_msp->n0,
			    aa1,seq->n1, 
			    aln_p,cur_ares_p,
			    ppst,
			    align_code_dyn,
			    m_msp->ann_arr,
			    m_msp->aa0a, m_msp->annot_p,
			    aa1_ann, seq->annot_p,
			    annot_str_dyn,
			    &score_delta,
			    f_str, m_msp->pstat_void,
			    m_msp->tot_show_code);

	cur_ares_p->aln_code_n = seq_code_len = strlen(seq_code);
	if (seq_code[1] == '0' && seq_code[0] == '=') {
	  fprintf(stderr," code begins with 0: %s\n", seq_code);
	}

	if (align_code_dyn != NULL) {
	  if (seq_code && cur_ares_p->aln_code == seq_code) {
	    free(seq_code);	/* free it since it is replaced below */
	  }
	  seq_code_len = strlen(align_code_dyn->string);
	  cur_ares_p->aln_code = (char *)calloc(seq_code_len+2,sizeof(char));
	  cur_ares_p->aln_code_n = seq_code_len+2;
	  SAFE_STRNCPY(cur_ares_p->aln_code,align_code_dyn->string, seq_code_len+2);
	  reset_dyn_string(align_code_dyn);
	}

	if (annot_str_dyn != NULL) {
	  annot_str_len = strlen(annot_str_dyn->string);
	  cur_ares_p->annot_code = (char *)calloc(annot_str_len+2,sizeof(char));
	  SAFE_STRNCPY(cur_ares_p->annot_code,annot_str_dyn->string, annot_str_len+2);
	  reset_dyn_string(annot_str_dyn);
	}
	else {annot_str_len = 0;}
	cur_ares_p->annot_code_n = annot_str_len;
	variant_calc_done = 1;
      }
    }

    if ((m_msp->tot_show_code & SHOW_CODE_IDD) == SHOW_CODE_IDD) {
      aln_p->lc=calc_idd(aa0,m_msp->n0,aa1,seq->n1,
			 aln_p, cur_ares_p,
			 ppst,
			 m_msp->annot_p, seq->annot_p,
			 &score_delta,
			 annot_str_dyn, f_str);
      variant_calc_done = 1;

      if (annot_str_dyn != NULL) {
	  annot_str_len = strlen(annot_str_dyn->string);
	  cur_ares_p->annot_var_idd = (char *)calloc(annot_str_len+2,sizeof(char));
	  SAFE_STRNCPY(cur_ares_p->annot_var_idd,annot_str_dyn->string, annot_str_len+2);
	}
	else {annot_str_len = 0;}

    }

    /* ensure that calc_id (or something else) is ALWAYS done to set score_delta */
    if (!variant_calc_done || (m_msp->tot_show_code & SHOW_CODE_ID) == SHOW_CODE_ID) {
      aln_p->lc=calc_id(aa0,m_msp->n0,aa1,seq->n1,
			aln_p, cur_ares_p,
			ppst,
			m_msp->annot_p, seq->annot_p,
			&score_delta,
			annot_str_dyn, f_str);
      
      if ((m_msp->tot_show_code & SHOW_CODE_ID)==SHOW_CODE_ID && (annot_str_dyn->string[0] != '\0')) {
	if ((cur_ares_p->annot_var_id = (char *)calloc(strlen(annot_str_dyn->string)+2, sizeof(char)))==NULL) {
	  fprintf(stderr,"*** ERROR *** [%s/%d] cannot allocate cur_ares_p->annot_var_s [%d]\n",
		  __FILE__, __LINE__, (int)strlen(annot_str_dyn->string)+2);
	}
	else {
	  strncpy(cur_ares_p->annot_var_id,annot_str_dyn->string,strlen(annot_str_dyn->string)+2);
	}
      }
    }

    if (score_delta > 0) {
      cur_ares_p->rst.score[0] += score_delta;
      cur_ares_p->rst.score[1] += score_delta;
      cur_ares_p->rst.score[2] += score_delta;
      cur_ares_p->sw_score += score_delta;
      cur_ares_p->score_delta = score_delta;
    }
    else {
      cur_ares_p->score_delta = 0;
    }

    /*
    if (annot_str_dyn->string[0] != '\0') {
      if ((cur_ares_p->annot_var_s = (char *)calloc(strlen(annot_str_dyn->string)+2, sizeof(char)))==NULL) {
	fprintf(stderr,"*** ERROR *** [%s/%d] cannot allocate cur_ares_p->annot_var_s [%d]\n",
		__FILE__, __LINE__, (int)strlen(annot_str_dyn->string)+2);
      }
      else {
	strncpy(cur_ares_p->annot_var_s,annot_str_dyn->string,strlen(annot_str_dyn->string)+2);
      }
    }
    */
	
    /* this should be all the information we need on the alignment */
  } /* end for (cur_ares_p;) */
  free_dyn_string(annot_str_dyn);
  free_dyn_string(align_code_dyn);
  return my_ares_p;
}
