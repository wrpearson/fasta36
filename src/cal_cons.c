/* cal_cons.c - routines for printing translated alignments for
   fasta, ssearch, ggsearch, glsearch  */

/* $Id: cal_cons.c 1280 2014-08-21 00:47:55Z wrp $ */
/* $Revision: 1280 $  */

/* copyright (c) 1998, 1999, 2007, 2014 by William R. Pearson and The
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


/* removed from dropgsw2.c, dropnfa.c April, 2007 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>

#include "defs.h"
#include "param.h"
#include "dyn_string.h"

#if defined(FASTA) || defined(TFASTA)
#include "dropnfa.h"
#endif

#if defined(SSEARCH) || defined(OSEARCH)
#include "dropgsw2.h"
#endif

#ifdef LALIGN
#include "dropgsw2.h"
#endif

#include "a_mark.h"

struct update_code_str {
  int p_op_idx;
  int p_op_cnt;
  int show_code;
  int cigar_order;
  int show_ext;
  char *op_map;
};

static char *ori_code = "=-+*x";
static char *cigar_code = "MDIMX";

static struct update_code_str *
init_update_data(int show_code);

static void
sprintf_code(char *tmp_str, struct update_code_str *, int op_idx, int op_cnt);

static void 
update_code(char *al_str, int al_str_max, 
	    struct update_code_str *update_data, int op, 
	    int sim_code, unsigned char sp0, unsigned char sp1);

static void
close_update_data(char *al_str, int al_str_max,
		  struct update_code_str *update_data);

extern void aancpy(char *to, char *from, int count, struct pstruct *ppst);
extern void *init_stack(int, int);
extern void push_stack(void *, void *);
extern void *pop_stack(void *);
extern void *free_stack(void *);

/* returns M_NEG, M_ZERO, M_POS, M_IDENT, M_DEL (a_mark.h)
   updates *aln->nsim, npos, nident, nmismatch */
extern int
align_type(int score, char sp0, char sp1, int nt_align, struct a_struct *aln, int pam_x_id_sim);

extern void
process_annot_match(int *itmp, int *pam2aa0v, 
		    long ip, long ia, char *sp1, char *sp1a, const unsigned char *sq,
		    struct annot_entry *annot_arr_p, char **ann_comment,
		    void *annot_stack, int *have_push_features, int *v_delta,
		    int *d_score_p, int *d_ident_p, int *d_alen_p, struct domfeat_link **left_domain_p,
		    long *left_end_p, int init_score);

extern int
next_annot_match(int *itmp, int *pam2aa0v, 
		 long ip, long ia, char *sp1, char *sp1a, const unsigned char *sq,
		 int i_annot, int n_annot, struct annot_entry **annot_arr, char **ann_comment,
		 void *annot_stack, int *have_push_features, int *v_delta,
		 int *d_score_p, int *d_ident_p, int *d_alen_p, struct domfeat_link **left_domain,
		 long *left_domain_end, int init_score);

extern void
close_annot_match (int ia, void *annot_stack, int *have_push_features,
		   int *d_score_p, int *d_ident_p, int *d_alen_p,
		   struct domfeat_link **left_domain_p,
		   long *left_end_p, int init_score);

extern void
comment_var(long i0, char sp0, long i1, char sp1, char o_sp1, char sim_char,
	    const char *ann_comment, struct dyn_string_str *annot_var_dyn,
	    int target, int d_type);

void
display_push_features(void *annot_stack, struct dyn_string_str *annot_var_dyn,
		      long i0_pos, char sp0, long i1_pos, char sp1, char sym, 
		      int score, double comp, int n0, int n1,
		      void *pstat_void, int d_type);

#define DP_FULL_FMT 1	/* Region: score: bits: id: ... */

extern int seq_pos(int pos, int rev, int off);

int calc_cons_a(const unsigned char *aa0, int n0,
		const unsigned char *aa1, int n1,
		int *nc,
		struct a_struct *aln,
		struct a_res_str *a_res,
		struct pstruct *ppst,
		char *seqc0, char *seqc1, char *seqca, int *cumm_seq_score,
		const unsigned char *ann_arr,
		const unsigned char *aa0a, const struct annot_str *annot0_p, char *seqc0a,
		const unsigned char *aa1a, const struct annot_str *annot1_p, char *seqc1a,
		int *score_delta, 
		struct dyn_string_str *annot_var_dyn,
		struct f_struct *f_str,
		void *pstat_void)
{
  int i0, i1, nn1;
  int op, lenc, nd, ns, itmp, p_match;
  const unsigned char *aa1p;
  char *sp0, *sp0a, *sp1, *sp1a, *spa, t_spa;
  int *i_spa;
  const unsigned char *sq;
  int *rp;
  int smins, mins, ntmp;
  int have_ann = 0;
  void *annot_stack;

  /* variables for variant changes */
  int *aa0_pam2_p;
  char *sim_sym = aln_map_sym[5];
  struct annot_entry **s_annot0_arr_p, **s_annot1_arr_p;

  char *ann_comment;
  int i0_annot, i1_annot;	/* i0_annot, i1_annot, count through
				   the list of annotations */
  long i0_left_end, i1_left_end; /* left-most coordinate of domain end */

  int v_delta, v_tmp;
  int d1_score, d1_ident, d1_alen;
  int d0_score, d0_ident, d0_alen;
  int have_push_features;
  struct domfeat_link *left_domain_list1, *left_domain_list0;

  /* variables for handling coordinate offsets */
  long q_offset, l_offset;

  *score_delta = 0;
  i0_left_end = i1_left_end = -1;
  left_domain_list0 = left_domain_list1 = NULL;
  d1_score = d1_ident = d1_alen = 0;
  d0_score = d0_ident = d0_alen = 0;

  NULL_dyn_string(annot_var_dyn);
  have_ann = (seqc0a != NULL); 

  if (ppst->ext_sq_set) {
    sq = ppst->sqx;
  }
  else {
    sq = ppst->sq;
  }

#ifndef TFASTA
  aa1p = aa1;
  nn1 = n1;
#else
  aa1p = f_str->aa1x;
  nn1 = f_str->n10;
#endif

  aln->amin0 = a_res->min0;
  aln->amax0 = a_res->max0;
  aln->amin1 = a_res->min1;
  aln->amax1 = a_res->max1;
  aln->calc_last_set = 1;

  q_offset = aln->q_offset;
  l_offset = aln->l_offset;

#ifndef LCAL_CONS
  /* will we show all the start ?*/
  if (min(a_res->min0,a_res->min1)<aln->llen || aln->showall==1)
    if (a_res->min0 >= a_res->min1) {              /* aa0 extends more to left */
      smins=0;
      if (aln->showall==1) mins = a_res->min0;
      else mins = min(a_res->min0,aln->llcntx);
      aancpy(seqc0,(char *)aa0+a_res->min0-mins,mins,ppst);
      aln->smin0 = a_res->min0-mins;
      if ((mins-a_res->min1)>0) {
	memset(seqc1,' ',mins-a_res->min1);
	aancpy(seqc1+mins-a_res->min1,(char *)aa1p,a_res->min1,ppst);
	aln->smin1 = 0;
      }
      else {
	aancpy(seqc1,(char *)aa1p+a_res->min1-mins,mins,ppst);
	aln->smin1 = a_res->min1-mins;
      }
    }
    else {
      smins=0;
      if (aln->showall == 1) mins=a_res->min1;
      else mins = min(a_res->min1,aln->llcntx);
      aancpy(seqc1,(char *)(aa1p+a_res->min1-mins),mins,ppst);
      aln->smin1 = a_res->min1-mins;
      if ((mins-a_res->min0)>0) {
	memset(seqc0,' ',mins-a_res->min0);
	aancpy(seqc0+mins-a_res->min0,(char *)aa0,a_res->min0,ppst);
	aln->smin0 = 0;
      }
      else {
	aancpy(seqc0,(char *)aa0+a_res->min0-mins,mins,ppst);
	aln->smin0 = a_res->min0-mins;
      }
    }
  else {
    mins= min(aln->llcntx,min(a_res->min0,a_res->min1));
    smins=mins;
    aln->smin0=a_res->min0 - smins;
    aln->smin1=a_res->min1 - smins;
    aancpy(seqc0,(char *)aa0+a_res->min0-mins,mins,ppst);
    aancpy(seqc1,(char *)aa1p+a_res->min1-mins,mins,ppst);
  }
  /* set the alignment code to zero for context */
  memset(seqca,0,mins);
  if (have_ann) {
    memset(seqc0a,' ',mins);
    memset(seqc1a,' ',mins);
  }
#else		/* no flanking context */
  smins = mins = 0;
  aln->smin0=a_res->min0;
  aln->smin1=a_res->min1;
#endif

  /* now get the middle */

  spa = seqca+mins;
  if (cumm_seq_score) i_spa = cumm_seq_score+mins;
  sp0 = seqc0+mins;
  sp1 = seqc1+mins;

  if (have_ann) {
    sp0a = seqc0a+mins;
    sp1a = seqc1a+mins;
  }

  rp = a_res->res;
  lenc = aln->nident = aln->nmismatch =
    aln->npos = aln->nsim = aln->ngap_q = aln->ngap_l = aln->nfs = op = 0;
  p_match = 1;
  i0 = a_res->min0;
  i1 = a_res->min1;

  v_delta = 0;
  i0_annot = i1_annot = 0;
  annot_stack = NULL;
  s_annot0_arr_p = s_annot1_arr_p = NULL;
  have_push_features=0;
  if (have_ann) {
    if ((annot1_p && annot1_p->n_annot>0) || (annot0_p && annot0_p->n_annot > 0)) {annot_stack = init_stack(64,64);}
    if (annot1_p && annot1_p->n_annot > 0) {
      s_annot1_arr_p = annot1_p->s_annot_arr_p;

      while (i1_annot < annot1_p->n_annot) {
	if (s_annot1_arr_p[i1_annot]->pos >= i1 + l_offset) {break;}
	if (s_annot1_arr_p[i1_annot]->end < i1 + l_offset) {i1_annot++; continue;}

	if (s_annot1_arr_p[i1_annot]->label == '-') {
	  process_annot_match(&itmp, aa0_pam2_p, l_offset+seq_pos(i1,aln->llrev,0), q_offset + seq_pos(i0,aln->qlrev,0),
			      sp1, sp1a, sq, s_annot1_arr_p[i1_annot],  &ann_comment, 
			      annot_stack, &have_push_features, &v_delta,
			      &d1_score, &d1_ident, &d1_alen, &left_domain_list1, &i1_left_end, 0);
	}
	i1_annot++;
      }
    }

    if (annot0_p && annot0_p->n_annot>0) {
      s_annot0_arr_p = annot0_p->s_annot_arr_p;

      while (i0_annot < annot0_p->n_annot && s_annot0_arr_p[i0_annot]->pos < i0 + q_offset) {
	if (s_annot0_arr_p[i0_annot]->pos >= i0 + q_offset) {break;}
	if (s_annot0_arr_p[i0_annot]->end < i0 + q_offset) {i0_annot++; continue;}

	if (s_annot0_arr_p[i0_annot]->label == '-') {
	  process_annot_match(&itmp, aa0_pam2_p, q_offset+seq_pos(i0,aln->qlrev,0), l_offset + seq_pos(i1,aln->llrev,0),
			      sp0, sp0a, sq, s_annot0_arr_p[i0_annot], &ann_comment, 
			      annot_stack, &have_push_features, &v_delta,
			      &d0_score, &d0_ident, &d0_alen, &left_domain_list0, &i0_left_end, 0);
	}
	i0_annot++;
      }
    }
  }

  while (i0 < a_res->max0 || i1 < a_res->max1) {
    /* match/mismatch (aligned residues */
    /* here, op is the "current" encoding, and *rp is the next one */
    if (op == 0 && *rp == 0) {
      op = *rp++;
      lenc++;

      if (ppst->pam_pssm) {aa0_pam2_p = ppst->pam2p[0][i0];}
      else {aa0_pam2_p = ppst->pam2[0][aa0[i0]];}

      itmp=aa0_pam2_p[aa1p[i1]];

      *sp0 = sq[aa0[i0]];
      *sp1 = sq[aa1p[i1]];

      if (have_ann) {
	have_push_features = 0;
	*sp0a = *sp1a = ' ';
	if (aa0a) {*sp0a = ann_arr[aa0a[i0]];}
	if (aa1a) {*sp1a = ann_arr[aa1a[i1]];}
	if (s_annot1_arr_p) {
	  if (i1+l_offset == s_annot1_arr_p[i1_annot]->pos || i1+l_offset == i1_left_end) {

	    i1_annot = next_annot_match(&itmp, aa0_pam2_p, l_offset+seq_pos(i1,aln->llrev,0),
					q_offset + seq_pos(i0,aln->qlrev,0), sp1, sp1a, sq, 
					i1_annot, annot1_p->n_annot, s_annot1_arr_p,
					&ann_comment, annot_stack, &have_push_features, &v_delta,
					&d1_score, &d1_ident, &d1_alen, &left_domain_list1, &i1_left_end, 0);

	    /* must be out of the loop to capture the last value */
	    if (sq[aa1p[i1]] != *sp1) {
	      t_spa = align_type(itmp, *sp0, *sp1, ppst->nt_align, NULL, ppst->pam_x_id_sim);

	      comment_var(q_offset+seq_pos(i0,aln->qlrev,0), *sp0,
			  l_offset+seq_pos(i1,aln->llrev,0), *sp1,
			  sq[aa1p[i1]], sim_sym[t_spa], ann_comment,
			  annot_var_dyn, 1, 1);
	    }
	  }
	  d1_score += itmp;
	}

	if (s_annot0_arr_p) {
	  if (i0 + q_offset == s_annot0_arr_p[i0_annot]->pos || i0 + q_offset == i0_left_end) {

	    i0_annot = next_annot_match(&itmp, aa0_pam2_p, q_offset+seq_pos(i0,aln->qlrev,0),
					l_offset + seq_pos(i1,aln->llrev,0), sp0, sp0a, sq, 
					i0_annot, annot0_p->n_annot, s_annot0_arr_p,
					&ann_comment, annot_stack, &have_push_features, &v_delta,
					&d0_score, &d0_ident, &d0_alen, &left_domain_list0, &i0_left_end, 0);

	    /* must be out of the loop to capture the last value */
	    if (sq[aa0[i0]] != *sp0) {
	      t_spa = align_type(itmp, *sp0, *sp1, ppst->nt_align, NULL, ppst->pam_x_id_sim);

	      comment_var(q_offset+seq_pos(i0,aln->qlrev,0), *sp0,
			  l_offset+seq_pos(i1,aln->llrev,0), *sp1,
			  sq[aa0[i0]], sim_sym[t_spa], ann_comment,
			  annot_var_dyn, 0, 1);
	    }
	  }
	  d0_score += itmp;
	}
	sp0a++; sp1a++;
      }

      *spa = align_type(itmp, *sp0, *sp1, ppst->nt_align, aln, ppst->pam_x_id_sim);

      d1_alen++;
      d0_alen++;
      if (*spa == M_IDENT) {
	d1_ident++;
	d0_ident++;
      }

      /* now we have done all the ?modified identity checks, display
	 potential site annotations */
      if (have_ann && have_push_features) {
	display_push_features(annot_stack, annot_var_dyn,
			      q_offset+seq_pos(i0,aln->qlrev,0), *sp0,
			      l_offset+seq_pos(i1,aln->llrev,0), *sp1,
			      sim_sym[*spa],
			      a_res->rst.score[ppst->score_ix], a_res->rst.comp, n0, n1,
			      pstat_void, DP_FULL_FMT);
	have_push_features = 0;
      }

      if (cumm_seq_score) *i_spa++ = itmp;
      i0++; i1++;
      sp0++; sp1++; spa++;
    }
    else {	/* indel */
      if (op==0) {
	op = *rp++;
	if (cumm_seq_score) *i_spa = ppst->gdelval;
	d1_score +=  ppst->gdelval;
	d0_score +=  ppst->gdelval;
      }
      if (cumm_seq_score) *i_spa++ += ppst->ggapval;
      d1_score +=  ppst->ggapval; d1_alen++;
      d0_score +=  ppst->ggapval; d0_alen++;

      if (op>0) {	/* insertion in aa0 */
	*sp1 = sq[aa1p[i1]];
	*sp0 = '-';
	*spa++ = M_DEL;
	if (have_ann) {
	  have_push_features = 0;
	  *sp0a++ = ' ';
	  if (aa1a) *sp1a = ann_arr[aa1a[i1]];
	  else *sp1a = ' ';
	  if (s_annot1_arr_p) {
	    if (i1+l_offset == s_annot1_arr_p[i1_annot]->pos || i1+l_offset == i1_left_end) {
	      i1_annot = next_annot_match(&itmp, aa0_pam2_p, l_offset+seq_pos(i1,aln->llrev,0),
					  q_offset+seq_pos(i0,aln->qlrev,0), sp1, sp1a, sq, 
					  i1_annot, annot1_p->n_annot, s_annot1_arr_p,
					  &ann_comment, annot_stack, &have_push_features, &v_delta,
					  &d1_score, &d1_ident, &d1_alen, &left_domain_list1, &i1_left_end,
					  ppst->ggapval+ppst->gdelval);
	    }
	  }

	  if (have_push_features) {
	    display_push_features(annot_stack, annot_var_dyn,
				  q_offset+seq_pos(i0,aln->qlrev,0), *sp0,
				  l_offset+seq_pos(i1,aln->llrev,0), *sp1,
				  sim_sym[*spa],
				  a_res->rst.score[ppst->score_ix], a_res->rst.comp, n0, n1,
				  pstat_void, DP_FULL_FMT);
	    have_push_features = 0;
	  }
	  sp1a++;
	}

	sp0++;
	sp1++;
	i1++;
	op--;
	lenc++;
	aln->ngap_q++;
      }
      else {		/* insertion in aa1 */
	*sp1 = '-';
	*spa++ = M_DEL;
	*sp0 = sq[aa0[i0]];
	if (have_ann) {
	  have_push_features = 0;
	  *sp1a++ = ' ';
	  if (aa0a) *sp0a = ann_arr[aa0a[i0]];
	  else *sp0a = ' ';
	  if (s_annot0_arr_p) {
	    if (i0+q_offset == s_annot0_arr_p[i0_annot]->pos || i0+q_offset == i0_left_end) {
	      i0_annot = next_annot_match(&itmp, ppst->pam2[0][aa1[i1]], q_offset+seq_pos(i0,aln->qlrev,0),
					  l_offset+seq_pos(i1,aln->llrev,0), sp0, sp0a, sq, 
					  i0_annot, annot0_p->n_annot, s_annot0_arr_p,
					  &ann_comment, annot_stack, &have_push_features, &v_delta,
					  &d0_score, &d0_ident, &d0_alen, &left_domain_list0, &i0_left_end,
					  ppst->ggapval+ppst->gdelval);

	    }
	  }

	  if (have_push_features) {
	    display_push_features(annot_stack, annot_var_dyn,
				  q_offset+seq_pos(i0,aln->qlrev,0), *sp0,
				  l_offset+seq_pos(i1,aln->llrev,0), *sp1,
				  sim_sym[*spa],
				  a_res->rst.score[ppst->score_ix], a_res->rst.comp, n0, n1, pstat_void, DP_FULL_FMT);
	    have_push_features = 0;
	  }
	  sp0a++;
	}

	i0++;
	sp0++;
	sp1++;
	op++;
	lenc++;
	aln->ngap_l++;
      }
    }
  }

  *score_delta = v_delta;

  *nc = lenc;
  if (have_ann) {
    *sp0a = *sp1a = '\0';
    have_push_features = 0;
    /* check for left ends after alignment */
    if (annot1_p && i1_left_end > 0) {
      close_annot_match(-1, annot_stack, &have_push_features,
			&d1_score, &d1_ident, &d1_alen, &left_domain_list1, &i1_left_end,
			0);
    }

    if (annot0_p && i0_left_end > 0) {
      close_annot_match(-1, annot_stack, &have_push_features,
			&d0_score, &d0_ident, &d0_alen, &left_domain_list0, &i0_left_end,
			0);
    }

    if (have_push_features) {
      display_push_features(annot_stack, annot_var_dyn,
			    a_res->max0-1 + q_offset, *sp0, 
			    a_res->max1-1 + l_offset, *sp1,
			    sim_sym[*spa],
			    a_res->rst.score[ppst->score_ix], a_res->rst.comp, n0, n1, pstat_void, DP_FULL_FMT);
    }
  }

  *spa = '\0';

#ifndef LCAL_CONS	/* have context around alignment */
  /*	now we have the middle, get the right end */
  if (!aln->llcntx_set) {
    ns = mins + lenc + aln->llen;	/* show an extra line? */
    ns -= (itmp = ns %aln->llen);	/* itmp = left over on last line */
    if (itmp>aln->llen/2) ns += aln->llen;  /* more than 1/2 , use another*/
    nd = ns - (mins+lenc);		/* this much extra */
  }
  else nd = aln->llcntx;

  if (nd > max(n0-a_res->max0,nn1-a_res->max1)) 
    nd = max(n0-a_res->max0,nn1-a_res->max1);
  
  if (aln->showall==1) {
    nd = max(n0-a_res->max0,nn1-a_res->max1);	/* reset for showall=1 */
    /* get right end */
    aancpy(seqc0+mins+lenc,(char *)aa0+a_res->max0,n0-a_res->max0,ppst);
    aancpy(seqc1+mins+lenc,(char *)aa1p+a_res->max1,nn1-a_res->max1,ppst);
    /* fill with blanks - this is required to use one 'nc' */
    memset(seqc0+mins+lenc+n0-a_res->max0,' ',nd-(n0-a_res->max0));
    memset(seqc1+mins+lenc+nn1-a_res->max1,' ',nd-(nn1-a_res->max1));
  }
  else {
    if ((nd-(n0-a_res->max0))>0) {
      aancpy(seqc0+mins+lenc,(char *)aa0+a_res->max0,(n0-a_res->max0),ppst);
      memset(seqc0+mins+lenc+n0-a_res->max0,' ',nd-(n0-a_res->max0));
    }
    else {
      aancpy(seqc0+mins+lenc,(char *)aa0+a_res->max0,nd,ppst);
    }

    if ((nd-(nn1-a_res->max1))>0) {
      aancpy(seqc1+mins+lenc,(char *)aa1p+a_res->max1,nn1-a_res->max1,ppst);
      memset(seqc1+mins+lenc+nn1-a_res->max1,' ',nd-(nn1-a_res->max1));
    }
    else {
      aancpy(seqc1+mins+lenc,(char *)aa1p+a_res->max1,nd,ppst);
    }
  }
  if (have_ann) {
    memset(seqc0a+mins+lenc,' ',nd);
    memset(seqc1a+mins+lenc,' ',nd);
    /*
      ntmp = nd-(n0-a_res->max0);
      if (ntmp > 0) memset(seqc0a+mins+lenc+n0-a_res->max0,' ',ntmp);
      ntmp = nd-(nn1-a_res->max1);
      if (ntmp > 0) memset(seqc1a+mins+lenc+nn1-a_res->max1,' ',ntmp);
    */
  }
#else
  nd = 0;
#endif

  /*  fprintf(stderr,"%d\n",mins+lenc+nd); */

  lenc = mins + lenc + nd;

  if (annot0_p || annot1_p) free_stack(annot_stack);

  return lenc;
}

void
calc_astruct(struct a_struct *aln_p, struct a_res_str *a_res_p, struct f_struct *f_str) {

  /* we do not pay attention to aln_p->calc_last_set, because all the
     functions (calc_astruct, calc_cons_a, calc_code) use exactly the same
     assignment */

  aln_p->amin0 = a_res_p->min0;
  aln_p->amax0 = a_res_p->max0;
  aln_p->amin1 = a_res_p->min1;
  aln_p->amax1 = a_res_p->max1;
}



static struct update_code_str *
init_update_data(show_code) {

  struct update_code_str *update_data_p;

  if ((update_data_p = (struct update_code_str *)calloc(1,sizeof(struct update_code_str)))==NULL) {
    fprintf(stderr,"*** error [%s:%d] - init_update_data(): cannot allocate update_code_str\n",
	      __FILE__, __LINE__);
    return NULL;
  }

  update_data_p->p_op_idx = -1;
  update_data_p->p_op_cnt = 0;
  update_data_p->show_code = show_code;

  if ((show_code & SHOW_CODE_MASK) == SHOW_CODE_CIGAR) {
    update_data_p->op_map = cigar_code;
    update_data_p->cigar_order = 1;
  }
  else {
    update_data_p->op_map = ori_code;
    update_data_p->cigar_order = 0;
  }

  if ((show_code & SHOW_CODE_EXT) == SHOW_CODE_EXT) {
    update_data_p->show_ext = 1;
  }
  else {
    update_data_p->show_ext = 0;
  }

  return update_data_p;
}

static void
close_update_data(char *al_str, int al_str_max, 
		  struct update_code_str *up_dp) {
  char tmp_cnt[MAX_SSTR];

  if (!up_dp) return;
  sprintf_code(tmp_cnt,up_dp, up_dp->p_op_idx, up_dp->p_op_cnt);
  strncat(al_str,tmp_cnt,al_str_max);

  free(up_dp);
}

/* update_code() has been modified to work more correctly with
   ggsearch/glsearch, which, because alignments can start with either
   insertions or deletions, can produce an initial code of "0=".  When
   that happens, it is ignored and no code is added.

   *al_str - alignment string [al_str_max] - not dynamic
   op -- encoded operation, currently 0=match, 1-delete, 2-insert, 3-term-match, 4-mismatch
   op_cnt -- length of run
   show_code -- SHOW_CODE_CIGAR uses cigar_code, otherwise legacy
*/

static void
sprintf_code(char *tmp_str, struct update_code_str *up_dp, int op_idx, int op_cnt) {
  if (up_dp->cigar_order) {
    sprintf(tmp_str,"%d%c",op_cnt,up_dp->op_map[op_idx]);
  }
  else {
    sprintf(tmp_str,"%c%d",up_dp->op_map[op_idx],op_cnt);
  }
}

static void
update_code(char *al_str, int al_str_max, 
	    struct update_code_str *up_dp, int op, 
	    int sim_code,  unsigned char sp0, unsigned char sp1)
{
  char tmp_cnt[MAX_SSTR];

  /* op == 0 : match state (could involve termination codons);
     op == 1 : deletion
     op == 2 : insertion
     op == 3 : *:*
     p_op == 5 : mismatch state
  */

  if (up_dp->p_op_cnt == 0) {
    up_dp->p_op_idx = op;
    up_dp->p_op_cnt = 1;
    return;
  }

  if (op == 1 || op == 2) {
    if (up_dp->p_op_idx == op) { up_dp->p_op_cnt++;}
    else {
      sprintf_code(tmp_cnt,up_dp, up_dp->p_op_idx,up_dp->p_op_cnt);
      strncat(al_str,tmp_cnt,al_str_max);
      up_dp->p_op_idx = op;
      up_dp->p_op_cnt = 1;
    }
  }
  else if (op==0 || op == 3) {
    if (sp0 != '*' && sp1 != '*') {	/* default case, not termination */
      if (up_dp->show_ext) {
	if (sim_code != M_IDENT) { op = 4;}
      }
    }
    else {	/* have a termination codon, output for !SHOW_CODE_CIGAR */
      if (!up_dp->cigar_order) {
	if (sp0 == '*' || sp1 == '*') { op = 3;}
      }
      else if (up_dp->show_ext && (sp0 != sp1)) { op = 4;}
    }

    if (op != up_dp->p_op_idx) {
      sprintf_code(tmp_cnt,up_dp, up_dp->p_op_idx,up_dp->p_op_cnt);
      strncat(al_str,tmp_cnt,al_str_max);
      up_dp->p_op_idx = op;
      up_dp->p_op_cnt = 1;
    }
    else {
      up_dp->p_op_cnt++;
    }
  }
}

/* build an array of match/ins/del - length strings */
/* 5-June-2014 - modified to split "match" encoding into identical (=)
   and mismatch (X)

   To support domain-based scoring, this function iterates through
   every aligned position, including insertions and deletions (which
   are encoded as runs).
 */

int calc_code(const unsigned char *aa0, int n0,
	      const unsigned char *aa1, int n1,
	      struct a_struct *aln,
	      struct a_res_str *a_res,
	      struct pstruct *ppst,
	      char *al_str, int al_str_n, 
	      const unsigned char *ann_arr, 
	      const unsigned char *aa0a,
	      const struct annot_str *annot0_p,
	      const unsigned char *aa1a,
	      const struct annot_str *annot1_p,
	      struct dyn_string_str *ann_code_dyn,
	      int *score_delta,
	      struct f_struct *f_str,
	      void *pstat_void,
	      int display_code)
{
  int i0, i1;
  int op, lenc;
  int p_op, op_cnt;
  int match, p_match, match_cnt;
  int mis, p_mis, mis_cnt;
  const unsigned char *aa1p;
  char sp0, sp1;
  struct update_code_str *update_data_p;
  unsigned char *sq;
  int *rp;
  int have_ann=0;
  char ann_ch0, ann_ch1;
  char tmp_astr[MAX_STR];
  int sim_code, t_spa;
  char *sim_sym = aln_map_sym[MX_ACC];
  int aa0c, aa1c, itmp;
  int show_code, annot_fmt, start_flag;

  /* variables for variant changes, regions */
  int *aa0_pam2_p;
  void *annot_stack;
  struct annot_entry **s_annot0_arr_p;
  struct annot_entry **s_annot1_arr_p;
  int i0_annot, i1_annot, v_delta, v_tmp;
  long i0_left_end, i1_left_end;
  int d1_score, d1_ident, d1_alen;
  int d0_score, d0_ident, d0_alen;
  struct domfeat_link *left_domain_list1, *left_domain_list0;
  int have_push_features;
  long q_offset, l_offset;

  *score_delta = 0;
  i0_left_end = i1_left_end = -1;
  left_domain_list0 = left_domain_list1 = NULL;
  d1_score = d1_ident = d1_alen = 0;
  d0_score = d0_ident = d0_alen = 0;

  show_code = (display_code & (SHOW_CODE_MASK+SHOW_CODE_EXT));	/* see defs.h; SHOW_CODE_ALIGN=2,_CIGAR=3,_CIGAR_EXT=4 */
  annot_fmt = 2;
  if (display_code & SHOW_ANNOT_FULL) {
    annot_fmt = 1;
  }
  
  if (aa0a != NULL && aa1a != NULL) { have_ann = 2;}
  else if (aa0a != NULL || aa1a != NULL) { have_ann = 1;}
  else {have_ann = 0;}

  if (ppst->ext_sq_set) { sq = ppst->sqx;}
  else {sq = ppst->sq;}

#ifndef TFASTA
  aa1p = aa1;
#else
  aa1p = f_str->aa1x;
#endif

  rp = a_res->res;
  lenc = aln->nident = aln->nmismatch = aln->nsim = aln->npos = aln->ngap_q = aln->ngap_l = aln->nfs = 0;

  update_data_p = init_update_data(show_code);

  i0 = a_res->min0;
  i1 = a_res->min1;

  q_offset = aln->q_offset;
  l_offset = aln->l_offset;

  op = p_op = 0;
  op_cnt = match_cnt = mis_cnt = 0;
  start_flag = 1;

  v_delta = 0;
  i0_annot = i1_annot = 0;
  annot_stack = NULL;
  s_annot0_arr_p = s_annot1_arr_p = NULL;
  if (have_ann) {
    have_push_features = 0;
    if (annot0_p || annot1_p) annot_stack = init_stack(64,64);

    if (annot1_p && annot1_p->n_annot > 0) {
      s_annot1_arr_p = annot1_p->s_annot_arr_p;
      while (i1_annot < annot1_p->n_annot) {
	if (s_annot1_arr_p[i1_annot]->pos >= i1 + l_offset) {break;}
	if (s_annot1_arr_p[i1_annot]->end < i1 + l_offset) {i1_annot++; continue;}

	if (s_annot1_arr_p[i1_annot]->label == '-') {
	  process_annot_match(&itmp, aa0_pam2_p, l_offset+seq_pos(i1,aln->llrev,0), q_offset + seq_pos(i0,aln->qlrev,0),
			      &sp1, NULL, sq, s_annot1_arr_p[i1_annot], NULL, 
			      annot_stack, &have_push_features, &v_delta,
			      &d1_score, &d1_ident, &d1_alen, &left_domain_list1, &i1_left_end, 0);
	}
	i1_annot++;
      }
    }

    if (annot0_p && annot0_p->n_annot > 0) {
      s_annot0_arr_p = annot0_p->s_annot_arr_p;
      while (i0_annot < annot0_p->n_annot && s_annot0_arr_p[i0_annot]->pos < i0+q_offset) {
	if (s_annot0_arr_p[i0_annot]->pos >= i0 + q_offset) {break;}
	if (s_annot0_arr_p[i0_annot]->end < i0 + q_offset) {i0_annot++; continue;}

	if (s_annot0_arr_p[i0_annot]->label == '-') {
	  process_annot_match(&itmp, aa0_pam2_p, q_offset+seq_pos(i0,aln->qlrev,0), l_offset + seq_pos(i1,aln->llrev,0),
			      &sp0, NULL, sq, s_annot0_arr_p[i0_annot], NULL, 
			      annot_stack, &have_push_features, &v_delta,
			      &d0_score, &d0_ident, &d0_alen, &left_domain_list0, &i0_left_end, 0);
	}
	i0_annot++;
      }
    }
  }

  while (i0 < a_res->max0 || i1 < a_res->max1) {
    /* match/mismatch (aligned residues */
    /* here, op is the "current" encoding, and *rp is the next one */

    /* there is an op==0 for every aligned/matched residue.
       insertions in aa0 > 0,
       deletions < 0

       for enhanced CIGAR, need run lengths for identities, mismatches
     */

    if (ppst->pam_pssm) {aa0_pam2_p = ppst->pam2p[0][i0];}
    else {aa0_pam2_p = ppst->pam2[0][aa0[i0]];}

    if (op == 0 && *rp == 0) {
      aa0c = aa0[i0];
      aa1c = aa1p[i1];
      itmp = aa0_pam2_p[aa1c];
      sp0 = sq[aa0c];
      sp1 = sq[aa1c];

      /* variant annot1_p annotations can cause substitution */
      if (s_annot1_arr_p) {
	if (i1+l_offset == s_annot1_arr_p[i1_annot]->pos || i1+l_offset == i1_left_end) {

	  i1_annot = next_annot_match(&itmp, aa0_pam2_p, l_offset+seq_pos(i1,aln->llrev,0),
				      q_offset+seq_pos(i0,aln->qlrev,0), &sp1, NULL, sq,
				      i1_annot, annot1_p->n_annot, s_annot1_arr_p,
				      NULL, annot_stack, &have_push_features, &v_delta,
				      &d1_score, &d1_ident, &d1_alen, &left_domain_list1, &i1_left_end, 0);

	  if (sq[aa1c] != sp1) {
	    t_spa = align_type(itmp, sp0, sp1, ppst->nt_align, NULL, ppst->pam_x_id_sim);
	    comment_var(q_offset+seq_pos(i0,aln->qlrev,0), sp0,
			l_offset+seq_pos(i1,aln->llrev,0), sp1,
			sq[aa1c], sim_sym[t_spa], NULL,
			ann_code_dyn, 1, annot_fmt);
	  }
	}
	d1_score += itmp;
      }

      if (s_annot0_arr_p) {
	if (i0+q_offset == s_annot0_arr_p[i0_annot]->pos || i0+q_offset == i0_left_end) {

	  i0_annot = next_annot_match(&itmp, aa0_pam2_p, q_offset+seq_pos(i0,aln->qlrev,0),
				      l_offset+seq_pos(i1,aln->llrev,0), &sp0, NULL, sq,
				      i0_annot, annot0_p->n_annot, s_annot0_arr_p,
				      NULL, annot_stack, &have_push_features, &v_delta,
				      &d0_score, &d0_ident, &d0_alen, &left_domain_list0, &i0_left_end, 0);

	  /* check for sequence change from variant */
	  if (sq[aa0c] != sp0) {
	    t_spa = align_type(itmp, sp0, sp1, ppst->nt_align, NULL, ppst->pam_x_id_sim);
	    comment_var(q_offset+seq_pos(i0,aln->qlrev,0), sp0,
			l_offset+seq_pos(i1,aln->llrev,0), sp1,
			sq[aa0c], sim_sym[t_spa], NULL,
			ann_code_dyn, 0, annot_fmt);
	  }
	}
	d0_score += itmp;
      }

      d0_alen++;
      d1_alen++;
      if ((sim_code == align_type(itmp, sp0, sp1, ppst->nt_align, aln, ppst->pam_x_id_sim)) == M_IDENT) {
	d0_ident++;
	d1_ident++;
      }

      update_code(al_str, al_str_n-strlen(al_str), update_data_p, op, sim_code, sp0, sp1);

      /* update op to the next encoded position */
      op = *rp++;
      lenc++;
      
      /* check for an annotation */
      if (have_ann) {
	ann_ch0 = ann_ch1 = '\0';
	/* conventional annotations */
	if (have_ann == 2 && (ann_arr[aa0a[i0]] != ' ' || ann_arr[aa1a[i1]] != ' ')) {
	  ann_ch0 = ann_arr[aa0a[i0]];
	  if (ann_ch0 == ' ') ann_ch0 = 'X';
	  ann_ch1 = ann_arr[aa1a[i1]];
	  if (ann_ch1 == ' ') ann_ch1 = 'X';
	}
	else if (aa0a != NULL && ann_arr[aa0a[i0]]!=' ') {
	  ann_ch0 = ann_arr[aa0a[i0]];
	  ann_ch1 = 'X';
	}
	else if (aa1a != NULL && ann_arr[aa1a[i1]]!=' ') {
	  ann_ch0 = 'X';
	  ann_ch1 = ann_arr[aa1a[i1]];
	}

	/* ann_ch0 only works below because ann_ch0=='X' if ann_ch1 */
	if ( ann_ch0 && !(ann_ch1 == '['  || ann_ch1 == ']' || ann_ch0 == '[' || ann_ch0 == ']')) {
	  sprintf(tmp_astr, "|%c%c:%ld%c%c%ld%c",
		  ann_ch0,ann_ch1, q_offset+seq_pos(i0,aln->qlrev,0)+1,sp0,
		  sim_sym[sim_code], l_offset+seq_pos(i1,aln->llrev,0)+1,sp1);
	  /* SAFE_STRNCAT(ann_code_s, tmp_astr, n_ann_code_s); */
	  dyn_strcat(ann_code_dyn, tmp_astr);
	}
	  
	if ((s_annot1_arr_p || s_annot0_arr_p) && have_push_features) {
	  display_push_features(annot_stack, ann_code_dyn,
				q_offset+seq_pos(i0,aln->qlrev,0), sp0,
				l_offset+seq_pos(i1,aln->llrev,0), sp1,
				sim_sym[sim_code], 
				a_res->rst.score[ppst->score_ix], a_res->rst.comp, n0, n1, pstat_void, annot_fmt);
	  have_push_features = 0;
	}
      }
      i0++; i1++;
    }
    else {	/* not in match run, in a gap */
      if (op == 0) {
	/* at a transition from match (previous) to indel (current) */
	d1_score += ppst->gdelval;
	d0_score += ppst->gdelval;
	op = *rp++;
      }
      d1_score += ppst->ggapval; d1_alen++;
      d0_score += ppst->ggapval; d0_alen++;

      if (op > 0) {
	update_code(al_str, al_str_n-strlen(al_str), update_data_p, 2, sim_code,'-','-');

	if (s_annot1_arr_p) {
	  if (i1+l_offset == s_annot1_arr_p[i1_annot]->pos || i1+l_offset == i1_left_end) {
	    i1_annot = next_annot_match(&itmp, aa0_pam2_p, l_offset+seq_pos(i1,aln->llrev,0),
					q_offset+seq_pos(i0,aln->qlrev,0), &sp1, NULL, sq, 
					i1_annot, annot1_p->n_annot, s_annot1_arr_p,
					NULL, annot_stack, &have_push_features, &v_delta,
					&d1_score, &d1_ident, &d1_alen, &left_domain_list1, &i1_left_end, 0);

	    if (have_push_features) {
	      display_push_features(annot_stack, ann_code_dyn,
				    q_offset+seq_pos(i0,aln->qlrev,0), sp0,
				    l_offset+seq_pos(i1,aln->llrev,0), sp1,
				    sim_sym[sim_code],
				    a_res->rst.score[ppst->score_ix], a_res->rst.comp, n0, n1, pstat_void, annot_fmt);
	      have_push_features = 0;
	    }
	  }
	}
	op--; lenc++; i1++; aln->ngap_q++;
      }
      else {  /* (op < 0) */

	update_code(al_str, al_str_n-strlen(al_str), update_data_p, 1, sim_code,'-','-');

	if (s_annot0_arr_p) {
	  if (i0+q_offset == s_annot0_arr_p[i0_annot]->pos || i0+q_offset == i0_left_end) {

	    i0_annot = next_annot_match(&itmp, aa0_pam2_p, q_offset+seq_pos(i0,aln->qlrev,0),
					l_offset+seq_pos(i1,aln->llrev,0), &sp0, NULL, sq, 
					i0_annot, annot0_p->n_annot, s_annot0_arr_p,
					NULL, annot_stack, &have_push_features, &v_delta,
					&d0_score, &d0_ident, &d0_alen, &left_domain_list0, &i0_left_end, 0);

	    if (have_push_features) {
	      display_push_features(annot_stack, ann_code_dyn,
				    q_offset+seq_pos(i0,aln->qlrev,0), sp0,
				    l_offset+seq_pos(i1,aln->llrev,0), sp1,
				    sim_sym[sim_code],
				    a_res->rst.score[ppst->score_ix], a_res->rst.comp, n0, n1, pstat_void, annot_fmt);
	      have_push_features = 0;
	    }
	  }
	}
	op++; lenc++; i0++; aln->ngap_l++;
      }
    }
  }

  /* all done, clean things up */

  close_update_data(al_str, al_str_n-strlen(al_str), update_data_p);

  if (have_ann) {
    have_push_features = 0;
    /* also check for regions after alignment */

    if (s_annot1_arr_p && i1_left_end > 0) {
      close_annot_match(-1, annot_stack, &have_push_features,
			&d1_score, &d1_ident, &d1_alen, &left_domain_list1, &i1_left_end,
			0);
    }
    if (s_annot0_arr_p && i0_left_end > 0) {
      close_annot_match(-1, annot_stack, &have_push_features,
			&d0_score, &d0_ident, &d0_alen, &left_domain_list0, &i0_left_end,
			0);
    }

    if (have_push_features) {
      display_push_features(annot_stack, ann_code_dyn,
			    q_offset+a_res->max0-1, sp0,
			    l_offset+a_res->max1-1, sp1,
			    sim_sym[sim_code],
			    a_res->rst.score[ppst->score_ix], a_res->rst.comp, n0, n1, pstat_void, annot_fmt);
    }
  }

  if (annot0_p || annot1_p) free_stack(annot_stack);

  *score_delta = v_delta;
  return lenc;
}

int calc_id(const unsigned char *aa0, int n0,
	    const unsigned char *aa1, int n1,
	    struct a_struct *aln, 
	    struct a_res_str *a_res,
	    struct pstruct *ppst,
	    const struct annot_str *annot0_p,
	    const struct annot_str *annot1_p,
	    int *score_delta,
	    struct dyn_string_str *annot_var_dyn,
	    struct f_struct *f_str)
{
  int i0, i1, nn1;
  int op, lenc;
  char sp0, sp1;
  char tmp_str[MAX_SSTR];
  const unsigned char *aa1p;
  int *rp;
  unsigned char *sq;

  /* variables for variant changes */
  int *aa0_pam2_p;
  struct annot_entry **s_annot0_arr_p;
  struct annot_entry **s_annot1_arr_p;
  int  itmp, i0_annot, i1_annot, v_delta, v_tmp;
  long q_offset, l_offset;
  long i0_left_end, i1_left_end;
  int d1_score, d1_ident, d1_alen;
  int d0_score, d0_ident, d0_alen;
  struct domfeat_link *left_domain_list1, *left_domain_list0;
  struct domfeat_link *this_dom, *next_dom;

  left_domain_list1 = left_domain_list0 = NULL;
  
  *score_delta = 0;
  i0_left_end = i1_left_end = -1;
  left_domain_list0 = left_domain_list1 = NULL;
  
  NULL_dyn_string(annot_var_dyn);

  if (ppst->ext_sq_set) { sq = ppst->sqx; }
  else { sq = ppst->sq; }

#ifndef TFASTA
  aa1p = aa1;
  nn1 = n1;
#else
  aa1p = f_str->aa1x;
  nn1 = f_str->n10;
#endif

  aln->amin0 = a_res->min0;
  aln->amax0 = a_res->max0;
  aln->amin1 = a_res->min1;
  aln->amax1 = a_res->max1;
  aln->calc_last_set = 1;

  q_offset = aln->q_offset;
  l_offset = aln->l_offset;

  rp = a_res->res;
  lenc = aln->nident = aln->nmismatch = aln->nsim = aln->ngap_q = aln->ngap_l = aln->nfs = op = 0;
  i0 = a_res->min0;
  i1 = a_res->min1;

  v_delta = 0;
  i0_annot = i1_annot = 0;

  d1_score = d1_ident = d1_alen = 0;
  d0_score = d0_ident = d0_alen = 0;

  if (annot1_p && annot1_p->n_annot > 0) s_annot1_arr_p = annot1_p->s_annot_arr_p;
  else s_annot1_arr_p = NULL;
  if (annot0_p && annot0_p->n_annot > 0) s_annot0_arr_p = annot0_p->s_annot_arr_p;
  else s_annot0_arr_p = NULL;

  while (i0 < a_res->max0 || i1 < a_res->max1) {

    if (ppst->pam_pssm) {
      aa0_pam2_p = ppst->pam2p[0][i0];
    }
    else {
      aa0_pam2_p = ppst->pam2[0][aa0[i0]];
    }

    if (op == 0 && *rp == 0) {
      /* op==0 -> we are in a match run, and current code is a match */
      op = *rp++;
      lenc++;

      itmp = ppst->pam2[0][aa0[i0]][aa1p[i1]];
      sp0 = sq[aa0[i0]];
      sp1 = sq[aa1p[i1]];

      if (s_annot1_arr_p && (i1 + l_offset == s_annot1_arr_p[i1_annot]->pos || i1+l_offset == i1_left_end)) {
	i1_annot = next_annot_match(&itmp, aa0_pam2_p, l_offset+seq_pos(i1,aln->llrev,0),
				    q_offset+seq_pos(i0,aln->qlrev,0), &sp1, NULL, sq,
				    i1_annot, annot1_p->n_annot, s_annot1_arr_p,
				    NULL, NULL, NULL, &v_delta, 
				    &d1_score, &d1_ident, &d1_alen, &left_domain_list1, &i1_left_end, itmp);

	/* must be out of the loop to capture the last value */
	if (sq[aa1p[i1]] != sp1) {
	  sprintf(tmp_str,"%c%d%c;",sq[aa1p[i1]],i1+1,sp1);
	  /*   SAFE_STRNCAT(annot_var_s,tmp_str,n_annot_var_s); */
	  dyn_strcat(annot_var_dyn, tmp_str);
	}
	d1_score += itmp;
      }

      if (s_annot0_arr_p && (i0 + q_offset == s_annot0_arr_p[i0_annot]->pos || i0+q_offset == i0_left_end)) {
	i0_annot = next_annot_match(&itmp, aa0_pam2_p, q_offset+seq_pos(i0,aln->qlrev,0),
				    l_offset+seq_pos(i1,aln->llrev,0), &sp0, NULL, sq,
				    i0_annot, annot0_p->n_annot, s_annot0_arr_p,
				    NULL, NULL, NULL, &v_delta, 
				    &d0_score, &d0_ident, &d0_alen, &left_domain_list0, &i0_left_end, itmp);



	if (sq[aa0[i0]] != sp0) {
	  sprintf(tmp_str,"q%c%d%c;",sq[aa0[i0]],i0+1,sp0);
	  /*  SAFE_STRNCAT(annot_var_s,tmp_str,n_annot_var_s); */
	  dyn_strcat(annot_var_dyn, tmp_str);
	}
	d0_score += itmp;
      }

      /* updates nident, nsim, npos */
      d0_alen++;
      d1_alen++;
      if (align_type(itmp, sp0, sp1, ppst->nt_align, aln, ppst->pam_x_id_sim) == M_IDENT) {
	d0_ident++;
	d1_ident++;
      }

      i0++; i1++;
    }
    else {
      if (op==0) op = *rp++;
      if (op>0) {	/* inserts in seq0 */
	op--; lenc++; i1++; aln->ngap_q++;
      }
      else {		/* inserts in seq 1 */
	op++; lenc++; i0++; aln->ngap_l++;
      }
    }
  }
  *score_delta = v_delta;

  return lenc;
}
