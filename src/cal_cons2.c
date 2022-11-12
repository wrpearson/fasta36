/* cal_cons.c - routines for printing alignments for
   fasta, ssearch, ggsearch, glsearch  */

/* $Id: cal_cons.c 1280 2014-08-21 00:47:55Z wrp $ */

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
  int btop_enc;
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
update_code(struct dyn_string_str *align_code_dyn,
	    struct update_code_str *update_data, int op, 
	    int sim_code, unsigned char sp0, unsigned char sp1);

static void
close_update_data(struct dyn_string_str *align_code_dyn,
		  struct update_code_str *update_data);

extern void aancpy(char *to, char *from, int count, struct pstruct *ppst);
extern void *init_stack(int, int);
extern void push_stack(void *, void *);
extern void *pop_stack(void *);
extern void *free_stack(void *);
extern struct domfeat_data * init_domfeat_data(const struct annot_str *annot_p);

/* returns M_NEG, M_ZERO, M_POS, M_IDENT, M_DEL (a_mark.h)
   updates *aln->nsim, npos, nident, nmismatch */
extern int
align_type(int score, char sp0, char sp1, int nt_align, struct a_struct *aln, int pam_x_id_sim);

extern void	/* in compacc2e.c */
process_annot_match(int *itmp, int *pam2aa0v, 
		    long ip, long ia, char *sp1, char *sp1a, const unsigned char *sq,
		    struct annot_entry *annot_arr_p, int n_domains, char **ann_comment,
		    void *annot_stack, int *have_push_features_p, int *v_delta,
		    int *d_score_p, int *d_ident_p, int *d_alen_p, int *d_gaplen_p,
		    struct domfeat_data **left_domain_head_p,
		    struct domfeat_data *left_domain_p,
		    long *left_end_p, int init_score);

extern int	/* in compacc2e.c */
next_annot_match(int *itmp, int *pam2aa0v, 
		 long ip, long ia, char *sp1, char *sp1a, const unsigned char *sq,
		 int i_annot, int n_annot, struct annot_entry **annot_arr, char **ann_comment,
		 void *annot_stack, int *have_push_features_p, int *v_delta,
		 int *d_score_p, int *d_ident_p, int *d_alen_p, int *d_gaplen_p,
		  struct domfeat_data **left_domain_head_p, struct domfeat_data *left_domain_p,
		 long *left_domain_end, int init_score);

extern void	/* in compacc2e.c */
close_annot_match (int ia, void *annot_stack, int *have_push_features_p,
		   int *d_score_p, int *d_ident_p, int *d_alen_p, int *d_gaplen_p,
		   struct domfeat_data **left_domain_p,
		   long *left_end_p, int init_score);

extern void	/* in compacc2e.c */
comment_var(long i0, char sp0, long i1, char sp1, char o_sp1, char sim_char,
	    const char *ann_comment, struct dyn_string_str *annot_var_dyn,
	    int target, int d_type);

extern void		/* in compacc2e.c */
display_push_features(void *annot_stack, struct dyn_string_str *annot_var_dyn,
		      long i0_pos, char sp0, long i1_pos, char sp1, char sym, 
		      int score, double comp, int sw_score, int n0, int n1,
		      void *pstat_void, int d_type);

#define DP_FULL_FMT 1	/* Region: score: bits: id: ... */

extern int seq_pos(int pos, int rev, int off);

/* values of calc_func_mode */
#define CALC_CONS 1    /* standard calculate consensus -- make label */
#define CALC_CODE 2    /* only calculate annotation code */
#define CALC_ID   3    /* only calculate identity -- no codes, so no display */
#define CALC_ID_DOM   4  /* like CALC_ID, but also looks at identity in domains?? (wrp 2022-07-12) */

int
pre_fill_cons(const unsigned char *aa0, const unsigned char *aa1p,
	      struct a_res_str *a_res,
	      struct pstruct *ppst,
	      struct a_struct *aln,
	      char *seqc0, char *seqc1, char *seqca, int *smins,
	      char *seqc0a, char *seqc1a) {
  int mins;

  /* will we show all the start ?*/
  if (min(a_res->min0,a_res->min1)<aln->llen || aln->showall==1)
    if (a_res->min0 >= a_res->min1) {              /* aa0 extends more to left */
      *smins=0;
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
      *smins=0;
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
    *smins=mins;
    aln->smin0=a_res->min0 - *smins;
    aln->smin1=a_res->min1 - *smins;
    aancpy(seqc0,(char *)aa0+a_res->min0-mins,mins,ppst);
    aancpy(seqc1,(char *)aa1p+a_res->min1-mins,mins,ppst);
  }
  /* set the alignment code to zero for context */
  memset(seqca,0,mins);
  if (seqc0a) {
    memset(seqc0a,' ',mins);
    memset(seqc1a,' ',mins);
  }
  return mins;
}

int 
post_fill_cons(const unsigned char *aa0, int n0,
	       const unsigned char *aa1p, int nn1,
	       struct a_res_str *a_res,
	       struct pstruct *ppst,
	       int mins, int lenc,
	       struct a_struct *aln,
	       char *seqc0, char *seqc1,
	       char *seqc0a, char *seqc1a) {

  int ns, nd, itmp;

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
  if (seqc0a) {
    memset(seqc0a+mins+lenc,' ',nd);
    memset(seqc1a+mins+lenc,' ',nd);
    /*
      ntmp = nd-(n0-a_res->max0);
      if (ntmp > 0) memset(seqc0a+mins+lenc+n0-a_res->max0,' ',ntmp);
      ntmp = nd-(nn1-a_res->max1);
      if (ntmp > 0) memset(seqc1a+mins+lenc+nn1-a_res->max1,' ',ntmp);
    */
  }
  return nd;
}

/* add_annot_code: adds annotation codes to struct dyn_string_str ann_code_dyn */
void
add_annot_code(int have_ann, char sp0, char sp1, 
	       char ann_aa0_i0, char ann_aa1_i1,
	       long q_off_pos, long l_off_pos, char sim_sym_code,
	       struct dyn_string_str *ann_code_dyn)
{
  char ann_ch0, ann_ch1;
  char tmp_astr[MAX_STR];

  ann_ch0 = ann_ch1 = '\0';

  /* conventional annotations */
  if (have_ann == 3 && (ann_aa0_i0 != ' ' || ann_aa1_i1 != ' ')) {
    ann_ch0 = ann_aa0_i0;
    if (ann_ch0 == ' ') ann_ch0 = 'X';
    ann_ch1 = ann_aa1_i1;
    if (ann_ch1 == ' ') ann_ch1 = 'X';
  }
  else if ((have_ann&2)==0 && ann_aa0_i0 != ' ') {
    ann_ch0 = ann_aa0_i0;
    ann_ch1 = 'X';
  }
  else if ((have_ann&1)==0 && ann_aa1_i1 != ' ') {
    ann_ch0 = 'X';
    ann_ch1 = ann_aa1_i1;
  }

  /* ann_ch0 only works below because ann_ch0=='X' if ann_ch1 */
  /*
  if ( ann_ch0 && !(ann_ch1 == '['  || ann_ch1 == ']' || ann_ch0 == '[' || ann_ch0 == ']')) {
    sprintf(tmp_astr, "|%c%c:%ld%c%c%ld%c",
	    ann_ch0,ann_ch1, q_offset+seq_pos(i0,aln->qlrev,0)+1,sp0,
	    sim_sym[sim_code], l_offset+seq_pos(i1,aln->llrev,0)+1,sp1);
  */
  /* excluding 'V' annotations because they were added by comment_var if they were used */

  if ( ann_ch0 && !(ann_ch1 == '['  || ann_ch1 == ']' || ann_ch0 == '[' || ann_ch0 == ']' || ann_ch1 == 'V' || ann_ch0 == 'V')) {
    sprintf(tmp_astr, "|%c%c:%ld%c%c%ld%c",
	    ann_ch0,ann_ch1, q_off_pos+1,sp0,
	    sim_sym_code, l_off_pos+1,sp1);
    dyn_strcat(ann_code_dyn, tmp_astr);
  }
}

/* calc_cons_u - combines calc_cons_a/calc_code/ calc_id */
int
calc_cons_u( /* inputs */
	    const unsigned char *aa0, int n0,
	    const unsigned char *aa1, int n1,
	    struct a_res_str *a_res,	/* alignment encoding */
	    struct pstruct *ppst,
	    struct f_struct *f_str,
	    void *pstat_void,
	    /* annotation stuff */
	    const unsigned char *ann_arr,
	    const unsigned char *aa0a, const struct annot_str *annot0_p,
	    const unsigned char *aa1a, const struct annot_str *annot1_p,
	    int calc_func_mode, 	/* CALC_CONS, CALC_CODE, CALC_ID */
	    int display_code, 	/* used only by CALC_CODE */
	    /* outputs */
	    int *nc,
	    char *seqc0, char *seqc1, char *seqca, int *cumm_seq_score,
	    char *seqc0a, char *seqc1a,
	    struct a_struct *aln,
	    int *score_delta, 
	    struct dyn_string_str *annot_var_dyn,
	    struct dyn_string_str *align_code_dyn)
{
  int i0, i1, nn1;
  int op, lenc, nd, ns, itmp;
  const unsigned char *aa1p;
  char *sp0_p, *sp0a_p, *sp1_p, *sp1a_p, *spa_p, t_spa;
  char sp0_c, sp1_c, spa_c;	/* used for CALC_ID, CALC_CODE */
  char sp0a_c, sp1a_c;		/* used for CALC_CODE */
  char tmp_str[MAX_SSTR];
  int *i_spa;
  const unsigned char *sq;
  int *rp;
  int smins, mins, ntmp;
  int have_ann;
  void *annot_stack = NULL;
  struct update_code_str *update_data_p;

  /* variables for variant changes */
  int *aa0_pam2_p;
  char *sim_sym = aln_map_sym[5];
  struct annot_entry **s_annot0_arr_p, **s_annot1_arr_p;

  char *ann_comment;
  int i0_annot, i1_annot;	/* i0_annot, i1_annot, count through
				   the list of annotations */
  long i0_left_end, i1_left_end; /* left-most coordinate of domain end */

  int show_code, annot_fmt;

  int v_delta, v_tmp;
  int d1_score, d1_ident, d1_alen, d1_gaplen;
  int d0_score, d0_ident, d0_alen, d0_gaplen;
  int have_push_features;
  int *have_push_features_p;

  /* struct domfeat_data is used to capture score, coordinate, and
     identity information for possibly overlapping sub-alignment
     scores -- each domfeat_data entry is associated with an
     annot_p->annot_arr_p entry */
  struct domfeat_data *left_domain_head1, *left_domain_head0;
  struct domfeat_data *left_domain_list1, *left_domain_list0;

  /* variables for handling coordinate offsets */
  long q_offset, l_offset;
  long i0_off, i1_off;

  if (ppst->ext_sq_set) { sq = ppst->sqx;}
  else {sq = ppst->sq;}

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

#ifndef LCAL_CONS	/* use for local context */
  if (calc_func_mode == CALC_CONS) {
    mins = pre_fill_cons(aa0, aa1p, a_res,ppst,  aln, seqc0, seqc1, seqca, &smins,
			 seqc0a, seqc1a);
  }
#else		/* no flanking context */
  smins = mins = 0;
  aln->smin0=a_res->min0;
  aln->smin1=a_res->min1;
#endif

  /* now get the middle */
  have_ann = 0;	/* default no annotation */
  left_domain_head0 = left_domain_head1 = NULL;
  left_domain_list0 = left_domain_list1 = NULL;

  /* have_ann encodes which sequences are annotated */
  if ((annot0_p && annot0_p->n_annot > 0) || (aa0a != NULL)) { have_ann |= 1;}
  if ((annot1_p && annot1_p->n_annot > 0) || (aa1a != NULL)) { have_ann |= 2;}

  if (calc_func_mode == CALC_CONS) {
    spa_p = seqca+mins;	/* pointer to alignment symbol */
    if (cumm_seq_score) i_spa = cumm_seq_score+mins;	/* set index for cumm_seq_score */
    sp0_p = seqc0+mins;
    sp1_p = seqc1+mins;
    /*     have_ann = (seqc0a != NULL);  */
    annot_fmt = DP_FULL_FMT;
  }
  else if (calc_func_mode == CALC_ID || calc_func_mode == CALC_ID_DOM) {
    spa_p = &spa_c;
    sp0_p = &sp0_c;
    sp1_p = &sp1_c;
    /* does not require aa0a/aa1a, only for variants */
    /*     have_ann = ((annot1_p && annot1_p->n_annot > 0) || (annot0_p && annot0_p->n_annot > 0)); */
    annot_fmt = 3;
  }
  else if (calc_func_mode == CALC_CODE) {
    spa_p = &spa_c;
    sp0_p = &sp0_c;
    sp1_p = &sp1_c;

    show_code = (display_code & (SHOW_CODE_MASK+SHOW_CODE_EXT)); /* see defs.h; SHOW_CODE_ALIGN=4,_CIGAR=8,_CIGAR_EXT=24, _BTOP_EXT=16 */
    annot_fmt = 2;
    if (display_code & SHOW_ANNOT_FULL) {
      annot_fmt = 1;
    }

    update_data_p = init_update_data(show_code);
  }
  else {
    fprintf(stderr,"*** ERROR [%s:%d] --- cal_cons_u() invalid calc_func_mode: %d\n",
	    __FILE__, __LINE__, calc_func_mode);
    exit(1);
  }

  have_push_features=0;

  if (have_ann) {  /* initialize annotation variables */
    if (calc_func_mode == CALC_CONS) {
      sp0a_p = seqc0a+mins;
      sp1a_p = seqc1a+mins;
      annot_stack = init_stack(64,64);
      have_push_features_p = &have_push_features;
    }
    else if (calc_func_mode == CALC_ID || calc_func_mode == CALC_ID_DOM) {
      sp0a_p = &sp0a_c;
      sp1a_p = &sp1a_c;
      have_push_features_p = &have_push_features;
      /*      ann_comment = NULL; */
      annot_stack = init_stack(64,64);
    }
    else if (calc_func_mode == CALC_CODE) {
      annot_stack = init_stack(64,64);
      sp0a_p = &sp0a_c;
      sp1a_p = &sp1a_c;
      have_push_features_p = &have_push_features;
    }

    *score_delta = 0;
    i0_left_end = i1_left_end = -1;
    NULL_dyn_string(annot_var_dyn);
  }
  /* always initialize, updated with no annotations */
  d1_score = d1_ident = d1_alen = d1_gaplen = 0;
  d0_score = d0_ident = d0_alen = d0_gaplen = 0;

  lenc = aln->nident = aln->nmismatch =
    aln->npos = aln->nsim = aln->ngap_q = aln->ngap_l = aln->nfs = op = 0;

  i0 = a_res->min0;	/* start in aa0[] */
  i1 = a_res->min1;	/* start in aa1[] */

  /* handle region annotations outside alignment */
  v_delta = 0;
  i0_annot = i1_annot = 0;
  s_annot0_arr_p = s_annot1_arr_p = NULL;
  if (have_ann) {
    i1_off = seq_pos(i1, aln->llrev,0) + l_offset;
    i0_off = seq_pos(i0, aln->qlrev,0) + q_offset;

    if (annot1_p && annot1_p->n_annot > 0) {

	left_domain_list1 = init_domfeat_data(annot1_p);
	s_annot1_arr_p = annot1_p->s_annot_arr_p;

	while (i1_annot < annot1_p->n_annot) {
	  if (s_annot1_arr_p[i1_annot]->pos >= i1_off) {break;}
	  if (s_annot1_arr_p[i1_annot]->end <= i1_off) {i1_annot++; continue;}

	  if (s_annot1_arr_p[i1_annot]->label == '-') {
	    process_annot_match(&itmp, aa0_pam2_p, i1_off, i0_off,
				sp1_p, sp1a_p, sq, s_annot1_arr_p[i1_annot],  annot1_p->n_domains, &ann_comment, 
				annot_stack, have_push_features_p, &v_delta,
				&d1_score, &d1_ident, &d1_alen, &d1_gaplen,
				&left_domain_head1,
				&left_domain_list1[i1_annot], &i1_left_end, 0);
	  }
	  i1_annot++;
	}
    }

    /* do not need have_ann here, because domain only */
    if (annot0_p && annot0_p->n_annot>0) {

      if (calc_func_mode == CALC_CONS || calc_func_mode == CALC_CODE) {

	/* inefficient -- the same initiation is done for every
	  query/subj alignment, even though it is always the same --
	  should be done in the build_ares() loop */
	left_domain_list0 = init_domfeat_data(annot0_p);
	s_annot0_arr_p = annot0_p->s_annot_arr_p;

	while (i0_annot < annot0_p->n_annot) {
	  if (s_annot0_arr_p[i0_annot]->pos >= i0_off) {break;}
	  if (s_annot0_arr_p[i0_annot]->end <= i0_off) {i0_annot++; continue;}

	  if (s_annot0_arr_p[i0_annot]->label == '-') {
	    process_annot_match(&itmp, NULL, i0_off, i1_off,
				sp0_p, sp0a_p, sq, s_annot0_arr_p[i0_annot], annot0_p->n_domains,  &ann_comment, 
				annot_stack, have_push_features_p, &v_delta,
				&d0_score, &d0_ident, &d0_alen, &d0_gaplen,
				&left_domain_head0,
				&left_domain_list0[i0_annot], &i0_left_end, 0);
	  }
	  i0_annot++;
	}
      }
    }
  }
  /* done with domains starting before alignment */

  /* handle alignment encoding */
  rp = a_res->res;	/* alignment encoding array */
  while (i0 < a_res->max0 || i1 < a_res->max1) {
    /* match/mismatch (aligned residues */
    /* here, op is the "current" encoding, and *rp is the next one */
    if (op == 0 && *rp == 0) {
      op = *rp++;
      lenc++;

      if (ppst->pam_pssm) {aa0_pam2_p = ppst->pam2p[0][i0];}
      else {aa0_pam2_p = ppst->pam2[0][aa0[i0]];}

      itmp=aa0_pam2_p[aa1p[i1]];

      *sp0_p = sq[aa0[i0]];
      *sp1_p = sq[aa1p[i1]];

      if (have_ann) {
	have_push_features = 0;
	*sp0a_p = *sp1a_p = ' ';
	if (aa0a) *sp0a_p = ann_arr[aa0a[i0]];
	if (aa1a) *sp1a_p = ann_arr[aa1a[i1]];

	if (s_annot1_arr_p) {
	  if (i1+l_offset == s_annot1_arr_p[i1_annot]->pos) {

	    i1_annot = next_annot_match(&itmp, aa0_pam2_p, l_offset+seq_pos(i1,aln->llrev,0),
					q_offset + seq_pos(i0,aln->qlrev,0), sp1_p, sp1a_p, sq, 
					i1_annot, annot1_p->n_annot, s_annot1_arr_p,
					&ann_comment, annot_stack, have_push_features_p, &v_delta,
					&d1_score, &d1_ident, &d1_alen, &d1_gaplen,
					&left_domain_head1, left_domain_list1, &i1_left_end, 0);

	    /* must be out of the loop to capture the last value */
	    if (sq[aa1p[i1]] != *sp1_p) {
	      t_spa = align_type(itmp, *sp0_p, *sp1_p, ppst->nt_align, NULL, ppst->pam_x_id_sim);
	      if (calc_func_mode != CALC_ID && calc_func_mode != CALC_ID_DOM) {
		comment_var(q_offset+seq_pos(i0,aln->qlrev,0), *sp0_p,
			    l_offset+seq_pos(i1,aln->llrev,0), *sp1_p,
			    sq[aa1p[i1]], sim_sym[t_spa], ann_comment,
			    annot_var_dyn, 1, annot_fmt);
	      }
	      else {
		sprintf(tmp_str,"%c%d%c;",sq[aa1p[i1]],i1+1,*sp1_p);
		/*   SAFE_STRNCAT(annot_var_s,tmp_str,n_annot_var_s); */
		dyn_strcat(annot_var_dyn, tmp_str);
	      }
	    }
	  }
	  d1_score += itmp;
	}

	if (s_annot0_arr_p) {
	  if (i0 + q_offset == s_annot0_arr_p[i0_annot]->pos) {

	    i0_annot = next_annot_match(&itmp, aa0_pam2_p, q_offset+seq_pos(i0,aln->qlrev,0),
					l_offset + seq_pos(i1,aln->llrev,0), sp0_p, sp0a_p, sq, 
					i0_annot, annot0_p->n_annot, s_annot0_arr_p,
					&ann_comment, annot_stack, have_push_features_p, &v_delta,
					&d0_score, &d0_ident, &d0_alen, &d0_gaplen,
					&left_domain_head0, left_domain_list0, &i0_left_end, 0);

	    /* must be out of the loop to capture the last value */
	    if (sq[aa0[i0]] != *sp0_p) {
	      t_spa = align_type(itmp, *sp0_p, *sp1_p, ppst->nt_align, NULL, ppst->pam_x_id_sim);
	      if (calc_func_mode != CALC_ID && calc_func_mode != CALC_ID_DOM) {

	      comment_var(q_offset+seq_pos(i0,aln->qlrev,0), *sp0_p,
			  l_offset+seq_pos(i1,aln->llrev,0), *sp1_p,
			  sq[aa0[i0]], sim_sym[t_spa], ann_comment,
			  annot_var_dyn, 0, annot_fmt);
	      }
	      else {
		sprintf(tmp_str,"q%c%d%c;",sq[aa0[i0]],i0+1,*sp0_p);
		/*  SAFE_STRNCAT(annot_var_s,tmp_str,n_annot_var_s); */
		dyn_strcat(annot_var_dyn, tmp_str);
	      }
	    }
	  }
	  d0_score += itmp;
	}
	if (calc_func_mode == CALC_CONS) {sp0a_p++; sp1a_p++;}
      }

      *spa_p = align_type(itmp, *sp0_p, *sp1_p, ppst->nt_align, aln, ppst->pam_x_id_sim);

      d1_alen++;
      d0_alen++;
      if (*spa_p == M_IDENT) {
	d1_ident++;
	d0_ident++;
      }

      if (s_annot1_arr_p && (i1 + l_offset == i1_left_end)) {
	i1_annot = next_annot_match(&itmp, aa0_pam2_p, l_offset+seq_pos(i1,aln->llrev,0),
				    q_offset + seq_pos(i0,aln->qlrev,0), sp1_p, sp1a_p, sq, 
				    i1_annot, annot1_p->n_annot, s_annot1_arr_p,
				    &ann_comment, annot_stack, have_push_features_p, &v_delta,
				    &d1_score, &d1_ident, &d1_alen, &d1_gaplen,
				    &left_domain_head1, left_domain_list1, &i1_left_end, 0);
      }

      if (s_annot0_arr_p && (i0 + q_offset == i0_left_end)) {
	i0_annot = next_annot_match(&itmp, aa0_pam2_p, q_offset+seq_pos(i0,aln->qlrev,0),
				    l_offset + seq_pos(i1,aln->llrev,0), sp0_p, sp0a_p, sq, 
				    i0_annot, annot0_p->n_annot, s_annot0_arr_p,
				    &ann_comment, annot_stack, have_push_features_p, &v_delta,
				    &d0_score, &d0_ident, &d0_alen, &d0_gaplen,
				    &left_domain_head0, left_domain_list0, &i0_left_end, 0);
      }

      if (calc_func_mode == CALC_CODE) {
	update_code(align_code_dyn, update_data_p, op, *spa_p, *sp0_p, *sp1_p);
      }

      /* now we have done all the ?modified identity checks, display
	 potential site annotations */
      if (have_ann && calc_func_mode == CALC_CODE) {
	  add_annot_code(have_ann, *sp0_p, *sp1_p, *sp0a_p, *sp1a_p,
			 q_offset + seq_pos(i0,aln->qlrev,0), l_offset+seq_pos(i1,aln->llrev,0),
			 sim_sym[*spa_p], annot_var_dyn);
      }

      if (have_push_features && calc_func_mode != CALC_ID) {
	display_push_features(annot_stack, annot_var_dyn,
			      q_offset+seq_pos(i0,aln->qlrev,0), *sp0_p,
			      l_offset+seq_pos(i1,aln->llrev,0), *sp1_p,
			      sim_sym[*spa_p],
			      a_res->rst.score[ppst->score_ix], a_res->rst.comp, a_res->sw_score,
			      n0, n1,
			      pstat_void, annot_fmt);
	have_push_features = 0;
      }

      if (cumm_seq_score) *i_spa++ = itmp;	/* update cummulative score */
      i0++; i1++;
      if (calc_func_mode == CALC_CONS) {
	sp0_p++; sp1_p++; spa_p++;
      }
    }
    else {	/* indel */
      /* include all calc_func_mode's, because i_annot must be incremented in indels */
      if (op==0) {
	op = *rp++;
	if (cumm_seq_score) *i_spa = ppst->gdelval;
	d1_score +=  ppst->gdelval;
	d0_score +=  ppst->gdelval;
      }
      if (cumm_seq_score) *i_spa++ += ppst->ggapval;
      d1_score +=  ppst->ggapval; d1_alen++; d1_gaplen++;
      d0_score +=  ppst->ggapval; d0_alen++; d0_gaplen++;

      if (op > 0) {	/* insertion in aa0 */
	*sp1_p = sq[aa1p[i1]];
	*sp0_p = '-';
	*spa_p = M_DEL;

	if (calc_func_mode == CALC_CODE) {
	  update_code(align_code_dyn, update_data_p, 2, *spa_p,*sp0_p,*sp1_p);
	}

	if (have_ann) {
	  have_push_features = 0;
	  *sp0a_p = ' ';
	  if (aa1a) *sp1a_p = ann_arr[aa1a[i1]];
	  else *sp1a_p = ' ';
	  if (s_annot1_arr_p) {
	    if (i1+l_offset == s_annot1_arr_p[i1_annot]->pos || i1+l_offset == i1_left_end) {
	      i1_annot = next_annot_match(&itmp, aa0_pam2_p, l_offset+seq_pos(i1,aln->llrev,0),
					  q_offset+seq_pos(i0,aln->qlrev,0), sp1_p, sp1a_p, sq, 
					  i1_annot, annot1_p->n_annot, s_annot1_arr_p,
					  &ann_comment, annot_stack, have_push_features_p, &v_delta,
					  &d1_score, &d1_ident, &d1_alen, &d1_gaplen,
					  &left_domain_head1, left_domain_list1, &i1_left_end,
					  ppst->ggapval+ppst->gdelval);
	    }
	  }

	  if (have_ann && calc_func_mode == CALC_CODE) {
	    add_annot_code(have_ann, *sp0_p, *sp1_p, *sp0a_p, *sp1a_p,
			   q_offset + seq_pos(i0,aln->qlrev,0), l_offset+seq_pos(i1,aln->llrev,0),
			   '-', annot_var_dyn);
	  }

	  if (have_push_features && calc_func_mode != CALC_ID) {
	    display_push_features(annot_stack, annot_var_dyn,
				  q_offset+seq_pos(i0,aln->qlrev,0), *sp0_p,
				  l_offset+seq_pos(i1,aln->llrev,0), *sp1_p,
				  sim_sym[*spa_p],
				  a_res->rst.score[ppst->score_ix], a_res->rst.comp, a_res->sw_score,
				  n0, n1,
				  pstat_void, annot_fmt);
	    have_push_features = 0;
	  }
	  if (calc_func_mode == CALC_CONS) {
	    sp0a_p++;
	    sp1a_p++;
	  }
	}
	if (calc_func_mode == CALC_CONS) {
	  spa_p++;
	  sp0_p++;
	  sp1_p++;
	}
	i1++;
	op--;
	lenc++;
	aln->ngap_q++;
      }
      else {		/* (op < 0),  insertion in aa1 */
	*sp1_p = '-';
	*spa_p = M_DEL;
	*sp0_p = sq[aa0[i0]];

	if (calc_func_mode == CALC_CODE) {
	  update_code(align_code_dyn, update_data_p, 1, *spa_p,*sp0_p,*sp1_p);
	}

	if (have_ann) {
	  have_push_features = 0;
	  *sp1a_p = ' ';
	  if (aa0a) *sp0a_p = ann_arr[aa0a[i0]];
	  else *sp0a_p = ' ';
	  if (s_annot0_arr_p) {
	    if (i0+q_offset == s_annot0_arr_p[i0_annot]->pos || i0+q_offset == i0_left_end) {
	      i0_annot = next_annot_match(&itmp, ppst->pam2[0][aa1[i1]], q_offset+seq_pos(i0,aln->qlrev,0),
					  l_offset+seq_pos(i1,aln->llrev,0), sp0_p, sp0a_p, sq, 
					  i0_annot, annot0_p->n_annot, s_annot0_arr_p,
					  &ann_comment, annot_stack, have_push_features_p, &v_delta,
					  &d0_score, &d0_ident, &d0_alen, &d0_gaplen,
					  &left_domain_head0, left_domain_list0, &i0_left_end,
					  ppst->ggapval+ppst->gdelval);

	    }
	  }

	  if (calc_func_mode == CALC_CODE) {
	    add_annot_code(have_ann, *sp0_p, *sp1_p, *sp0a_p, *sp1a_p,
			   q_offset + seq_pos(i0,aln->qlrev,0), l_offset+seq_pos(i1,aln->llrev,0),
			   '-', annot_var_dyn);
	  }

	  if (have_push_features && calc_func_mode != CALC_ID) {
	    display_push_features(annot_stack, annot_var_dyn,
				  q_offset+seq_pos(i0,aln->qlrev,0), *sp0_p,
				  l_offset+seq_pos(i1,aln->llrev,0), *sp1_p,
				  sim_sym[*spa_p],
				  a_res->rst.score[ppst->score_ix], a_res->rst.comp, a_res->sw_score,
				  n0, n1, pstat_void, annot_fmt);
	    have_push_features = 0;
	  }
	  if (calc_func_mode == CALC_CONS) {
	    sp0a_p++;
	    sp1a_p++;
	  }
	}

	i0++;
	if (calc_func_mode == CALC_CONS) {
	  spa_p++;
	  sp0_p++;
	  sp1_p++;
	}
	op++;
	lenc++;
	aln->ngap_l++;
      }
    }
  }

  if (calc_func_mode == CALC_CODE) {
    close_update_data(align_code_dyn, update_data_p);
  }

  *score_delta = v_delta;

  *nc = lenc;
  if (have_ann) {
    have_push_features = 0;
    /* check for left ends after alignment */
    if (annot1_p && i1_left_end > 0) {
      close_annot_match(-1, annot_stack, have_push_features_p,
			&d1_score, &d1_ident, &d1_alen, &d1_gaplen,
			&left_domain_head1, &i1_left_end, 0);
    }

    if (annot0_p && i0_left_end > 0) {
      close_annot_match(-1, annot_stack, have_push_features_p,
			&d0_score, &d0_ident, &d0_alen, &d0_gaplen,
			&left_domain_head0, &i0_left_end, 0);
    }

    if (have_push_features && calc_func_mode != CALC_ID) {
      display_push_features(annot_stack, annot_var_dyn,
			    a_res->max0-1 + q_offset, *sp0_p,
			    a_res->max1-1 + l_offset, *sp1_p,
			    sim_sym[*spa_p],
			    a_res->rst.score[ppst->score_ix], a_res->rst.comp, a_res->sw_score,
			    n0, n1, pstat_void, annot_fmt);
    }
  }

  *spa_p = '\0';
  if (have_ann) {
    *sp0a_p = *sp1a_p = '\0';
  }
  if (calc_func_mode == CALC_CONS) {
#ifndef LCAL_CONS	/* have context around alignment */
    nd = post_fill_cons(aa0, n0, aa1p, nn1,
			a_res, ppst, mins, lenc, aln,
			seqc0, seqc1, seqc0a, seqc1a);
#else
    nd = 0;
#endif
    lenc = mins + lenc + nd;
  }

  if (have_ann) {
    if (left_domain_list0) free(left_domain_list0);
    if (left_domain_list1) free(left_domain_list1);
    annot_stack = free_stack(annot_stack);
  }

  return lenc;
}

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
  return calc_cons_u(
		     aa0, n0, aa1, n1,
		     a_res, ppst, f_str, pstat_void,
		     ann_arr, aa0a, annot0_p, aa1a, annot1_p, CALC_CONS, 0,
		     nc, seqc0, seqc1, seqca, cumm_seq_score,
		     seqc0a, seqc1a, aln, score_delta, annot_var_dyn, NULL
		     );
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
init_update_data(int show_code) {

  struct update_code_str *update_data_p;

  if ((update_data_p = (struct update_code_str *)calloc(1,sizeof(struct update_code_str)))==NULL) {
    fprintf(stderr,"*** ERROR [%s:%d] - init_update_data(): cannot allocate update_code_str\n",
	      __FILE__, __LINE__);
    return NULL;
  }

  update_data_p->p_op_idx = -1;
  update_data_p->p_op_cnt = 0;
  update_data_p->show_code = show_code;
  update_data_p->btop_enc = 0;

  if ((show_code & SHOW_CODE_CIGAR) == SHOW_CODE_CIGAR) {	/* CIGAR enc */
    update_data_p->op_map = cigar_code;
    update_data_p->cigar_order = 1;
  }
  else if ((show_code & SHOW_CODE_BTOP) == SHOW_CODE_BTOP) {	/* btop_enc */
    update_data_p->op_map = ori_code;
    update_data_p->cigar_order = 0;
    update_data_p->btop_enc = 1;
  }    
  else {   /* orig (ALIGN) enc */
    update_data_p->op_map = ori_code;
    update_data_p->cigar_order = 0;
  }

  if ((show_code & SHOW_CODE_EXT) == SHOW_CODE_EXT) {	/* set for CIGAR/ALIGN, BTOP already set */
    update_data_p->show_ext = 1;
  }
  else {
    update_data_p->show_ext = 0;
  }

  return update_data_p;
}


static void
close_update_data(struct dyn_string_str *align_code_dyn,
		  struct update_code_str *up_dp) {
  char tmp_cnt[MAX_SSTR];
  tmp_cnt[0] = '\0';

  if (!up_dp) return;

  if (up_dp->p_op_cnt) {
    if (up_dp->btop_enc) { /* btop_enc always has a p_opt_cnt == 0 unless in run of identical match */
      sprintf(tmp_cnt,"%d",up_dp->p_op_cnt);
      up_dp->p_op_cnt = 0;
    }
    else {
      sprintf_code(tmp_cnt,up_dp, up_dp->p_op_idx, up_dp->p_op_cnt);
    }
    dyn_strcat(align_code_dyn,tmp_cnt);
  }

  free(up_dp);
}

/* update_code() has been modified to work more correctly with
   ggsearch/glsearch, which, because alignments can start with either
   insertions or deletions, can produce an initial code of "0=".  When
   that happens, it is ignored and no code is added.

   *align_code_dyn - alignment string (dynamic)
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

/* only called for btop alignment encoding, for identity, update
   count, otherwise, print previous count and current difference.
   assumes that up_dp->p_op_cnt only tracks identity
*/

static void
sprintf_btop(char *tmp_str, 
	     struct update_code_str *up_dp, 
	     int op, int sim_code,
	     unsigned char sp0, unsigned char sp1)
{
  char local_str[MAX_SSTR];
  local_str[0]='\0';

  tmp_str[0] = '\0';

  if (op==0 && sim_code == M_IDENT) {
    up_dp->p_op_cnt++;
    return;
  }
  else {
    if (up_dp->p_op_cnt > 0) {
      sprintf(local_str,"%d",up_dp->p_op_cnt);
    }
    up_dp->p_op_cnt = 0;
    sprintf(tmp_str,"%s%c%c",local_str,sp0,sp1);
  }
}

static void
update_code(struct dyn_string_str *align_code_dyn, 
	    struct update_code_str *up_dp, int op, 
	    int sim_code,  unsigned char sp0, unsigned char sp1)
{
  char tmp_cnt[MAX_SSTR];
  tmp_cnt[0]='\0';

  /* op == 0 : match state (could involve termination codons);
     op == 1 : deletion
     op == 2 : insertion
     op == 3 : *:*
     p_op == 5 : mismatch state
  */

  if (up_dp->btop_enc) {
    sprintf_btop(tmp_cnt, up_dp, op, sim_code, sp0, sp1);
    dyn_strcat(align_code_dyn, tmp_cnt);
    return;
  }

  /* not btop_enc */
  if (up_dp->p_op_cnt == 0) {
    up_dp->p_op_idx = op;
    up_dp->p_op_cnt = 1;
    return;
  }

  if (op == 1 || op == 2) {
    if (up_dp->p_op_idx == op) { up_dp->p_op_cnt++;}
    else {
      sprintf_code(tmp_cnt,up_dp, up_dp->p_op_idx,up_dp->p_op_cnt);
      dyn_strcat(align_code_dyn, tmp_cnt);
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
      dyn_strcat(align_code_dyn, tmp_cnt);
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
	      struct dyn_string_str *align_code_dyn,
	      /* char *al_str, int al_str_n,  */
	      const unsigned char *ann_arr, 
	      const unsigned char *aa0a,
	      const struct annot_str *annot0_p,
	      const unsigned char *aa1a,
	      const struct annot_str *annot1_p,
	      struct dyn_string_str *annot_code_dyn,
	      int *score_delta,
	      struct f_struct *f_str,
	      void *pstat_void,
	      int display_code)
{
  int nc;

  return calc_cons_u(
		     aa0, n0, aa1, n1,
		     a_res, ppst, f_str, pstat_void,
		     ann_arr, aa0a, annot0_p, aa1a, annot1_p, CALC_CODE, 
		     display_code,
		     &nc, NULL, NULL, NULL, NULL,
		     NULL, NULL, aln, score_delta, annot_code_dyn,
		     align_code_dyn
		     );
}

/* calc_id never looks at domains or features, only variation */

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
  int nc;

  return calc_cons_u(
		     aa0, n0, aa1, n1,
		     a_res, ppst, f_str, NULL,
		     NULL, NULL, annot0_p, NULL, annot1_p, CALC_ID, 0,
		     &nc, NULL, NULL, NULL, NULL,
		     NULL, NULL, aln, score_delta, annot_var_dyn,
		     NULL
		     );
}

/* calc_id never looks at domains or features, only variation */

int calc_idd(const unsigned char *aa0, int n0,
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
  int nc;

  return calc_cons_u(
		     aa0, n0, aa1, n1,
		     a_res, ppst, f_str, NULL,
		     NULL, NULL, annot0_p, NULL, annot1_p, CALC_ID_DOM, 0,
		     &nc, NULL, NULL, NULL, NULL,
		     NULL, NULL, aln, score_delta, annot_var_dyn,
		     NULL
		     );
}
