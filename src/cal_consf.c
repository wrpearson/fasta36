/* cal_consf.c - routines for printing translated alignments for [t]fast[sf] */

/* copyright (c) 1998, 1999, 2007, 2014 by William R. Pearson and The
   Rector & Visitors of the University of Virginia */

/*  $Id: cal_consf.c 1263 2014-06-25 10:40:39Z wrp $ */

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

/* removed from dropfs2.c, dropff2.c April, 2007 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>

#include "defs.h"
#include "param.h"

#include "tatstats.h"

#include "a_mark.h"

#define DROP_INTERN
#include "drop_func.h"

void update_code(struct dyn_string_str *align_code_dyn, int op, int op_cnt, int fnum, int show_code);
extern void aancpy(char *to, char *from, int count, const struct pstruct *ppst);

int
calc_cons_a(const unsigned char *aa0, int n0,
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
  int i0, i1, nn1, n0t;
  int op, lenc, len_gap, nd, ns, itmp, p_ac, fnum, o_fnum;
  int *i_spa;
  const unsigned char *aa1p;
  const unsigned char *aa0ap;
  char *sp0, *sp0a, *sp1, *sp1a, *spa;
  int *rp;
  int mins, smins, ntmp;
  int have_ann = 0;

  *score_delta = 0;
  NULL_dyn_string(annot_var_dyn);
  have_ann = (seqc0a != NULL && (aa0a != NULL || aa1a != NULL));
  
#ifndef TFAST
  aa1p = aa1;
  nn1 = n1;
#else
  aa1p = f_str->aa1x;
  nn1 = f_str->n10;
#endif

  aln->amin0 = a_res->min0 + f_str->aa0t_off;
  aln->amax0 = a_res->max0 + f_str->aa0t_off;;
  aln->amin1 = a_res->min1;
  aln->amax1 = a_res->max1;
  aln->calc_last_set = 1;

  /* first fill in the ends */
  n0 -= (f_str->nm0-1);

  if (min(a_res->min0,a_res->min1)<aln->llen || aln->showall==1)
    			/* will we show all the start ?*/
    if (a_res->min0>=a_res->min1) {                        /* aa0 extends more to left */
      smins=0;
      if (aln->showall==1) mins=a_res->min0;
      else mins = min(a_res->min0,aln->llen/2);
      aancpy(seqc0,(char *)f_str->aa0t+a_res->min0-mins,mins,ppst);
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
      else mins = min(a_res->min1,aln->llen/2);
      aancpy(seqc1,(char *)(aa1p+a_res->min1-mins),mins,ppst);
      aln->smin1 = a_res->min1-mins;
      if ((mins-a_res->min0)>0) {
	memset(seqc0,' ',mins-a_res->min0);
	aancpy(seqc0+mins-a_res->min0,(char *)f_str->aa0t,a_res->min0,ppst);
	aln->smin0 = 0;
      }
      else {
	aancpy(seqc0,(char *)f_str->aa0t+a_res->min0-mins,mins,ppst);
	aln->smin0 = a_res->min0-mins;
      }
    }
  else {
    mins= min(aln->llen/2,min(a_res->min0,a_res->min1));
    smins=mins;
    aln->smin0=a_res->min0;
    aln->smin1=a_res->min1;
    aancpy(seqc0,(char *)f_str->aa0t+a_res->min0-mins,mins,ppst);
    aancpy(seqc1,(char *)aa1p+a_res->min1-mins,mins,ppst);
  }

  memset(seqca,M_BLANK,mins);
  if (have_ann) {
    /* pad annotation before alignment - this strategy means no
       annotation before alignment */
    memset(seqc0a,' ', mins);
    memset(seqc1a,' ', mins);
  }

/* now get the middle */

  spa = seqca+mins;
  if (cumm_seq_score) i_spa = cumm_seq_score+mins;
  sp0 = seqc0+mins;
  sp0a = seqc0a+mins;
  sp1 = seqc1+mins;
  sp1a = seqc1a+mins;
  rp = a_res->res;
  n0t=lenc=len_gap=aln->nident=aln->nmismatch=aln->nsim=aln->npos=aln->ngap_q=aln->ngap_l=op=p_ac= 0;
  i0 = a_res->min0;
  i1 = a_res->min1;
  
  /* op is the previous "match/insert" operator; *rp is the current
     operator or repeat count */

#if defined(FASTS) || defined(FASTM)
  o_fnum = f_str->aa0ti[i0];
  if (aa0a) aa0ap = &aa0a[f_str->nmoff[o_fnum]+i0];
#endif

  while (i0 < a_res->max0 || i1 < a_res->max1) {
#if defined(FASTS) || defined(FASTM)
    fnum = f_str->aa0ti[i0];
#endif

    if (op == 0 && *rp == 0) {	/* previous was match (or start), current is match */

#if defined(FASTS) || defined(FASTM)
      if (p_ac == 0) { /* previous code was a match */
	if (fnum != o_fnum) { /* continuing a match, but with a different fragment */
	  if (have_ann) { aa0ap = &aa0a[f_str->nmoff[fnum]];}
	  o_fnum = fnum;
	}
      }
      else {
	p_ac = 0; o_fnum = fnum = f_str->aa0ti[i0];
	if (have_ann) {aa0ap = &aa0a[f_str->nmoff[fnum]];}
      }
#endif
      op = *rp++;		/* get the next match/insert operator */

      /* get the alignment symbol */
      if ((itmp=ppst->pam2[0][f_str->aa0t[i0]][aa1p[i1]])<0) { *spa = M_NEG; }
      else if (itmp == 0) { *spa = M_ZERO;}
      else {*spa = M_POS;}
      if (*spa == M_POS)  { aln->npos++;}
      if (*spa == M_ZERO || *spa == M_POS) { aln->nsim++;}

      if (cumm_seq_score) *i_spa++ += itmp;

      *sp0 = ppst->sq[f_str->aa0t[i0++]];	/* get the residues for the consensus */

      if (have_ann) {
	if (aa0a) {*sp0a++ = ann_arr[*aa0ap++];}
	else {*sp0a++ = ' ';}
	if (aa1a) {*sp1a++ = ann_arr[aa1a[i1]];}
	else {*sp1a++ = ' ';}
      }

      *sp1 = ppst->sq[aa1p[i1++]];
      n0t++;
      lenc++;
      if (toupper(*sp0) == toupper(*sp1)) {aln->nident++; *spa = M_IDENT;}
      else {aln->nmismatch++;}
      sp0++; sp1++; spa++;
    }
    else {	/* either op != 0 (previous was insert) or *rp != 0
		   (current is insert) */
      if (op==0) { op = *rp++;}	/* previous was match, start insert */
      				/* previous was insert - count through gap */
#if defined(FASTS) || defined(FASTM)
      if (p_ac != 1) {
	p_ac = 1;
	fnum = f_str->aa0ti[i0];
      }
#endif
      if (have_ann) {
	*sp0a++ = ' ';
	if (aa1a) {*sp1a++ = ann_arr[aa1a[i1]];}
	else {*sp1a++ = ' ';}
      }
      *sp0++ = '-';
      *sp1++ = ppst->sq[aa1p[i1++]];
      *spa++ = M_DEL;
      op--;
      len_gap++;
      lenc++;
    }
  }	/* end alignment while() */

  if (have_ann) {*sp0a = *sp1a = '\0';}
  *spa = '\0';
  *nc = lenc-len_gap;

  /* now we have the middle, get the right end */

  /* ns should be the length of alignmnet display in seqc0/0a/1/1a */
  ns = mins + lenc + aln->llen;
  /* adjust for the last line */
  ns -= (itmp = ns %aln->llen);
  /* add another full line */
  if (itmp>aln->llen/2) ns += aln->llen;
  /* the amount to display at the end, now (after lenc) */
  nd = ns - (mins+lenc);
  if (nd > max(n0t-a_res->max0,nn1-a_res->max1)) nd = max(n0t-a_res->max0,nn1-a_res->max1);
  
  if (aln->showall==1) {
    nd = max(n0t-a_res->max0,nn1-a_res->max1);	/* reset for showall=1 */
    /* get right end */
    aancpy(seqc0+mins+lenc,(char *)f_str->aa0t+a_res->max0,n0t-a_res->max0,ppst);
    aancpy(seqc1+mins+lenc,(char *)aa1p+a_res->max1,nn1-a_res->max1,ppst);
    /* fill with blanks - this is required to use one 'nc' */
    memset(seqc0+mins+lenc+n0t-a_res->max0,' ',nd-(n0t-a_res->max0));
    memset(seqc1+mins+lenc+nn1-a_res->max1,' ',nd-(nn1-a_res->max1));
  }
  else {
    if ((nd-(n0t-a_res->max0))>0) {
      /* finish copying out the sequence */
      aancpy(seqc0+mins+lenc,(char *)f_str->aa0t+a_res->max0,
	     n0t-a_res->max0,ppst);
      /* add blanks to pad */
      memset(seqc0+mins+lenc+n0t-a_res->max0,' ',nd-(n0t-a_res->max0));
    }
    else {
      /* just use up some sequence */
      aancpy(seqc0+mins+lenc,(char *)f_str->aa0t+a_res->max0,nd,ppst);
    }
    if ((nd-(nn1-a_res->max1))>0) {
      /* finish copying out the sequence */
      aancpy(seqc1+mins+lenc,(char *)aa1p+a_res->max1,nn1-a_res->max1,ppst);
      /* add blanks to pad */
      memset(seqc1+mins+lenc+nn1-a_res->max1,' ',nd-(nn1-a_res->max1));
    }
    else {
      /* just use up some sequence */
      aancpy(seqc1+mins+lenc,(char *)aa1p+a_res->max1,nd,ppst);
    }
  }
  if (have_ann) {
    /* also pad the annotation -- this strategy means no annotations
       in the unaligned region*/
    memset(seqc0a+mins+lenc,' ',nd);
    memset(seqc1a+mins+lenc,' ',nd);
    /* 
    ntmp = nd-(n0t-a_res->max0);
    if (ntmp > 0) memset(seqc0a+mins+lenc+n0-a_res->max0,' ',ntmp);
    ntmp = nd-(nn1-a_res->max1);
    if (ntmp > 0) memset(seqc1a+mins+lenc+nn1-a_res->max1,' ',ntmp);
    */
  }

  aln->smin0 = f_str->aa0t_off;

  return mins+lenc+nd;
}

void
calc_astruct(struct a_struct *aln_p, struct a_res_str *a_res_p, struct f_struct *f_str) {

  /* we do not pay attention to aln_p->calc_last_set, because all the
     functions (calc_astruct, calc_cons_a, calc_code) use exactly the same
     assignment */

  aln_p->amin0 = a_res_p->min0 + f_str->aa0t_off;
  aln_p->amax0 = a_res_p->max0 + f_str->aa0t_off;
  aln_p->amin1 = a_res_p->min1;
  aln_p->amax1 = a_res_p->max1;
}

/* build an array of match/ins/del - length strings */
int
calc_code(const unsigned char *aa0, const int n0,
	  const unsigned char *aa1, const int n1,
	  struct a_struct *aln,
	  struct a_res_str *a_res,
	  struct pstruct *ppst,
	  struct dyn_string_str *align_code_dyn,
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
  int i0, i1, nn1;
  int op, lenc, len_gap;
  int p_ac, op_cnt;
  const unsigned char *aa1p;
  const unsigned char *aa0ap;
  char tmp_cnt[20];
  char sp0, sp1, spa;
  unsigned char *sq;
  int *rp;
  int mins, smins;
  int o_fnum,fnum = 0;

  int have_ann = 0;
  char ann_ch0, ann_ch1;
  char tmp_astr[MAX_STR];
  int sim_code;
  int show_code, annot_fmt;
  char *sim_sym= aln_map_sym[MX_ACC];

  *score_delta = 0;

  show_code = (display_code & SHOW_CODE_MASK);
  annot_fmt = 2;
  if (display_code & SHOW_ANNOT_FULL) {
    annot_fmt = 1;
  }

  if (aa0a != NULL && aa1a != NULL) { have_ann = 2;}
  else if (aa0a != NULL || aa1a != NULL) { have_ann = 1;}
  else {have_ann = 0;}

  if (ppst->ext_sq_set) {sq = ppst->sqx;}
  else {sq = ppst->sq;}

#ifndef TFAST
  aa1p = aa1;
  nn1 = n1;
#else
  aa1p = f_str->aa1x;
  nn1 = f_str->n10;
#endif

  /* 
  aln->amin0 = a_res->min0;
  aln->amin1 = a_res->min1;
  aln->amax0 = a_res->max0;
  aln->amax1 = a_res->max1;
  */

  rp = a_res->res;
  lenc = len_gap =aln->nident=aln->nmismatch=aln->nsim=aln->npos=aln->ngap_q=aln->ngap_l=aln->nfs=op=p_ac = 0;
  op_cnt = 0;

  i0 = a_res->min0;	/* start in aa0 (f_str->aa0t) */
  i1 = a_res->min1;	/* start in aa1 */
  tmp_cnt[0]='\0';
  
#if defined(FASTS) || defined(FASTM)
  o_fnum = f_str->aa0ti[i0]+1;
  if (aa0a) { aa0ap = &aa0a[f_str->nmoff[o_fnum]+i0]; }
#endif

  while (i0 < a_res->max0 || i1 < a_res->max1) {
    fnum = f_str->aa0ti[i0]+1;
    if (op == 0 && *rp == 0) {	/* previous was match, this is match */
#if defined(FASTS) || defined(FASTM)
      if (p_ac == 0) {	/* previous code was a match */
	if (fnum == o_fnum) { op_cnt++;	}
	else {		/* continuing a match, but with a different fragment */
	  update_code(align_code_dyn, p_ac, op_cnt, o_fnum, show_code);
	  if (have_ann) aa0ap = &aa0a[f_str->nmoff[fnum]];
	  o_fnum = fnum;
	  op_cnt=1;
	}
      }
      else {
	update_code(align_code_dyn,p_ac,op_cnt,o_fnum, show_code);
	op_cnt = 1; p_ac = 0; o_fnum = fnum = f_str->aa0ti[i0] + 1;
	if (have_ann) {aa0ap = &aa0a[f_str->nmoff[fnum]];}
      }
#endif
      op = *rp++;
      lenc++;
      sim_code = M_NEG;
      if (ppst->pam2[0][f_str->aa0t[i0]][aa1p[i1]]>0) {
	sim_code = M_POS;
	aln->npos++;
	aln->nsim++;
      }
      else if (ppst->pam2[0][f_str->aa0t[i0]][aa1p[i1]]==0) {
	sim_code = M_ZERO;
	aln->nsim++;
      }

      sp0 = ppst->sq[f_str->aa0t[i0]];
      sp1 = ppst->sq[aa1p[i1]];
      if (toupper(sp0) == toupper(sp1)) {
	sim_code = M_IDENT;
	aln->nident++;
      }
      else {aln->nmismatch++;}

      /* check for an annotation */
      if (have_ann) {
	ann_ch0 = ann_ch1 = '\0';
	if (have_ann == 2 && (ann_arr[*aa0ap] != ' ' || ann_arr[aa1a[i1]] != ' ')) {
	  ann_ch0 = ann_arr[*aa0ap];
	  if (ann_ch0 == ' ') ann_ch0 = 'X';
	  ann_ch1 = ann_arr[aa1a[i1]];
	  if (ann_ch1 == ' ') ann_ch1 = 'X';
	}
	else if (aa0a != NULL && ann_arr[*aa0ap]!=' ') {
	  ann_ch0 = ann_arr[*aa0ap];
	  ann_ch1 = 'X';
	}
	else if (aa1a != NULL && ann_arr[aa1a[i1]]!=' ') {
	  ann_ch0 = 'X';
	  ann_ch1 = ann_arr[aa1a[i1]];
	}
	aa0ap++;
	if (ann_ch0) {
	  sprintf(tmp_astr, "|%ld:%ld:%c%c:%c%c%c",
		  aln->q_offset+i0+1,aln->l_offset+i1+1,
		  ann_ch0,ann_ch1,sim_sym[sim_code],sp0,sp1);
	  /* strncat(ann_str, tmp_astr, ann_str_n - strlen(ann_str) - 1); */
	  dyn_strcat(annot_code_dyn, tmp_astr);
	}
      }
      i0++;
      i1++;
    }
    else {
      if (op==0) op = *rp++;
      if (p_ac == 1) { op_cnt++;}
      else {
	update_code(align_code_dyn,p_ac,op_cnt,o_fnum, show_code);
#if defined(FASTS) || defined(FASTM)
	p_ac = 1;
	fnum = f_str->aa0ti[i0];
#endif
	op_cnt = 1; fnum = f_str->aa0ti[i0] + 1;
      }
      op--; lenc++; i1++; len_gap++;
    }
  }
  update_code(align_code_dyn,p_ac,op_cnt,o_fnum, show_code);

  return lenc - len_gap;
}

/* update_code(): if "op" == 0, this is the end of a match of length
   "op_cnt" involving fragment "fnum"
   otherwise, this is an insertion (op==1) or deletion (op==2)
*/

void
update_code(struct dyn_string_str *align_code_dyn, int op, int op_cnt, int fnum, int show_code) {

  char align_char[4]={"=-+"};
  char cigar_char[4]={"MDI"};
  char tmp_cnt[20];

  if (op_cnt == 0) return;

  if (show_code == SHOW_CODE_CIGAR) {
    sprintf(tmp_cnt,"%d%c",op_cnt,cigar_char[op]);
  }
  else {
    if (op == 0)
      sprintf(tmp_cnt,"%c%d[%d]",align_char[op],op_cnt,fnum);
    else
      sprintf(tmp_cnt,"%c%d",align_char[op],op_cnt);
  }
  dyn_strcat(align_code_dyn, tmp_cnt);
}

int
calc_id(const unsigned char *aa0, int n0,
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
  int op, lenc, len_gap;
  const unsigned char *aa1p;
  int sp0, sp1;
  int *rp;
  int mins, smins;
  
  *score_delta = 0;
  NULL_dyn_string(annot_var_dyn);

#ifndef TFAST
  aa1p = aa1;
  nn1 = n1;
#else
  aa1p = f_str->aa1x;
  nn1 = f_str->n10;
#endif

  aln->amin0 = a_res->min0 + f_str->aa0t_off;
  aln->amax0 = a_res->max0 + f_str->aa0t_off;
  aln->amin1 = a_res->min1;
  aln->amax1 = a_res->max1;
  aln->calc_last_set = 1;

  /* first fill in the ends */
  n0 -= (f_str->nm0-1);

  /* now get the middle */
  rp = a_res->res;
  lenc=len_gap=aln->nident=aln->nmismatch=aln->nsim=aln->npos=aln->ngap_q = aln->ngap_l = aln->nfs = op = 0;
  i0 = a_res->min0;
  i1 = a_res->min1;
  
  while (i0 < a_res->max0 || i1 < a_res->max1) {
    if (op == 0 && *rp == 0) {
      op = *rp++;

      if (ppst->pam2[0][f_str->aa0t[i0]][aa1p[i1]]>0) {
	aln->nsim++;
	aln->npos++;
      }
      else if (ppst->pam2[0][f_str->aa0t[i0]][aa1p[i1]]==0) {
	aln->nsim++;
      }

      sp0 = ppst->sq[f_str->aa0t[i0++]];
      sp1 = ppst->sq[aa1p[i1++]];
      lenc++;
      if (toupper(sp0) == toupper(sp1)) aln->nident++;
      else {aln->nmismatch++;}
    }
    else {
      if (op==0) { op = *rp++;}
      i1++;
      op--;
      len_gap++;
      lenc++;
    }
  }
  return lenc-len_gap;
}

int
calc_idd(const unsigned char *aa0, int n0,
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
  return calc_id(aa0,n0,aa1,n1,aln, a_res, ppst, annot0_p, annot1_p, score_delta, annot_var_dyn, f_str);
}
