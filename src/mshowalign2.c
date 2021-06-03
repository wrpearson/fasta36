/* $Id: mshowalign2.c 1269 2014-07-29 21:24:25Z wrp $ */

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

/* mshowalign.c - show sequence alignments in pvcomplib */

/* 
   In the serial and current threaded versions of the programs,
   showalign gets a list of high scoring sequences and must
   re_getlib() the sequence, do_walign(), and then calculate the
   alignment.

*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "defs.h"
#include "msg.h"
#include "structs.h"
#include "param.h"

#include "mm_file.h"

/* best_stats.h must come after mm_file.h */
#include "best_stats.h"

/* used to position the library sequence for re_getlib - also gets
   description */
#define RANLIB (m_fptr->ranlib)

extern struct lmf_str *
re_openlib(struct lmf_str *, int outtty);

int
re_getlib(unsigned char *aa1, struct annot_str **annot_p,
	  int maxn, int maxt,
	  int loff, int cont, int term_code,
	  long *l_offset, long *l_off, 
	  struct lmf_str *m_fptr);

#include "drop_func.h"
/* drop_func.c includes dyn_string.h */

extern void calc_astruct(struct a_struct *aln_p, struct a_res_str *a_res_p, void *f_str);

extern void calc_coord(int n0, int n1, long qoffset, long loffset,
		       struct a_struct *aln);

void initseq(char **, char **, char **, int);
void initseq_ann(char **, char **, int);

void freeseq(char **, char **, char **);
void freeseq_ann(char **, char **);

void do_show(FILE *fp, int n0, int n1, int score,
	     char *name0, char *name1, int nml, char *link_name,
	     const struct mngmsg *m_msp, const struct pstruct *ppst,
	     char *seqc0, char *seqc0a,  char *seqc1, char *seqc1a,
	     char *seqca, int *cumm_seq_score, int nc,
	     float percent, float gpercent, int lc,
	     struct a_struct *aln, const char *annot_var_s,
	     const struct annot_str *q_annot_p,
	     const struct annot_str *l_annot_p);

void
do_lav(FILE *fp, struct a_struct *aln, char *seqc, float percent, int is_mirror);

void
buf_align_seq(unsigned char **aa0, int n0,
	      struct beststr **bestp_arr, int nbest,
	      struct pstruct *ppst, struct mngmsg *m_msp,
	      const struct mng_thr *m_bufi_p
#if !defined(COMP_THR) && !defined(PCOMPLIB)
	      , void **f_str
#endif
	      );

/* pre-alignment */
extern void 
pre_load_best(unsigned char *aa1, int maxn,struct beststr **bbp_arr,
	      int nbest, struct mngmsg *m_msp, int debug);

float
calc_fpercent_id(float scale, int n_ident, int n_alen, int tot_ident, float fail);

extern int E1_to_s(double e_val, int n0, int n1, int db_size, void *pu);

extern void discons(FILE *fd, const struct mngmsg *m_msg,
		    char *seqc0, char *seqc0a,
		    char *seqc1, char *seqc1a,
		    char *seqca, int *cumm_seq_score, int nc, 
		    int n0, int n1, char *name0, char *name1, int nml,
		    struct a_struct *aln);

extern void disgraph(FILE *fd, int n0, int n1,
		     float percent, int score,
		     int min0, int min1, int max0, int max1, long sq0off,
		     char *name0, char *name1, int nml, int llen, int markx);

extern double find_z(int score, double escore, int length, double comp,void *);
extern double zs_to_bit(double, int, int);
extern double s_to_bit(int score, int n0, int  n1, void *pu);
extern double bit_to_E (double bit, int n0, int n1, long db_size, void *pu);
extern double zs_to_E(double zs, int n1, int dnaseq, long db_size, struct db_str db);

extern void
do_url1(FILE *, const struct mngmsg *, const struct pstruct *, char *, int,
	const struct a_struct *, const char *,
	const struct annot_str *, const struct annot_str *);

#ifndef A_MARK
#define A_MARK ">>"
#endif

/* this version does not check for m_msg->e_cut because nshow/nbest has
   already been set to limit on e_cut */

void showalign (FILE *fp, unsigned char **aa0, unsigned char *aa1save, int maxn,
		struct beststr **bptr, int nbest, int qlib, 
		struct mngmsg *m_msp, struct pstruct *ppst, 
		char *info_gstring2
		, void **f_str, struct mng_thr *m_bufi_p
		)
{
  unsigned char *aa1, *aa1a;
  char tmp_str[20];
  char info_str[200];
  char bline[2048], *qline_p, *bline_p, *bl_ptr, *bp, fmt[40];
#ifdef LALIGN
  char *bp1;
#endif
  struct dyn_string_str *annot_var_dyn, *align_code_dyn;
  char *annot_var_s10;
  int tmp_len, ttmp_len, l_llen, desc_llen, ranlib_done;
  char name0[80], name0s[80], name1[200];
  char l_name[128], link_name[140];	/* link name */
  int istop, i = 0, ib, nml, first_line;
  int l_ashow;
  int  n1tot;
  struct beststr *bbp;
  struct a_res_str *cur_ares_p;
  struct rstruct *rst_p;
  int nc, lc, maxc;
  double lzscore, lzscore2, lbits;
  struct a_struct l_aln, *l_aln_p;
  float percent, gpercent, ng_percent, disp_percent, disp_similar;
  int disp_alen;
  /* strings, lengths for conventional alignment */
  char *seqc0, *seqc0a, *seqc1, *seqc1a, *seqca;
  int *cumm_seq_score;
  /* strings, lengths, for encoded alignment for MX10 */
  char *seq_code=NULL, *annot_code=NULL;
  int seq_code_len=0, annot_code_len=0;
  long loffset, l_off;
  long qt_offset, lt_offset;
  int lsw_score, l_score0;
  char html_pre_E[120], html_post_E[120];
  int disp_dna_align = ((m_msp->qdnaseq>0) && (m_msp->ldb_info.ldnaseq > 0));
  int score_delta = 0;
#ifdef LALIGN
  int lalign_repeat_thresh_done = 0;
#endif

  int n1;
  struct lmf_str *m_fptr;
  int ngap;

  align_code_dyn = init_dyn_string(4096, 4096);
  annot_var_dyn = init_dyn_string(4096, 4096);

  qline_p = m_msp->qtitle;
  if (!m_msp->gi_save && !strncmp(m_msp->qtitle,"gi|",3)) {
    qline_p = strchr(qline_p+4,'|');
    /* check for additional '|'s associated with NCBI gi|12346|db|acc entry */
    if (!qline_p || strchr(qline_p+1,'|')==NULL) {
      qline_p = m_msp->qtitle;
    }
    else { qline_p += 1;}
  }

  memcpy(&l_aln, &(m_msp->aln),sizeof(struct a_struct));
  l_aln_p = &l_aln;  /*   aln_p = &m_msp->aln; */

  /* set the name0,1 label length */
  if (m_msp->markx & (MX_M10FORM+MX_MBLAST)) nml = 12;
  else if (m_msp->markx & MX_M11OUT) nml = MAX_UID;
  else nml = m_msp->nmlen;

  if (strlen(qline_p) > 0) {
    if (qline_p[0]=='>') {SAFE_STRNCPY(name0s,qline_p+1,sizeof(name0s));}
    else {SAFE_STRNCPY(name0s,qline_p,sizeof(name0s));}
  }
  else {
    SAFE_STRNCPY(name0s,m_msp->tname,sizeof(name0s));
  }

  if ((bp=strchr(name0s,' '))!=NULL) *bp='\0';

  if (m_msp->revcomp) name0[nml-1]='-';

  if (m_msp->markx & MX_HTML) {
    SAFE_STRNCPY(html_pre_E,"<font color=\"darkred\">",sizeof(html_pre_E));
    SAFE_STRNCPY(html_post_E,"</font>",sizeof(html_post_E));

  }
  else {
    html_pre_E[0] = html_post_E[0] = '\0';
  }

  desc_llen = l_llen = m_msp->aln.llen;
  if ((m_msp->markx & MX_M9SUMM) && (m_msp->show_code != SHOW_CODE_ID && m_msp->show_code != SHOW_CODE_IDD)) {
    l_llen += 40;
    if (l_llen > 200) l_llen=200;
  }

  if (m_msp->markx & MX_MBLAST) {
    sprintf(fmt,">%%-%ds\n%%sLength=%%d\n",l_llen+15);
    desc_llen = l_llen+25;
  }
  else {
    sprintf(fmt,"%s%%-%ds (%%d %s)\n",A_MARK,l_llen-5,m_msp->sqnam);
  }

  if (m_msp->std_output && !(m_msp->markx&MX_M10FORM)) fprintf (fp,"\n");

  l_ashow = m_msp->ashow;
  if (l_ashow < 0) l_ashow = m_msp->nshow;
  istop = min(min(nbest,l_ashow),m_msp->nshow);

  tmp_len = sizeof(bline)-1;
  if (!(m_msp->markx & MX_M10FORM) && !m_msp->long_info) {tmp_len = l_llen-5;}
  
  /* don't call pre_load_best if we already have sequences and alignments */
  if (!m_msp->align_done) {
    if (!m_msp->pre_load_done) { pre_load_best(aa1save, maxn, bptr, istop, m_msp, ppst->debug_lib); }

    /* don't call buf_align_seq if the algorithm does not support
       pre-alignment */
    if (ppst->can_pre_align) {
#ifdef LALIGN
      for (ib=0; ib<istop; ib++) { 
	bbp = bptr[ib];
	bbp->repeat_thresh = 
	  min(E1_to_s(ppst->e_cut_r, m_msp->n0, bbp->seq->n1,ppst->zdb_size, m_msp->pstat_void),
	      bbp->rst.score[ppst->score_ix]);
      }
      lalign_repeat_thresh_done = 1;
#endif

      buf_align_seq(aa0, m_msp->n0, bptr, istop, ppst, m_msp, m_bufi_p
#if !defined(COMP_THR) && !defined(PCOMPLIB)
		    , f_str
#endif
		    );
    }
  }

  for (ib=0; ib<istop; ib++) {
    bbp = bptr[ib];

#ifdef LALIGN
    if (!lalign_repeat_thresh_done) {
      bbp->repeat_thresh = 
	min(E1_to_s(ppst->e_cut_r, m_msp->n0, bbp->seq->n1,ppst->zdb_size, m_msp->pstat_void),
	    bbp->rst.score[ppst->score_ix]);
    }
#endif
    /* preload stuff guaranteed to be in bbp->seq */
    n1 = bbp->seq->n1;
    aa1 = bbp->seq->aa1b;
    if (bbp->seq->annot_p != NULL) aa1a = bbp->seq->annot_p->aa1_ann;
    else aa1a = NULL;
    l_off = bbp->seq->l_off;
    loffset = bbp->seq->l_offset;

    /* make sure we have a description */
    if (bbp->mseq->bline == NULL || bbp->mseq->bline_max < tmp_len) {
      if ((m_fptr=re_openlib(bbp->mseq->m_file_p,!m_msp->quiet))==NULL)
	exit(1);
      RANLIB(bline,tmp_len,bbp->mseq->lseek,bbp->mseq->libstr,bbp->mseq->m_file_p);
      bline[tmp_len]='\0';
      ranlib_done = 1;
    }
    else {
      ranlib_done = 0;
      SAFE_STRNCPY(bline, bbp->mseq->bline, sizeof(bline));
    }

    /* make sure we have a sequence */
    if (bbp->seq->aa1b == NULL || (m_msp->ann_flg==1 && &(bbp->seq->annot_p)==NULL)) {
      if (!ranlib_done) {
	if ((m_fptr=re_openlib(bbp->mseq->m_file_p,!m_msp->quiet))==NULL)
	  exit(1);
	RANLIB(bline,tmp_len,bbp->mseq->lseek,bbp->mseq->libstr,bbp->mseq->m_file_p);
	bline[tmp_len]='\0';
	ranlib_done = 1;
      }

      n1 = re_getlib(aa1save, (m_msp->ann_flg==1) ? &(bbp->seq->annot_p) : NULL, maxn,
		     m_msp->ldb_info.maxt3,
		     m_msp->ldb_info.l_overlap,bbp->mseq->cont,m_msp->ldb_info.term_code,
		     &loffset,&l_off,bbp->mseq->m_file_p);
      aa1 = aa1save;
      if (m_msp->ann_flg==1 && bbp->seq->annot_p->aa1_ann) {aa1a = bbp->seq->annot_p->aa1_ann;}
    }

#ifdef DEBUG
    if (n1 != bbp->seq->n1) {
      fprintf(stderr," library sequence: %s lengths differ: %d != %d\n",
	      bline,bbp->seq->n1, n1);
      fprintf(stderr, "offset is: %lld\n",bbp->mseq->lseek);
    }
#endif

    /* make sure we have an alignment encoding */
    if (!(bbp->have_ares & 0x1)) {

      bbp->a_res = do_walign(aa0[bbp->frame],m_msp->n0, aa1, n1,
			     bbp->frame, bbp->repeat_thresh, ppst,
			     f_str[bbp->frame], &bbp->have_ares);
    }
    else {
      pre_cons(aa1,n1,bbp->frame,f_str[bbp->frame]);
    }

    cur_ares_p = bbp->a_res;
    /* current do_walign()'s provide valid rst */
    /*
    memcpy(&cur_ares_p->rst,&bbp->rst, sizeof(struct rstruct));
    */
    aln_func_vals(bbp->frame, l_aln_p);

    if (strlen(bline)==0) {
      SAFE_STRNCPY(bline,">",sizeof(bline));
      SAFE_STRNCAT(bline,m_msp->lname,l_llen-5);
    }

    bline_p = bline;
    /* always remove "gi|" for alignments */
    if (!m_msp->gi_save && !strncmp(bline,"gi|",3)) {
      bline_p = strchr(bline+4,'|');
      if (!bline_p || !strchr(bline_p+1,'|')) {bline_p = bline;}
      else bline_p += 1;
    }

    /* re-format bline */
    while ((bp=strchr(bline_p,'\n'))!=NULL) *bp=' ';
    if (m_msp->long_info) {
      ttmp_len = strlen(bline_p);
      bl_ptr = bline_p;
      if (!(m_msp->markx & MX_M10FORM)) {
	while (ttmp_len > desc_llen) {
	  for (i=desc_llen; i>10; i--)
	    if (bl_ptr[i]==' ') {
	      bl_ptr[i]='\n';
	      break;
	    }
	  if (i <= 10) break;
	  ttmp_len -= i;
	  bl_ptr += i;
	}
      }
      bline[tmp_len]='\0';
    }

    n1tot = (bbp->mseq->n1tot_p) ? *bbp->mseq->n1tot_p : bbp->seq->n1;

    /* name1 is used to label the display */
    /* bline_p does not have gi|12345, but could have pf26|12345 or sp|P09488 */
    SAFE_STRNCPY(name1,bline_p,sizeof(name1));

    if (!(m_msp->markx & MX_M10FORM)) name1[nml]='\0';
    if ((bp = strchr(name1,' '))!=NULL) *bp = '\0';

    /* l_name is used to build an HTML link from the bestscore line to
       the alignment.  It can also be used to discriminate multiple hits
       from the same long sequence.  Text must match that in mshowbest.c */

    SAFE_STRNCPY(l_name,bline_p,sizeof(l_name));
    l_name[sizeof(l_name)-1]='\0';
    if ((bp = strchr(l_name,' '))!=NULL) *bp = '\0';
    if ((bp=strchr(&l_name[6],'|'))!=NULL) *bp='\0';  	/* increase to [6] from [3] to allow longer db names "ref", "unk", */
    if (m_msp->nframe > 2) sprintf(&l_name[strlen(l_name)],"_%d",bbp->frame+1);
    else if (m_msp->qframe >= 0 && bbp->frame == 1) {
      SAFE_STRNCAT(l_name,"_r",sizeof(l_name));
    }
    if (bbp->mseq->cont-1 > 0) {
      sprintf(tmp_str,":%d",bbp->mseq->cont-1);
      SAFE_STRNCAT(l_name,tmp_str,sizeof(l_name)-strlen(l_name));
    }

    if (m_msp->markx & MX_MBLAST) { SAFE_STRNCPY(name1,"Sbjct",sizeof(name1));}

    if (!(m_msp->markx & MX_M10FORM)) name1[nml]='\0';

    /* print out score information; */

    if (m_msp->markx & MX_HTML ) {
      SAFE_STRNCPY(link_name, l_name, sizeof(link_name));
      fprintf (fp,"<a name=\"%s\"><pre>",link_name);
    }
    SAFE_STRNCPY(name0,name0s,nml+1);
    if (m_msp->markx & MX_MBLAST) { SAFE_STRNCPY(name0,"Query",sizeof(name0));}
    name0[nml]='\0';

    if (ppst->zsflag%10 == 6) {
      sprintf(info_str," comp: %.5f H: %.5f",bbp->rst.comp,bbp->rst.H);
    }
    else info_str[0]='\0';

    if (m_msp->markx & MX_M11OUT) {
      qt_offset = m_msp->q_offset + (m_msp->q_off-1)+(m_msp->sq0off);
      lt_offset = loffset + (bbp->seq->l_off-1) + (m_msp->sq1off);
      fprintf (fp, "s {\n   \"%s\" %ld %ld \n   \"%s\" %ld %ld\n}\n",
	       name0, qt_offset, qt_offset + m_msp->n0 - 1,
	       name1, lt_offset, lt_offset + bbp->seq->n1 - 1);
      fprintf (fp, "h {\n   \"%s\"\n   \"%s\"\n}\n", qline_p, bline_p);
    }


    /* enables >>seq_acc seq_description length for first alignment, >- after */
    first_line = 1;

    while (cur_ares_p != NULL && cur_ares_p->nres > 0) {

      /* estimate space for alignment consensus */
      if (m_msp->aln.showall==1) {
	maxc = cur_ares_p->nres + max(cur_ares_p->min0,cur_ares_p->min1)+
	  max((m_msp->n0-cur_ares_p->max0),(n1-cur_ares_p->max1))+4;
      }
      else {
	maxc = cur_ares_p->nres + 4*m_msp->aln.llen+4;
      }

      /* get space to put the sequence alignment consensus */
      initseq(&seqc0, &seqc1, &seqca, maxc);
      cumm_seq_score = NULL;
      if (m_msp->markx & MX_RES_ALIGN_SCORE) {
	if ((cumm_seq_score = (int *)calloc(maxc,sizeof(int)))==NULL) {
	  fprintf(stderr,"***error*** [%s:%d] cannot allocate cumm_seq_score[%d]\n",
		  __FILE__, __LINE__, (int)(maxc*sizeof(int)));
	}
      }
      if (m_msp->ann_flg && (m_msp->aa0a != NULL || aa1a!=NULL || m_msp->annot_p)) {
	initseq_ann(&seqc0a, &seqc1a, maxc);
      }
      else { seqc0a = seqc1a = NULL;}

      calc_astruct(l_aln_p, cur_ares_p, f_str[bbp->frame]);

      calc_coord(m_msp->n0,bbp->seq->n1,
		 m_msp->q_offset+(m_msp->q_off-1)+(m_msp->sq0off-1),
		 loffset+(l_off-1)+(m_msp->sq1off-1),
		 l_aln_p);

#ifdef LALIGN
      if ((m_msp->markx & MX_M11OUT) == MX_M11OUT) {	/* lav output - skip lots of stuff */
	lsw_score = cur_ares_p->sw_score;
	lzscore = find_z(lsw_score, 0.0, bbp->seq->n1, 0.0, m_msp->pstat_void);
	lzscore2 = find_z(lsw_score, 0.0, bbp->seq->n1, 0.0, m_msp->pstat_void2);
	lbits = zs_to_bit(lzscore, m_msp->n0, bbp->seq->n1);

	NULL_dyn_string(annot_var_dyn);
	NULL_dyn_string(align_code_dyn);

	lc=calc_code(aa0[bbp->frame],m_msp->n0,
		     aa1,n1, 
		     l_aln_p, cur_ares_p,
		     ppst, 
		     align_code_dyn,
		     /*	seqc0, maxc, */
		     m_msp->ann_arr,
		     m_msp->aa0a, m_msp->annot_p,
		     aa1a, bbp->seq->annot_p,
		     annot_var_dyn,
		     &score_delta,
		     f_str[bbp->frame],m_msp->pstat_void,m_msp->show_code + SHOW_ANNOT_FULL);

	if (lc > 0) {
	  percent = (100.0*(float)l_aln_p->nident)/(float)lc;
	  ng_percent = (100.0*(float)l_aln_p->nident)/(float)(lc-(l_aln_p->ngap_q + l_aln_p->ngap_l));
	}
	else { percent = ng_percent = -1.00; }

	fprintf (fp, "a {\n");
	if (annot_var_dyn->string[0]) {
	  bp = annot_var_dyn->string;
	  while ((bp1=strchr(bp, '\n'))) {
	    *bp1 = '\0';
	    fprintf (fp, "# %s\n", bp);
	    *bp1 = '\n';
	    bp = bp1 + 1;
	  }
	}
	fprintf (fp, "  s %d %.1f\n", lsw_score, lbits);
	do_lav(fp, l_aln_p, align_code_dyn->string, percent, 0);

	if (ppst->nseq == 1) {
	  fprintf (fp, "a {\n");
	  fprintf (fp, "  s %d %.1f\n", lsw_score, lbits);
	  do_lav(fp, l_aln_p, align_code_dyn->string, percent, 1);
	}

	cur_ares_p = cur_ares_p->next;
	continue;
      }
#endif		/* ifdef LALIGN */

      NULL_dyn_string(annot_var_dyn);
      NULL_dyn_string(align_code_dyn);

      nc=calc_cons_a(aa0[bbp->frame],m_msp->n0, aa1, n1,
		     &lc,l_aln_p, cur_ares_p, ppst, 
		     seqc0, seqc1, seqca, cumm_seq_score,
		     m_msp->ann_arr,
		     m_msp->aa0a, m_msp->annot_p, seqc0a,
		     aa1a,  bbp->seq->annot_p, seqc1a,
		     &score_delta,
		     annot_var_dyn,
		     f_str[bbp->frame],
		     m_msp->pstat_void
		     );

      if (cur_ares_p->score_delta > 0) score_delta -= cur_ares_p->score_delta;

      disp_percent = percent = calc_fpercent_id(100.0, l_aln_p->nident,lc,m_msp->tot_ident, -1.0);
      disp_similar = calc_fpercent_id(100.0, l_aln_p->nsim, lc, m_msp->tot_ident, -1.0);
      disp_alen = lc;

      ngap = l_aln_p->ngap_q + l_aln_p->ngap_l;
      ng_percent = calc_fpercent_id(100.0, l_aln_p->nident,lc-ngap,m_msp->tot_ident, -1.0);
      if (m_msp->blast_ident) { 
	disp_percent = ng_percent;
	disp_similar = calc_fpercent_id(100.0, l_aln_p->npos, lc-ngap, m_msp->tot_ident, -1.0);
	disp_alen = lc - ngap;
      }

#ifndef SHOWSIM
      gpercent = ng_percent;
#else
      gpercent = disp_similar;
#endif

      lsw_score = cur_ares_p->sw_score + score_delta;
      /* removed 'first_line &&' so that LALIGN shows subject name/description */
      if (first_line && !(m_msp->markx&MX_M11OUT )) {
	if ((m_msp->markx & MX_ATYPE)!=7 && !(m_msp->markx & MX_M10FORM)) {
	  if (m_msp->markx & MX_MBLAST) {
	    /* provides >>id  description (length) line */
	    fprintf (fp, fmt,bline_p,annot_var_dyn->string,n1tot);
	  }
	  else {
	    fprintf (fp, fmt,bline_p,n1tot);
	  }
	}
	else if (m_msp->markx & MX_M10FORM) {
	  if (annot_var_dyn->string[0]) {	/* have annotation with '\n', replace with ';' in copy */
	    if ((annot_var_s10=(char *)calloc(strlen(annot_var_dyn->string)+1,sizeof(char)))==NULL) {
	      fprintf(stderr," ***error*** [%s:%d] -m 10 cannot allocate annot_var_s10[%d]\n",
		      __FILE__,__LINE__,(int)strlen(annot_var_dyn->string));
	      fprintf (fp,">>%s\n",bline_p);
	    }
	    else {
	      SAFE_STRNCPY(annot_var_s10,annot_var_dyn->string,strlen(annot_var_dyn->string));
	      bp = annot_var_s10;
	      while ((bp = strchr(bp+1,'\n'))) {
		if (bp[-1] != ';') {*bp = ';';}
		else {*bp = ' ';}
	      }
	      fprintf (fp,">>%s;%s\n",bline_p, annot_var_s10);
	      free(annot_var_s10); annot_var_s10 = NULL;
	    }
	  }
	  else {
	    fprintf (fp,">>%s\n",bline_p);
	  }
	}
      }

      /* this is required because cur_ares_p can carry scores from an
	 alignment without -S low-complexity re-scored using low
	 complexity */
#ifndef LALIGN
      if (first_line) {
	rst_p = &bbp->rst;
	first_line = 0;
      }
      else {
	rst_p = &cur_ares_p->rst;
      }
#else
      /* ensures that LALIGN alignments do not report 100% match
	 values */
      rst_p = &cur_ares_p->rst;
      first_line = 0;
#endif

      l_score0 = rst_p->score[ppst->score_ix] + score_delta;

      if (max(strlen(seqc0),strlen(seqc1)) > nc) {
	fprintf(stderr," mshowalign: nc/maxc: %d/%d seqc0/1: %lu/%lu\n",
		nc,maxc,strlen(seqc0),strlen(seqc1));
      }

#ifdef DEBUG
      /*
	if (lsw_score < bbp->rst.score[ppst->score_ix]) {
	fprintf(stderr," *** Warning - SW score=%d < opt score=%d ***\n",
	lsw_score, bbp->rst.score[ppst->score_ix]);
	}
      */
#endif
      calc_coord(m_msp->n0,bbp->seq->n1,
		m_msp->q_offset+(m_msp->q_off-1)+(m_msp->sq0off-1),
		bbp->seq->l_offset+(bbp->seq->l_off-1)+(m_msp->sq1off-1),
		l_aln_p);

      lzscore = find_z(l_score0, rst_p->escore, bbp->seq->n1, rst_p->comp, m_msp->pstat_void);
      if (ppst->zsflag > 20) {
	lzscore2 = find_z(l_score0, rst_p->escore, bbp->seq->n1, rst_p->comp, m_msp->pstat_void2);
      }
      lbits = zs_to_bit(lzscore, m_msp->n0, bbp->seq->n1);

      if (m_msp->markx & MX_MBLAST) {
	fprintf(fp, "\n Score = %.1f bits (%d),  Expect = %.1g\n",
		lbits, rst_p->score[ppst->score_ix] + score_delta,
		zs_to_E(lzscore, bbp->seq->n1, ppst->dnaseq, ppst->zdb_size, m_msp->db));

	fprintf(fp, " Identities = %d/%d (%d%%)", l_aln_p->nident, lc-ngap,
		(int)((100.0*(float)l_aln_p->nident+0.5)/(float)(lc-ngap)));

	if (!disp_dna_align) {
	  fprintf(fp, ", Positives = %d/%d (%d%%)", l_aln_p->npos, lc-ngap, 
		  (int)((100.0*(float)l_aln_p->npos+0.5)/(float)(lc-ngap)));
	}

	fprintf(fp, ", Gaps = %d/%d (%d%%)\n",ngap, lc-ngap,
		(int)((100.0*(float)(l_aln_p->ngap_q+l_aln_p->ngap_l)+0.5)/(float)(lc-ngap)));

	if (disp_dna_align) {
	  if (m_msp->qframe > 1 && bbp->frame > 0) {
	    fprintf (fp, " Strand=Minus/Plus\n");
	  }
	  else {
	    fprintf (fp, " Strand=Plus/Plus\n");
	  }
	}
	else {
	  if (m_msp->nframe > 2) {
	    fprintf (fp, " Frame= %d\n",bbp->frame+1);
	  }
	  else if (m_msp->nframe > 1) {
	    fprintf (fp, " Frame = %s\n",(bbp->frame>0 ? "Reverse" : "Forward"));
	  }
	  else if (m_msp->qframe > 1) {
	    fprintf (fp, " Frame = %s\n",(bbp->frame>0 ? "Reverse" : "Forward"));
	  }
	}
      }
      else if ((m_msp->markx & MX_ATYPE)!=7 && !(m_msp->markx & MX_M10FORM)) {
	if (annot_var_dyn->string[0]) {
	  if (m_msp->markx & MX_HTML) {
	    fprintf(fp,"<!-- ANNOT_START \"%s\" -->",link_name);}
	  /* ensure that last character is "\n" */
	  if (!m_msp->m8_show_annot) {
	    if (annot_var_dyn->string[strlen(annot_var_dyn->string)-1] != '\n') {
	      annot_var_dyn->string[strlen(annot_var_dyn->string)-1] = '\n';
	    }
	    fputs(annot_var_dyn->string, fp);
	  }
	  else { fputs("\n",fp);}

	  if (m_msp->markx & MX_HTML) {fputs("<!-- ANNOT_STOP -->",fp);}
	}

	/* this code makes sense for library searches, but not for
	   multiple non-intersecting alignments */

#ifndef LALIGN
	if (m_msp->nframe > 2) 
	  fprintf (fp, "Frame: %d",bbp->frame+1);
	else if (m_msp->nframe > 1) 
	  fprintf (fp, "Frame: %c",(bbp->frame? 'r': 'f'));
	else if (m_msp->qframe >= 0 && bbp->frame > 0 ) {
	  fputs("rev-comp",fp);
	  name0[nml-1]='\0';
	  if (!(m_msp->markx & MX_MBLAST)) SAFE_STRNCAT(name0,"-",sizeof(name0)-1);
	}

	if (m_msp->arelv > 0)
	  fprintf (fp, " %s: %3d", m_msp->alab[0],rst_p->score[0] + score_delta);
	if (m_msp->arelv > 1)
	  fprintf (fp, " %s: %3d", m_msp->alab[1],rst_p->score[1] + score_delta);
	if (m_msp->arelv > 2)
	  fprintf (fp, " %s: %3d", m_msp->alab[2],rst_p->score[2] + score_delta);
	fprintf (fp,"%s",info_str);
	if (ppst->zsflag>=0) {
	  fprintf (fp, "  Z-score: %4.1f  bits: %3.1f %sE(%ld): %4.2g%s", 
		   lzscore, lbits, 
		   html_pre_E, ppst->zdb_size,
		   zs_to_E(lzscore, bbp->seq->n1, ppst->dnaseq, ppst->zdb_size, m_msp->db),
		   html_post_E);
	  if (ppst->zsflag > 20) {
	    fprintf(fp," E2(): %4.2g",zs_to_E(lzscore2, bbp->seq->n1, ppst->dnaseq, ppst->zdb_size, m_msp->db));
	  }
	}
	fprintf (fp, "\n");

#else /* LALIGN */
	if ((m_msp->markx & MX_M11OUT) == 0) {
	  fprintf (fp, " %s score: %d; ", m_msp->alabel, lsw_score);
	  fprintf (fp," %3.1f bits; E(%ld) <  %.2g\n", lbits, ppst->zdb_size,
		   zs_to_E(lzscore, bbp->seq->n1, ppst->dnaseq, ppst->zdb_size, m_msp->db));
	}
#endif
      }
      else if (m_msp->markx & MX_M10FORM) {
#ifndef LALIGN
	if (m_msp->qframe > -1) {
	  if (m_msp->nframe > 2) {
	    fprintf(fp,"; %s_frame: %d\n",m_msp->f_id0,bbp->frame+1);
	  }
	  else {
	    fprintf(fp,"; %s_frame: %c\n",m_msp->f_id0,(bbp->frame > 0? 'r':'f'));
	  }
	}
	fprintf (fp, "; %s_%s: %3d\n", m_msp->f_id0,m_msp->alab[0],cur_ares_p->rst.score[0]+score_delta);
	if (m_msp->arelv > 1)
	  fprintf (fp,"; %s_%s: %3d\n", m_msp->f_id0,m_msp->alab[1],cur_ares_p->rst.score[1]+score_delta);
	if (m_msp->arelv > 2)
	  fprintf (fp,"; %s_%s: %3d\n", m_msp->f_id0,m_msp->alab[2],cur_ares_p->rst.score[2]+score_delta);
	if (info_str[0]) fprintf(fp,"; %s_info: %s\n",m_msp->f_id0,info_str);
	if (ppst->zsflag>=0) 
	  fprintf (fp,"; %s_z-score: %4.1f\n; %s_bits: %3.1f\n; %s_expect: %6.2g\n",
		   m_msp->f_id0,lzscore,
		   m_msp->f_id0,zs_to_bit(lzscore, m_msp->n0, bbp->seq->n1),
		   m_msp->f_id0,zs_to_E(lzscore, bbp->seq->n1, ppst->dnaseq, ppst->zdb_size, m_msp->db));
#else
	if ((m_msp->markx & MX_M11OUT) == 0) {
	  fprintf (fp,"; %s_%s: %d\n", m_msp->f_id0, m_msp->alab[0], lsw_score);
	  fprintf (fp,"; %s_z-score: %4.1f\n; %s_bits: %3.1f\n; %s_expect: %6.2g\n",
		   m_msp->f_id0, lzscore, m_msp->f_id0, lbits, m_msp->f_id0, 
		   zs_to_E(lzscore, bbp->seq->n1, ppst->dnaseq, ppst->zdb_size, m_msp->db));
	}
#endif	
      }

      do_show(fp, m_msp->n0, bbp->seq->n1, lsw_score, name0, name1, nml,
	      link_name, 
	      m_msp, ppst, seqc0, seqc0a, seqc1, seqc1a, seqca, cumm_seq_score,
	      nc, disp_percent, gpercent, disp_alen, l_aln_p, annot_var_dyn->string,
	      m_msp->annot_p, bbp->seq->annot_p);

      /* display the encoded alignment left over from showbest()*/

      if ((m_msp->markx & MX_M10FORM) &&
	  (m_msp->markx & MX_M9SUMM) && 
	  ((m_msp->show_code & SHOW_CODE_ALIGN) == SHOW_CODE_ALIGN)) {

	seq_code = cur_ares_p->aln_code;
	seq_code_len = cur_ares_p->aln_code_n;
	annot_code = cur_ares_p->annot_code;
	annot_code_len = cur_ares_p->annot_code_n;

	if (seq_code_len > 0 && seq_code != NULL) {
	  fprintf(fp,"; al_code: %s\n",seq_code);
	  /* 	  free(seq_code);  -- this is now freed in comp_lib2.c */
	  if (annot_code_len > 0 && annot_code != NULL) {
	    fprintf(fp,"; al_code_ann: %s\n",annot_code);
	    /* 	    free(ann_code);  -- this is now freed in comp_lib2.c */
	  }
	}
      }

      if (m_msp->markx & MX_HTML) fprintf(fp,"</pre><hr />");
      fflush(fp);

      freeseq(&seqc0,&seqc1, &seqca);
      freeseq_ann(&seqc0a, &seqc1a);

      cur_ares_p = cur_ares_p->next;

      if (cur_ares_p != NULL) {
	if (m_msp->markx & MX_HTML) {
	  sprintf(link_name,"%s_%d",l_name, cur_ares_p->index);
	  fprintf (fp,"<a name=\"%s\"><pre>",link_name);
	}
	else {
	  if (!(m_msp->markx & MX_MBLAST)) fprintf(fp,">--\n");
	}
      }		/* done finishing up  */
    }		/* while (cur_ares_p) */
    /* we are done displaying the alignment - be sure to free a_res memory */
  }

  free_dyn_string(annot_var_dyn);

  if (!(m_msp->markx & (MX_M8OUT+MX_HTML))) fprintf(fp,"\n");
}

void do_show(FILE *fp, int n0,int n1, int score,
	     char *name0, char *name1, int nml, char *link_name,
	     const struct mngmsg *m_msp, const struct pstruct *ppst,
	     char *seqc0, char *seqc0a,  char *seqc1, char *seqc1a,
	     char *seqca, int *cumm_seq_score, int nc,
	     float percent, float gpercent, int lc,
	     struct a_struct *aln, const char *annot_var_s,
	     const struct annot_str * q_annot_p,
	     const struct annot_str * l_annot_p)
{
  int tmp;

  if (m_msp->markx & MX_AMAP && (m_msp->markx & MX_ATYPE)==7)
    /* show text graphic of alignment (very rarely used) */
    disgraph(fp, n0, n1, percent, score,
	     aln->amin0, aln->amin1, aln->amax0, aln->amax1, m_msp->sq0off,
	     name0, name1, nml, aln->llen, m_msp->markx);
  else if (m_msp->markx & MX_M10FORM) {
    /* old tagged/parse-able format */
    if (ppst->sw_flag && m_msp->arelv>0)
      fprintf(fp,"; %s_score: %d\n",m_msp->f_id1,score);
    fprintf(fp,"; %s_ident: %5.3f\n",m_msp->f_id1,percent/100.0);
#ifndef SHOWSIM
    fprintf(fp,"; %s_gident: %5.3f\n",m_msp->f_id1,gpercent/100.0);
#else
    fprintf(fp,"; %s_sim: %5.3f\n",m_msp->f_id1,gpercent/100.0);
#endif

    fprintf(fp,"; %s_overlap: %d\n",m_msp->f_id1,lc);
    discons(fp, m_msp,
	    seqc0, seqc0a, seqc1, seqc1a, seqca, cumm_seq_score, nc,
	    n0, n1, name0, name1, nml, aln);
  }
  else  {    /* all "normal" alignment formats */
    if (!(m_msp->markx & MX_MBLAST)) {
#ifndef LALIGN
      fprintf(fp,"%s score: %d; ",m_msp->alabel, score);
#endif
#ifndef SHOWSIM
      fprintf(fp,"%4.1f%% identity (%4.1f%% ungapped) in %d %s overlap (%ld-%ld:%ld-%ld)\n",
	      percent,gpercent,lc,m_msp->sqnam,aln->d_start0,aln->d_stop0,
	      aln->d_start1,aln->d_stop1);
#else
      fprintf(fp,"%4.1f%% identity (%4.1f%% similar) in %d %s overlap (%ld-%ld:%ld-%ld)\n",
	      percent,gpercent,lc,m_msp->sqnam,aln->d_start0,aln->d_stop0,
	      aln->d_start1,aln->d_stop1);
#endif
    }

    if (m_msp->markx & MX_HTML) {
      do_url1(fp, m_msp, ppst, link_name,n1, aln,
	      annot_var_s, q_annot_p, l_annot_p);
    }

    if ((m_msp->markx & MX_AMAP) && ((m_msp->markx & MX_ATYPE)!=MX_ATYPE)) {
      fputc('\n',fp);
      tmp = n0;

      if (m_msp->qdnaseq == SEQT_DNA && m_msp->ldb_info.ldnaseq== SEQT_PROT)
	tmp /= 3;

      disgraph(fp, tmp, n1, percent, score,
	       aln->amin0, aln->amin1,
	       aln->amax0, aln->amax1,
	       m_msp->sq0off,
	       name0, name1, nml, aln->llen,m_msp->markx);
    }

    if (m_msp->markx & MX_HTML) {
      fprintf(fp, "<!-- ALIGN_START \"%s\" -->",link_name);
    }
    discons(fp, m_msp,
	    seqc0, seqc0a, seqc1, seqc1a, seqca, cumm_seq_score, nc,
	    n0, n1, name0, name1, nml, aln);
    if (m_msp->markx & MX_HTML) {fputs("<!-- ALIGN_STOP -->",fp);}
    fputc('\n',fp);

  }
}

void 
do_lav(FILE *fp, struct a_struct *aln_p, char *seqc,
       float percent, int is_mirror) {
  int cur_b0, cur_b1, cur_e0, cur_e1;
  int ipercent;
  long len;
  char *seqc_p, *num_e;

  ipercent = (int)(percent+0.5);

  cur_b0 = aln_p->d_start0;
  cur_b1 = aln_p->d_start1;
  cur_e0 = aln_p->d_stop0;
  cur_e1 = aln_p->d_stop1;

  if (!is_mirror) {
    fprintf (fp, "  b %d %d\n  e %d %d\n", 
	     cur_b0, cur_b1, cur_e0, cur_e1);
  }
  else  {
    fprintf (fp, "  b %d %d\n  e %d %d\n", 
	     cur_b1, cur_b0, cur_e1, cur_e0);
  }

  seqc_p = seqc;
  
  while (*seqc_p) {
    if (*seqc_p == '=') {	/* extend match in both sequences */
      len = strtol(seqc_p+1, &num_e, 10);
      cur_e0 = cur_b0 + len - 1;
      cur_e1 = cur_b1 + len - 1;
      if (!is_mirror) {
	fprintf(fp, "  l %d %d %d %d %d\n",
		cur_b0, cur_b1, cur_e0, cur_e1,
		ipercent);
      }
      else  {
	fprintf(fp, "  l %d %d %d %d %d\n",
		cur_b1, cur_b0, cur_e1, cur_e0,
		ipercent);
      }
      cur_b0 = cur_e0 + 1;
      cur_b1 = cur_e1 + 1;
    }
    else if (*seqc_p == '+') {	/* extend insertion in seq0  by incrementing seq1 */
      len = strtol(seqc_p+1, &num_e, 10);
      cur_b1 += len;
    }
    else {			/* extend insertion in seq1 by incrementing seq0 */
      len = strtol(seqc_p+1, &num_e, 10);
      cur_b0 += len;
    }
    seqc_p = num_e;
  }

  fprintf (fp, "}\n");
}

void	/* initialize consensus arrays */
initseq(char **seqc0, char **seqc1, char **seqca, int seqsiz)
{
  *seqc0=(char *)calloc((size_t)seqsiz*3,sizeof(char));
  if (*seqc0==NULL)
    {fprintf(stderr,"cannot allocate consensus arrays %d\n",seqsiz);
     exit(1);}
  *seqc1=*seqc0 + seqsiz;
  *seqca=*seqc1 + seqsiz;
}

void freeseq(char **seqc0, char **seqc1, char **seqca)
{
  free(*seqc0);
  *seqc0 = *seqc1 = *seqca = NULL;
}

void	/* initialize consensus annotation arrays */
initseq_ann(char **seqc0a, char **seqc1a, int seqsiz)
{
  *seqc0a=(char *)calloc((size_t)seqsiz*5,sizeof(char));
  if (*seqc0a==NULL)
    {fprintf(stderr,"cannot allocate consensus arrays %d\n",seqsiz);
     exit(1);}
  *seqc1a=*seqc0a + seqsiz;
}

void freeseq_ann(char **seqc0a, char **seqc1a)
{
  if (*seqc0a != NULL) {
    free(*seqc0a);
    *seqc0a = *seqc1a = NULL;
  }
}

#include <math.h>
/* calculates percentages, optionally ensuring that 100% is completely
   identical*/
float calc_fpercent_id(float scale, int n_ident, int n_alen, int tot_ident, float fail) {
  float f_id, f_decimal;
  /* int n_sig; */

  if (n_alen <= 0) { return fail;}

  /*
  n_sig = 3;
  n_sig = tot_ident;
  if (tot_ident==1) { n_sig = 3;}
  */

  f_id = (float)n_ident/(float)n_alen;

  if (tot_ident && n_ident != n_alen) {
    f_decimal = 0.999;
    /*
    if (n_sig == 4) {f_decimal = 0.9999;}
    else if (n_sig == 3) { f_decimal = 0.999;}
    else if (n_sig == 5) { f_decimal = 0.99999;}
    else {f_decimal = 1.0 - powf(0.1, n_sig);}
    */
    if (f_id > f_decimal) f_id = f_decimal;
  }

  return scale*f_id;
}
