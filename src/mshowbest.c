/* $Id: mshowbest.c 1281 2014-08-21 17:32:06Z wrp $ */

/* copyright (c) 1996, 1997, 1998, 1999, 2014 by William R. Pearson and
   The Rector and Visitors of the University of Virginia */

/* Licensed under the Apache License, Version 2.0 (the "License"); you
   may not use this file except in compliance with the License.  You
   may obtain a copy of the License at

   http://www.apache.org/licenses/LICENSE-2.0

   Unless required by applicable law or agreed to in writing,
   software distributed under this License is distributed on an "AS
   IS" BASIS, WITHOUT WRRANTIES OR CONDITIONS OF ANY KIND, either
   express or implied.  See the License for the specific language
   governing permissions and limitations under the License. 
*/

/* 2-April-2009 changes to simplify interactive display logic.  Coming
   into showbest(), things are interactive (quiet==0) or use
   m_msg.nshow */

/*   29-Oct-2003 - changes so that bbp->mseq->cont < 0 => aa1 sequence is
     already in aa1, no re_openlib or re_getlib required
*/

/*   14-May-2003 Changes to use a more consistent coordinate numbering
     system for displays.  aln->d_start[01] is now consistently used
     to report the start of the alignment in all functions, and
     mshowbest.c has been modified to use d_start[01] instead of
     d_start[01]-1.  aln->min[01] now starts at 0 for all functions;
     instead of 1 for some functions (dropnfa.c, dropgsw.c, dropfs2.c
     earlier).
*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "defs.h"
#include "structs.h"
#include "param.h"
#include "mm_file.h"
#include "best_stats.h"

#define MAX_BLINE 256

/* function calls necessary to re_getlib() the sequence and, do
   alignments, if necessary
*/

#define RANLIB (m_fptr->ranlib)


int
re_getlib(unsigned char *, struct annot_str **, 
	  int, int, int, int, int, long *, long *, 
	  struct lmf_str *m_fptr);

#include "drop_func.h"

void
s_annot_to_aa1a(int n1, struct annot_str *annot_p, unsigned char *ann_arr);

struct a_res_str *
build_ares_code(unsigned char *aa0, int n0,
		unsigned char *aa1, struct seq_record *seq,
		int frame, int *have_ares, int repeat_thresh, 
		struct mngmsg *m_msp, struct pstruct *ppst,
		void *f_str
		);

struct lmf_str *re_openlib(struct lmf_str *, int outtty);

extern void calc_coord(int n0, int n1, long qoffset, long loffset,
		      struct a_struct *aln);

extern float
calc_fpercent_id(float scale, int n_ident, int n_alen, int tot_ident, float fail);

extern int
get_annot(char *sname, struct mngmsg *m_msp, char *bline, int n1, struct annot_str **annot_p,
 	  int target, int debug);
extern double find_z(int score, double escore, int length, double comp,void *);
extern double zs_to_E(double zs, int n1, int dnaseq, long db_size, struct db_str db);
extern double zs_to_bit(double, int, int);
extern int E1_to_s(double e_val, int n0, int n1, int db_size, void *pu);

void header_aux(FILE *);
void show_aux(FILE *, struct beststr *);
void w_abort (char *p, char *p1);

extern double zs_to_bit(double, int, int);

/* showbest() shows a list of high scoring sequence descriptions, and
   their rst.scores.  If -m 9, then an additional complete set of
   alignment information is provided.

   If PCOMPLIB or m_msg.quiet then the number of high scores to be
   shown is pre-determined by m_msg.mshow before showbest is called.

   The comp_lib2.c version re_getlib()'s the sequence for its
   discription, and then does another alignment for -m 9 (Thus, it
   needs an f_str.  The PCOMPLIB version has everything available in
   beststr before showbest() is called.
*/

void showbest (FILE *fp, unsigned char **aa0, unsigned char *aa1save, int maxn,
	       struct beststr **bptr,int nbest,
	       int qlib, struct mngmsg *m_msp,
	       struct pstruct *ppst, struct db_str db,
	       char **info_gstring2
	       ,void **f_str
)
{
  unsigned char *aa1;
  int best_align_done = 0;
  int ntmp = 0;
  char bline[MAX_BLINE], fmt[40], pad[MAX_BLINE], fmt2[40], rline[40];
  char l_name[128], link_name[140];
  int istart = 0, istop, ib;
  int nshow;		/* number of sequences shown before prompt,
			   and ultimately displayed */
  int first_line, link_shown;
  int quiet;
  int r_margin;
  struct beststr *bbp;
  int n1tot;
  char *bp, *bline_p;
  char rel_label[12];
  char score_label[120];
  char tmp_str[20], *seq_code, *annot_str;
  int seq_code_len, annot_str_len;
  long loffset;		/* loffset is offset from beginning of real sequence */
  long l_off;		/* l_off is the the virtual coordinate of residue 1 */
  int n1, ranlib_done;
  struct rstruct rst;
  int l_score0, ngap;
  double lzscore, lzscore2, lbits;
  float percent, gpercent, ng_percent;
  struct a_struct *aln_p;
  struct a_res_str *cur_ares_p;
  struct rstruct *rst_p;
  int gi_num;
  char html_pre_E[120], html_post_E[120];
  int have_lalign = 0;

  struct lmf_str *m_fptr;

  /* for lalign alignments, only show stuff when -m != 11 */

  if (m_msp->markx & MX_M11OUT) return;
  if (strcmp(m_msp->label,"ls-w")==0) {
    have_lalign = 1;
    if ((m_msp->markx & MX_M9SUMM) == 0) return;
  }

  rel_label[0]='\0';
  SAFE_STRNCPY(score_label,"scores", sizeof(score_label));

  quiet = m_msp->quiet;

  if (m_msp->aln.llen > MAX_BLINE) m_msp->aln.llen = MAX_BLINE;

  if (ppst->zsflag < 0) r_margin = 10;
  else if (ppst->zsflag>=0  && m_msp->srelv > 1 ) r_margin = 19;
  else r_margin = 10;

  if (m_msp->markx & MX_M9SUMM && (m_msp->show_code == SHOW_CODE_ID || m_msp->show_code == SHOW_CODE_IDD)) {
#ifdef SHOWSIM
    r_margin += 15;
#else
    r_margin += 10;
#endif
  }
  else if (m_msp->markx & MX_MBLAST2) {
    r_margin -= 10;
  }
  else if (m_msp->markx & (MX_M9SUMM + MX_M8OUT)) {
    r_margin = 0;
  }

  if (m_msp->markx & MX_HTML) {
    strncpy(html_pre_E,"<font color=\"darkred\">",sizeof(html_pre_E));
    strncpy(html_post_E,"</font>",sizeof(html_post_E));

  }
  else {
    html_pre_E[0] = html_post_E[0] = '\0';
  }

  if (m_msp->nframe < 0) {
    sprintf(fmt,"%%-%ds (%%4d)",m_msp->aln.llen-r_margin);
  }
  else {
    sprintf(fmt,"%%-%ds (%%4d)",m_msp->aln.llen-(r_margin+4));
  }
  sprintf(fmt2,"%%-%ds",m_msp->aln.llen-r_margin+8);

  memset(pad,' ',m_msp->aln.llen-(r_margin+6));
  pad[m_msp->aln.llen-(r_margin+12)]='\0';
  if (have_lalign) {
    if (ppst->show_ident) {
      SAFE_STRNCPY(score_label,"alignments", sizeof(score_label));
      pad[m_msp->aln.llen-(r_margin+16)]='\0';
    }
    else {
      SAFE_STRNCPY(score_label,"non-identical alignments", sizeof(score_label));
      pad[m_msp->aln.llen-(r_margin+30)]='\0';
    }
  }

  nshow = min(m_msp->nshow,nbest);

  if ((bp = strchr (m_msp->qtitle, '\n')) != NULL) *bp = '\0';
  if (m_msp->markx & MX_M8OUT) {
    if ((bp = strchr (m_msp->qtitle, ' ')) != NULL) *bp = '\0';
  }

/*   fprintf (fp, "%3d %s\n", qlib,m_msp->qtitle); */

  if (m_msp->markx & MX_HTML) fprintf(fp,"<pre>");

  /* **************************************************************** */
  /* done with display format */
  /* **************************************************************** */

  /* **************************************************************** */
  /* prompt for number of best scores if quiet == 0 */
  /* **************************************************************** */

  if (quiet == 0) {	/* interactive */
    nshow = min(m_msp->nshow, nbest);
    printf(" How many scores would you like to see? [%d] ",nshow);
    fflush(stdout);
    if (fgets(rline,20,stdin)==NULL) exit(0);
    if (rline[0]!='\n' && rline[0]!=0) sscanf(rline,"%d",&nshow);
    if (nshow > nbest) nshow=nbest;
    if (nshow<=0) nshow = min(20,nbest);
  }

  /* display number of hits for -m 8C (Blast Tab-commented format) */
  if (m_msp->markx & MX_M8COMMENT) {
    /* line below copied from BLAST+ output */
    fprintf(fp,"# Fields: query id, subject id, %% identity, alignment length, mismatches, gap opens, q. start, q. end, s. start, s. end, evalue, bit score");
    if (ppst->zsflag > 20) {fprintf(fp,", eval2");}
    if (m_msp->show_code & (SHOW_CODE_ALIGN+SHOW_CODE_CIGAR)) { fprintf(fp,", aln_code");}
    else if ((m_msp->show_code & SHOW_CODE_BTOP)==SHOW_CODE_BTOP) { fprintf(fp,", BTOP");}

    fprintf(fp,"\n");
    fprintf(fp,"# %d hits found\n",nshow);
  }

  /* **************************************************************** */
  /* have number of scores in interactive or quiet mode */
  /* display "The best scores are" */
  /* **************************************************************** */

  if (m_msp->markx & MX_MBLAST2) {
    fprintf(fp, "%81s\n"," Score     E");
    fprintf(fp, "Sequences producing significant alignments:                          (Bits)  Value\n\n");
  }
  else if (!(m_msp->markx & MX_M8OUT)) {
    if (ppst->zsflag >= 0) {
      if (m_msp->z_bits==1) {/* show bit score */
	fprintf(fp,"\nThe best%s %s are:%s%s bits %sE(%ld)%s",
		rel_label,score_label,pad,m_msp->label,html_pre_E,ppst->zdb_size,html_post_E);
	if (ppst->zsflag > 20) {
	  fprintf(fp," E2()");
	}
      }
      else {/* show z-score */
	fprintf(fp,"\nThe best%s %s are:%s%s z-sc %sE(%ld)%s",
		rel_label,score_label,pad,m_msp->label,html_pre_E,ppst->zdb_size,html_post_E);
	if (ppst->zsflag > 20) {
	  fprintf(fp," E2()");
	}
      }
      header_aux(fp);
      if (m_msp->markx & MX_M9SUMM) {
	if (m_msp->show_code == SHOW_CODE_ID || m_msp->show_code == SHOW_CODE_IDD) {
#ifdef SHOWSIM
	  fprintf(fp," %%_id  %%_sim  alen");
#else
	  fprintf(fp," %%_id  alen");
#endif
	}
	else {
	  if (m_msp->markx & MX_HTML && m_msp->show_code != SHOW_CODE_ID && m_msp->show_code != SHOW_CODE_IDD) { fprintf(fp,"<!-- ");}
#ifndef SHOWSIM
	  fprintf(fp,"\t%%_id  %%_gid %4s  alen  an0  ax0  pn0  px0  an1  ax1 pn1 px1 gapq gapl  fs ",m_msp->f_id1);
#else
	  fprintf(fp,"\t%%_id  %%_sim %4s  alen  an0  ax0  pn0  px0  an1  ax1 pn1 px1 gapq gapl  fs ",m_msp->f_id1);
#endif
	}
	if (m_msp->show_code & (SHOW_CODE_ALIGN+SHOW_CODE_CIGAR)) { fprintf(fp," aln_code"); }
	if (m_msp->markx & MX_HTML && m_msp->show_code != SHOW_CODE_ID && m_msp->show_code != SHOW_CODE_IDD) { fprintf(fp," -->");}
      }
      fprintf(fp,"\n");
    }
    else {
      fprintf(fp,"\nThe best%s %s are:%s%s",rel_label,score_label,pad,m_msp->label);
      header_aux(fp);
      if (m_msp->markx & MX_M9SUMM) {
	if (m_msp->show_code == SHOW_CODE_ID || m_msp->show_code == SHOW_CODE_IDD) {
#ifdef SHOWSIM
	  fprintf(fp," %%_id  %%_sm  alen");
#else
	  fprintf(fp," %%_id  alen");
#endif
	}
	else {
#ifndef SHOWSIM
	  fprintf(fp,"\t%%_id  %%_gid %4s  alen  an0  ax0  pn0  px0  an1  ax1 pn1 px1 gapq gapl  fs ",m_msp->f_id1);
#else
	  fprintf(fp,"\t%%_id  %%_sim %4s  alen  an0  ax0  pn0  px0  an1  ax1 pn1 px1 gapq gapl  fs ",m_msp->f_id1);
#endif	/* SHOWSIM */
	}
      }
      if (m_msp->show_code & (SHOW_CODE_ALIGN+SHOW_CODE_CIGAR)) { fprintf(fp," aln_code"); }
      fprintf(fp,"\n");
    }
  }	/* !(m_msp->markx & MX_M8OUT) */

  istart = 0;
l1:
  istop = min(nshow, nbest);

  for (ib=istart; ib<istop; ib++) {
    bbp = bptr[ib];
    if (ppst->do_rep) {
      bbp->repeat_thresh = 
	min(E1_to_s(ppst->e_cut_r, m_msp->n0, bbp->seq->n1,ppst->zdb_size, m_msp->pstat_void),
	    bbp->rst.score[ppst->score_ix]);
    }

#ifdef DEBUG
    if (bbp->seq->n1 != bbp->n1 ) {
      fprintf(stderr, " *** lib len error [%d!=%d] *** %s score %d\n",
	      bbp->seq->n1,bbp->n1, bbp->mseq->libstr, bbp->rst.score[0]);
    }
#endif

    /* this gets us a valid bline[] and the library for searching if necessary
       do not read if we have a long enough bline or we don't need a sequence 
    */
    if (bbp->mseq->bline != NULL && bbp->mseq->bline_max >= m_msp->aln.llen) {
      ranlib_done = 0;

      /* copy m_msp->aln.llen, not llen-r_margin, because the r_margin
	 will be set later, possibly after the gi|12345 is removed */
      strncpy(bline,bbp->mseq->bline,m_msp->aln.llen);
      bline[m_msp->aln.llen]='\0';
    }
    else {
      if ((m_fptr=re_openlib(bbp->mseq->m_file_p,!m_msp->quiet))==NULL) {
	fprintf(stderr,"*** cannot re-open %s\n",bbp->mseq->m_file_p->lb_name);
	exit(1);
      }
      RANLIB(bline,m_msp->aln.llen,bbp->mseq->lseek,bbp->mseq->libstr,m_fptr);
      ranlib_done = 1;
    }

    /* get a valid cur_ares_p chain and put it in bbp->ares */
    if (!m_msp->align_done && (m_msp->stages>1 || (m_msp->markx & MX_M9SUMM))) {	/* we need a sequence */
      if (bbp->seq->aa1b == NULL || (m_msp->ann_flg==1 && bbp->seq->annot_p==NULL)) {
	if (!ranlib_done) {	/* we didn't open the library already */
	  if ((m_fptr=re_openlib(bbp->mseq->m_file_p,!m_msp->quiet))==NULL) {
	    fprintf(stderr,"*** cannot re-open %s\n",bbp->mseq->m_file_p->lb_name);
	    exit(1);
	  }
	  RANLIB(bline,m_msp->aln.llen,bbp->mseq->lseek,bbp->mseq->libstr,m_fptr);
	  ranlib_done = 1;
	}
	n1 = re_getlib(aa1save,
		       (m_msp->ann_flg==1) ? &(bbp->seq->annot_p) : NULL, 
		       maxn,m_msp->ldb_info.maxt3,
		       m_msp->ldb_info.l_overlap,bbp->mseq->cont,m_msp->ldb_info.term_code,
		       &bbp->seq->l_offset,&bbp->seq->l_off,bbp->mseq->m_file_p);

	aa1 = aa1save;

	if (m_msp->ann_flg==2 && bbp->seq->annot_p==NULL ) {
	  /* get information about this sequence from bline */
	  if (get_annot(m_msp->annot1_sname, m_msp, bline, bbp->seq->n1, &(bbp->seq->annot_p), 1, ppst->debug_lib) > 0) {
	    /* do something with annotation */
	    s_annot_to_aa1a(bbp->n1, bbp->seq->annot_p, m_msp->ann_arr);
	  }
	}
      }
      else {
	n1 = bbp->seq->n1;
	aa1 = bbp->seq->aa1b;
      }

      if (n1 != bbp->n1) {
	fprintf(stderr," *** sequence length conflict %d != %d: %s\n", n1, bbp->n1, bline);
	continue;
      }

      if ( m_msp->stages > 1 && bbp->rst.score[2] == -BIGNUM) { 
	/* this is not typically done unless m_msp->stages > 1 */
	do_opt (aa0[bbp->frame], m_msp->n0, aa1, n1, bbp->frame, ppst, f_str[bbp->frame], &rst);
	bbp->rst.score[2]=rst.score[2];
      }

      if (!bbp->have_ares & 0x1) {
	bbp->a_res = build_ares_code(aa0[bbp->frame], m_msp->n0, aa1, bbp->seq,
				     bbp->frame, &bbp->have_ares,
				     bbp->repeat_thresh, m_msp, ppst, f_str[bbp->frame] );
	best_align_done = 1;
      }
    }	/* end stages > 1 || MX_M9SUMM9 */

    n1tot = (bbp->mseq->n1tot_p) ? *bbp->mseq->n1tot_p : bbp->seq->n1;

    bline_p = bline;
    if (!(m_msp->markx & (MX_M8OUT)) && !strncmp(bline,"gi|",3)) {
      bline_p = strchr(bline+4,'|')+1;
      *(bline_p-1) = 0;
      gi_num = atoi(bline+3);
    }

  /* l_name is used to build an HTML link from the bestscore line to
     the alignment.  It can also be used to discriminate multiple hits
     from the same long sequence.  This requires that fast_pan use -m 6.

     (6-April-2013) Add ability to specify additional alignments with
     link_name;
  */

    SAFE_STRNCPY(l_name,bline_p,sizeof(l_name)); /* get rid of text after second "|" */
    if ((bp=strchr(l_name,' '))!=NULL) *bp=0;
    if ((bp=strchr(&l_name[6],'|'))!=NULL) *bp='\0'; 	/* increase to [6] from [3] to allow longer db names "ref", "unk", */
    if (m_msp->nframe > 2) sprintf(&l_name[strlen(l_name)],"_%d",bbp->frame+1);
    else if (m_msp->nframe > 0 && bbp->frame == 1)
      SAFE_STRNCAT(l_name,"_r",sizeof(l_name));
    if (bbp->mseq->cont-1 > 0) {
      sprintf(tmp_str,":%d",bbp->mseq->cont-1);
      SAFE_STRNCAT(l_name,tmp_str,sizeof(l_name));
    }

    if (m_msp->markx & MX_M8OUT) {
      if ((bp=strchr(bline_p,' '))!=NULL) *bp = '\0';
    }
    else {
      bline_p[m_msp->aln.llen-r_margin]='\0';
      /* check for translated frame info */
      if (m_msp->nframe > -1) bline_p[m_msp->aln.llen-(r_margin+4)]='\0';
    }
    /* now its time to report the summary numbers for all the alignments */

    /* in the next loop, cur_ares_p could be NULL if we haven't done do_walign() */
    cur_ares_p = bbp->a_res;

    first_line = 1;
    do {
      /* if cur_res_p != NULL, then we get rst from a_res->rst
	 Otherwise, it comes from bbp->rst
      */

      if ((!first_line || (have_lalign && !ppst->show_ident)) && cur_ares_p ) {
	rst_p = &cur_ares_p->rst;
      }
      else {
	rst_p = &bbp->rst;
      }

      n1 = bbp->seq->n1;
      l_score0 = rst_p->score[ppst->score_ix];
      lzscore = find_z(l_score0, rst_p->escore, n1, rst_p->comp, m_msp->pstat_void);
      if (ppst->zsflag > 20) {
	lzscore2 = find_z(l_score0, rst_p->escore, n1, rst_p->comp, m_msp->pstat_void2);
      }
      lbits = zs_to_bit(lzscore, m_msp->n0, n1);

      /* *********************************** */
      /* standard "The best scores are" here */
      /* *********************************** */

      if (!(m_msp->markx & (MX_M8OUT + MX_MBLAST2))) {
	if (first_line) {
	  first_line = 0;
	  fprintf (fp, fmt,bline_p,n1tot);
	  if (m_msp->nframe > 2) fprintf (fp, " [%d]", bbp->frame+1);
	  else if (m_msp->nframe >= 0) fprintf(fp," [%c]",(bbp->frame > 0 ?'r':'f'));
	}
	else {
	  fprintf (fp, fmt2,"\n+-");
	}

	if (m_msp->srelv == 1) fprintf (fp, " %4d", rst_p->score[ppst->score_ix]);
	else {
	  if (m_msp->srelv-1 > 0) fprintf (fp, " %4d", rst_p->score[0]);
	  if (m_msp->srelv-1 > 1 || m_msp->stages>1)
	    fprintf (fp, " %4d", rst_p->score[1]);
	  fprintf (fp, " %4d", rst_p->score[ppst->score_ix]);
	}

	if (ppst->zsflag>=0) { 
	  if (m_msp->z_bits==1) {
	    fprintf (fp, " %.1f %s%7.2g%s",lbits,html_pre_E,
		     zs_to_E(lzscore, n1, ppst->dnaseq, ppst->zdb_size, m_msp->db),
		     html_post_E);
	    if (ppst->zsflag > 20) {
	      fprintf (fp, " %7.2g",zs_to_E(lzscore2, n1, ppst->dnaseq, ppst->zdb_size, m_msp->db));
	    }
	  }
	  else {
	    fprintf (fp, " %.1f %s%7.2g%s",lzscore,html_pre_E,
		     zs_to_E(lzscore, n1, ppst->dnaseq, ppst->zdb_size, m_msp->db),
		     html_post_E);
	    if (ppst->zsflag > 20) {
	      fprintf (fp, " %7.2g",zs_to_E(lzscore2, n1, ppst->dnaseq, ppst->zdb_size, m_msp->db));
	    }
	  }
	}
	show_aux(fp,bbp);
      }
      else if (m_msp->markx & MX_M8OUT) {	/* MX_M8OUT -- provide query, library */
	if (first_line) {first_line = 0;}
	fprintf (fp,"%s\t%s",m_msp->qtitle,bline_p);
      }
      else if (m_msp->markx & MX_MBLAST2) {	/* blast "Sequences producing" */ 
	if (first_line) {first_line = 0;}
	fprintf (fp,"%-67s %6.1f    %.1g", bline_p, lbits,
		    zs_to_E(lzscore,n1,ppst->dnaseq,ppst->zdb_size,m_msp->db));
      }

      if (m_msp->markx & MX_M9SUMM || m_msp->markx & MX_M8OUT) {
	loffset = bbp->seq->l_offset;
	l_off = bbp->seq->l_off;
	aln_p = &cur_ares_p->aln;
	seq_code = cur_ares_p->aln_code;
	seq_code_len = cur_ares_p->aln_code_n;
	annot_str = cur_ares_p->annot_code;
	annot_str_len = cur_ares_p->annot_code_n;

	ngap = cur_ares_p->aln.ngap_q + cur_ares_p->aln.ngap_l;
        percent = calc_fpercent_id(100.0,aln_p->nident,aln_p->lc, m_msp->tot_ident, -100.0);
        ng_percent = calc_fpercent_id(100.0,aln_p->nident,aln_p->lc-ngap, m_msp->tot_ident, -100.0);

#ifndef SHOWSIM
	gpercent = calc_fpercent_id(100.0, aln_p->nident, aln_p->lc-ngap, m_msp->tot_ident, -100.0);
#else
	gpercent = calc_fpercent_id(100.0, cur_ares_p->aln.nsim, aln_p->lc, m_msp->tot_ident, -100.0);
#endif	/* SHOWSIM */

	if (m_msp->show_code != SHOW_CODE_ID && m_msp->show_code != SHOW_CODE_IDD) {	/* show more complete info than just identity */

	  /*  	calc_astruct(aln_p, cur_ares_p); -- this function
		should not be used after calc_code or any other
		alignment that calculates amax0/amax1 */

	  /* we need the coordinates for annotated SHOW_CODE_ALIGN */
	  calc_coord(m_msp->n0,bbp->seq->n1,
		     m_msp->q_offset + (m_msp->q_off-1) + (m_msp->sq0off-1),
		     loffset + (l_off-1) + (m_msp->sq1off-1),
		     aln_p);

	  /* if (m_msp->markx & MX_HTML) fprintf(fp,"<!-- "); */
	  /*            %_id  %_sim s-w alen an0  ax0  pn0  px0  an1  ax1  pn1  px1 gapq gapl fs  */
	  /*                    alignment    min  max            min  max */
	  /*                    sequence coordinate    min  max            min  max */
	  if (!(m_msp->markx & MX_M8OUT)) {
	    fprintf(fp,"\t%5.3f %5.3f %4d %4d %4ld %4ld %4ld %4ld %4ld %4ld %4ld %4ld %3d %3d %3d",
		    percent/100.0,gpercent/100.0, 
		    cur_ares_p->sw_score,
		    aln_p->lc,
		    aln_p->d_start0,aln_p->d_stop0,
		    aln_p->q_start_off, aln_p->q_end_off,
		    aln_p->d_start1,aln_p->d_stop1,
		    aln_p->l_start_off, aln_p->l_end_off,
		    aln_p->ngap_q,aln_p->ngap_l,aln_p->nfs);
	    if ((m_msp->show_code & (SHOW_CODE_ALIGN+SHOW_CODE_CIGAR+SHOW_CODE_BTOP))
		&& seq_code_len > 0 && seq_code != NULL) {
	      fprintf(fp,"\t%s",seq_code);
	      if (annot_str_len > 0 && annot_str != NULL) {
		fprintf(fp,"\t%s",annot_str);
	      }
	    }
	  }
	  else {	/* MX_M8OUT -- blast order, tab separated */
	    fprintf(fp,"\t%.2f\t%d\t%d\t%d\t%ld\t%ld\t%ld\t%ld\t%.2g\t%.1f",
		    ng_percent,aln_p->lc,aln_p->nmismatch,
		    aln_p->ngap_q + aln_p->ngap_l+aln_p->nfs,
		    aln_p->d_start0, aln_p->d_stop0,
		    aln_p->d_start1, aln_p->d_stop1,
		    zs_to_E(lzscore,n1,ppst->dnaseq,ppst->zdb_size,m_msp->db),
		    lbits);
	    if (ppst->zsflag > 20) {
	      fprintf(fp,"\t%.2g",zs_to_E(lzscore2, n1, ppst->dnaseq, ppst->zdb_size, m_msp->db));
	    }
	    if ((m_msp->show_code & (SHOW_CODE_ALIGN+SHOW_CODE_CIGAR+SHOW_CODE_BTOP)) && seq_code_len > 0 && seq_code != NULL) {
	      fprintf(fp,"\t%s",seq_code);
	      if (annot_str_len > 0 && annot_str != NULL) {
		fprintf(fp,"\t%s",annot_str);
	      }
	    }
	    fprintf(fp,"\n");
	  }
	}
	else {	/* !SHOW_CODE -> SHOW_ID or SHOW_IDD*/
#ifdef SHOWSIM
	  fprintf(fp," %5.3f %5.3f %4d", 
		  percent/100.0,
		  (float)aln_p->nsim/(float)aln_p->lc,aln_p->lc);
#else
	  fprintf(fp," %5.3f %4d", percent/100.0,aln_p->lc);
#endif
	  if (m_msp->markx & MX_HTML) {
	    if (cur_ares_p->index > 0) {
	      sprintf(link_name,"%s_%d",l_name, cur_ares_p->index);
	    }
	    else {
	      SAFE_STRNCPY(link_name, l_name, sizeof(l_name));
	    }
	    fprintf(fp," <a href=\"#%s\">align</a>",link_name);
	    link_shown = 1;
	  }
	  else { link_shown = 0;}

	  if ((m_msp->show_code & SHOW_CODE_ID) == SHOW_CODE_ID) {
	    annot_str = cur_ares_p->annot_var_id;
	  }
	  else if ((m_msp->show_code & SHOW_CODE_IDD) == SHOW_CODE_IDD) {
	    annot_str = cur_ares_p->annot_var_idd;
	  }
	  else {
	    annot_str = NULL;
	  }
	  if (annot_str && annot_str[0]) {
	    fprintf(fp," %s",annot_str);
	  }
	}
      }
    } while ( cur_ares_p && (cur_ares_p = cur_ares_p->next));

    /*    if ((m_msp->markx & MX_HTML) && !link_shown) fprintf(fp," <a href=\"#%s\">align</a>",l_name); */
    if (!(m_msp->markx & MX_M8OUT)) fprintf(fp, "\n");
    fflush(fp);
  }

  if (quiet==0) {
    printf(" More scores? [0] ");
    fflush(stdout);
    if (fgets(rline,20,stdin)==NULL) exit(0);
    ntmp = 0;
    if (rline[0]!='\n' && rline[0]!=0) sscanf(rline,"%d",&ntmp);
    if (ntmp<=0) ntmp = 0;
    if (ntmp>0) {
      istart = istop;
      nshow = min(nshow+ntmp, nbest);
      goto l1;
    }
  }	/* end of for (ib) loop */

  if (m_msp->markx & MX_MBLAST2) {fprintf(fp, "\n\n");}

  m_msp->nshow = nshow;	/* save the number of hits displayed for showalign */

  if (best_align_done) { m_msp->align_done = 1;}	/* note that alignments are done */

  if (m_msp->markx & MX_HTML) fprintf(fp,"</pre><hr>\n");
}
