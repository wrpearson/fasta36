/* $Id: compacc2.c 1280 2014-08-21 00:47:55Z wrp $ */

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

/* Concurrent read version */

#include <stdio.h>
#include <stdlib.h>
#if defined(UNIX)
#include <unistd.h>
#endif
#if defined(UNIX) || defined(WIN32)
#include <sys/types.h>
#endif

#include <limits.h>
#include <ctype.h>
#include <float.h>

#include <string.h>
#include <time.h>
#include <math.h>

#include "defs.h"
#include "param.h"
#include "structs.h"

#include "mm_file.h"
#include "best_stats.h"

#define XTERNAL
#include "uascii.h"
#include "upam.h"
#undef XTERNAL

#ifdef DEBUG
extern char ext_qtitle[];
#endif

extern void abort ();

#include "drop_func.h"	/* get init_work() */
/* drop_func.h includes dyn_string.h */

void revcomp(unsigned char *seq, int n, int *c_nt);
extern void qshuffle(unsigned char *aa0, int n0, int nm0, void *);
#ifdef DEBUG
unsigned long adler32(unsigned long, const unsigned char *, unsigned int);
#endif

void
s_annot_to_aa1a(long offset, int n1, struct annot_str *annot_p, unsigned char *ann_arr, char *tmp_line);

extern void add_annot_def(struct mngmsg *m_msp, char *line, int qa_flag);
int add_annot_char(unsigned char *ann_arr, char ctmp_label);

int get_annot(char *sname, struct mngmsg *, char *bline, long offset, int n1, 
	      struct annot_str **annot_p,int target, int debug);
int
get_annot_list(char *sname, struct mngmsg *m_msp, struct beststr **bestp_arr,
	       int nbest,int target, int debug);
void
print_sum(FILE *fd, struct db_str *qtt, struct db_str *ntt, int in_mem, long mem_use);
int
check_seq_range(unsigned char *aa1b, int n1, int nsq, char *str);
/* print timing information */
extern void ptime (FILE *, long);

/* this function consolidates code in comp_lib4.c for non-threaded, and in
   work_thr2.c (threads) and work_comp2.c (worker nodes)
*/

void
init_aa0(unsigned char **aa0, int n0, int nm0,
	 unsigned char **aa0s, unsigned char **aa1s,
	 int qframe, int qshuffle_flg, int max_tot,
	 struct pstruct *ppst, void **f_str, void **qf_str,
	 void *my_rand_state) {
  int id;

  /* note that aa[5,4,3,2] are never used, but are provided so that frame
     can range from 0 .. 5; likewise for f_str[5..2] */

  aa0[5] = aa0[4] = aa0[3] = aa0[2] = aa0[1] = aa0[0];

  /* zero out for SSE2/ALTIVEC -- make sure this is ALWAYS done */
  for (id=0; id < SEQ_PAD; id++) aa0[0][n0+id] = '\0';

  init_work (aa0[0], n0, ppst, &f_str[0]);
  f_str[5] = f_str[4] = f_str[3] = f_str[2] = f_str[1] = f_str[0];

  if (qframe == 2) {
    if ((aa0[1]=(unsigned char *)calloc((size_t)n0+2+SEQ_PAD,sizeof(unsigned char)))==NULL) {
      fprintf(stderr,"*** error [%s:%d] - cannot allocate aa01[%d]\n", __FILE__, __LINE__, n0);
    }
    *aa0[1]='\0';
    aa0[1]++;
    memcpy(aa0[1],aa0[0],n0+1);
    /* for ALTIVEC/SSE2, must pad with 16 NULL's */
    for (id=0; id<SEQ_PAD; id++) {aa0[1][n0+id]=0;}
    revcomp(aa0[1],n0,ppst->c_nt);
    init_work (aa0[1], n0, ppst, &f_str[1]);
  }

  if (qshuffle_flg) {
    if ((*aa0s=(unsigned char *)calloc(n0+2+SEQ_PAD,sizeof(char)))==NULL) {
      fprintf(stderr,"*** error [%s:%d] - cannot allocate aa0s[%d]\n",__FILE__, __LINE__, n0+2);
      exit(1);
    }
    **aa0s='\0';
    (*aa0s)++;
    memcpy(*aa0s,aa0[0],n0);
    qshuffle(*aa0s,n0,nm0, my_rand_state);
    /* for SSE2/ALTIVEC, must pad with 16 NULL's */
    for (id=0; id<SEQ_PAD; id++) {(*aa0s)[n0+id]=0;}
    init_work (*aa0s, n0, ppst, qf_str);
  }

  /* always allocate shuffle space */
  if((*aa1s=calloc(max_tot+1,sizeof(char))) == NULL) {
    fprintf(stderr,"*** error [%s:%d] - unable to allocate shuffled library sequence [%d]\n", __FILE__, __LINE__, max_tot);
    exit(1);
  }
  else {
    **aa1s=0;
    (*aa1s)++;
  }
}

/* because it is used to pre-allocate space, maxn has various
   constraints.  For "simple" comparisons, it is simply the length of
   the longest library sequence.  But for translated comparisons, it
   must be 3 or 6X the length of the query sequence.

   In addition, however, it can be reduced to make certain that
   sequences are read in smaller chunks.  And, maxn affect how large
   overlaps must be when sequences are read in chunks.
*/

int
reset_maxn(struct mngmsg *m_msp, int over_len, int maxn) {

  /* reduce maxn if requested */
  if (m_msp->ldb_info.maxn > 0 && m_msp->ldb_info.maxn < maxn) maxn = m_msp->ldb_info.maxn;

  if (m_msp->qdnaseq==m_msp->ldb_info.ldnaseq || m_msp->qdnaseq==SEQT_DNA ||
      m_msp->qdnaseq == SEQT_RNA) {/* !TFAST - either FASTA or FASTX */

    if (m_msp->n0 > m_msp->max_tot - m_msp->ldb_info.dupn) {
      fprintf(stderr,"*** error [%s:%d] -  query sequence is too long %d > %d - %d %s\n",
	      __FILE__, __LINE__,
	      m_msp->n0,
	      m_msp->max_tot, m_msp->ldb_info.dupn,
	      m_msp->sqnam);
      exit(1);
    }

    m_msp->ldb_info.l_overlap = over_len;
    m_msp->ldb_info.maxt3 = maxn-m_msp->ldb_info.l_overlap;
  }
  else {	/* is TFAST */
    if (m_msp->n0 > MAXTST) {
      fprintf(stderr,"*** error [%s:%d] -  query sequence is too long %d %s\n",
	      __FILE__, __LINE__, m_msp->n0,m_msp->sqnam);
      exit(1);
    }

    if (m_msp->n0*3 > maxn ) {	/* n0*3 for the three frames - this
				   will only happen if maxn has been
				   set low manually */

      if (m_msp->n0*4+2 < m_msp->max_tot) { /* m_msg0*3 + m_msg0 */
	fprintf(stderr,
		"*** error [%s:%d] - query sequence too long for library segment: %d - resetting to %d\n",
		__FILE__, __LINE__,
		maxn,m_msp->n0*3);
	maxn = m_msp->ldb_info.maxn = m_msp->n0*3;
      }
      else {
	fprintf(stderr,"*** error [%s:%d] -  query sequence too long for translated search: %d * 4 > %d %s\n",
		__FILE__, __LINE__, m_msp->n0,maxn, m_msp->sqnam);
	exit(1);
      }
    }

    /* set up some constants for overlaps */
    m_msp->ldb_info.l_overlap = 3*over_len;
    m_msp->ldb_info.maxt3 = maxn-m_msp->ldb_info.l_overlap-3;
    m_msp->ldb_info.maxt3 -= m_msp->ldb_info.maxt3%3;
    m_msp->ldb_info.maxt3++;

    maxn = maxn - 3; maxn -= maxn%3; maxn++;
  }
  return maxn;
}


int
scanseq(unsigned char *seq, int n, char *str) {
  int tot,i;
  char aaray[128];		/* this must be set > nsq */

  for (i=0; i<128; i++)  aaray[i]=0;
  for (i=0; i < (int)strlen(str); i++) aaray[qascii[str[i]]]=1;
  for (i=tot=0; i<n; i++) tot += aaray[seq[i]];
  return tot;
}

/* subs_env takes a string, possibly with ${ENV}, and looks up all the
   potential environment variables and substitutes them into the
   string */

void subs_env(char *dest, char *src, int dest_size) {
  char *last_src, *bp, *bp1;

  last_src = src;

  if ((bp = strchr(src,'$'))==NULL) {
    strncpy(dest, src, dest_size);
    dest[dest_size-1] = '\0';
  }
  else {
    *dest = '\0';
    while (strlen(dest) < dest_size-1 && bp != NULL ) {
      /* copy stuff before ${*/
      *bp = '\0';
      strncpy(dest, last_src, dest_size);
      *bp = '$';

      /* copy ENV */
      if (*(bp+1) != '{') {
	strncat(dest, "$", dest_size - strlen(dest) -1);
	dest[dest_size-1] = '\0';
	bp += 1;
      }
      else {	/* have  ${ENV} - put it in */
	if ((bp1 = strchr(bp+2,'}'))==NULL) {
	  fprintf(stderr, "*** error [%s:%d] - Unterminated ENV: %s\n",
		  __FILE__, __LINE__, src);
	  break;
	}
	else {
	  *bp1 = '\0';
	  if (getenv(bp+2)!=NULL) {
	    strncat(dest, getenv(bp+2), dest_size - strlen(dest) - 1);
	    dest[dest_size-1] = '\0';
	    *bp1 = '}';
	  }
	  bp = bp1+1;	/* bump bp even if getenv == NULL */
	}
      }
      last_src = bp;

      /* now get the next ${ENV} if present */
      bp = strchr(last_src,'$');
    }
    /* now copy the last stuff */
    strncat(dest, last_src, dest_size - strlen(dest) - 1);
    dest[dest_size-1]='\0';
  }
}


void
selectbest(struct beststr **bptr, int k, int n)	/* k is rank in array */
{
  int v, i, j, l, r;
  struct beststr *tmptr;

  l=0; r=n-1;

  while ( r > l ) {
    v = bptr[r]->rst.score[0];
    i = l-1;
    j = r;
    do {
      while (bptr[++i]->rst.score[0] > v) ;
      while (bptr[--j]->rst.score[0] < v) ;
      tmptr = bptr[i]; bptr[i]=bptr[j]; bptr[j]=tmptr;
    } while (j > i);
    bptr[j]=bptr[i]; bptr[i]=bptr[r]; bptr[r]=tmptr;
    if (i>=k) r = i-1;
    if (i<=k) l = i+1;
  }
}

void
selectbestz(struct beststr **bptr, int k, int n)	/* k is rank in array */
{
  int i, j, l, r;
  struct beststr *tmptr;
  double v;

  l=0; r=n-1;

  while ( r > l ) {
    v = bptr[r]->zscore;
    i = l-1;
    j = r;
    do {
      while (bptr[++i]->zscore > v) ;
      while (bptr[--j]->zscore < v) ;
      tmptr = bptr[i]; bptr[i]=bptr[j]; bptr[j]=tmptr;
    } while (j > i);
    bptr[j]=bptr[i]; bptr[i]=bptr[r]; bptr[r]=tmptr;
    if (i>=k) r = i-1;
    if (i<=k) l = i+1;
  }
}

/* improved shellsort with high-performance increments */
/*
shellsort(itemType a[], int l, int r)
{ int i, j, k, h; itemType v;
 int incs[16] = { 1391376, 463792, 198768, 86961, 33936,
		  13776, 4592, 1968, 861, 336,
		  112, 48, 21, 7, 3, 1 };
 for ( k = 0; k < 16; k++)
   for (h = incs[k], i = l+h; i <= r; i++) {
       v = a[i]; j = i;
       while (j > h && a[j-h] > v) {
         a[j] = a[j-h]; j -= h;
       }
       a[j] = v;
     }
}
*/

/* ?improved? version of sortbestz using optimal increments and fewer
   exchanges */
void sortbestz(struct beststr **bptr, int nbest)
{
  int gap, i, j, k;
  struct beststr *tmp;
  double v;
  int incs[14] = { 198768, 86961, 33936,
		   13776, 4592, 1968, 861, 336,
		   112, 48, 21, 7, 3, 1 };

  for ( k = 0; k < 14; k++) {
    gap = incs[k];
    for (i=gap; i < nbest; i++) {
      tmp = bptr[i];
      j = i;
      v = bptr[i]->zscore;
      while ( j >= gap && bptr[j-gap]->zscore < v) {
	bptr[j] = bptr[j - gap];
	j -= gap;
      }
      bptr[j] = tmp;
    }
  }
}


/* sort based on sequence index */
void sortbesti(struct beststr **bptr, int nbest)
{
  int gap, i, j, k;
  struct beststr *tmp;
  double v;
  int incs[12] = { 33936, 13776, 4592, 1968, 861, 336,
		   112, 48, 21, 7, 3, 1 };

  for ( k = 0; k < 12; k++) {
    gap = incs[k];
    for (i=gap; i < nbest; i++) {
      tmp = bptr[i];
      j = i;
      v = bptr[i]->seq->index;
      while ( j >= gap && bptr[j-gap]->seq->index < v) {
	bptr[j] = bptr[j - gap];
	j -= gap;
      }
      bptr[j] = tmp;
    }
  }
}

void
sortbeste(struct beststr **bptr, int nbest)
{
  int gap, i, j, k;
  struct beststr *tmp;
  double v;
  int incs[14] = { 198768, 86961, 33936,
		   13776, 4592, 1968, 861, 336,
		   112, 48, 21, 7, 3, 1 };

  for ( k = 0; k < 14; k++) {
    gap = incs[k];
    for (i=gap; i < nbest; i++) {
      j = i;
      tmp = bptr[i];
      v = tmp->rst.escore;
      while ( j >= gap && bptr[j-gap]->rst.escore > v) {
	bptr[j] = bptr[j - gap];
	j -= gap;
      }
      bptr[j] = tmp;
    }
  }

  /* sometimes there are many high scores with E()==0.0, sort
     those by z() score */

  j = 0;
  while (j < nbest && bptr[j]->rst.escore <= 2.0*DBL_MIN ) {j++;}
  if (j > 1) sortbestz(bptr,j);
}

extern char *prog_func;
extern char *verstr, *iprompt0, *refstr, *mp_verstr;
extern long tstart, tscan, tprev, tdone;	/* Timing */
#ifdef COMP_MLIB
extern long ttscan, ttdisp;
#endif
extern time_t tdstart, tddone;

/* ****************************************************************
   print command line arguments (argv_line)
   possibly HTML header
   !BLAST
     please cite
     version
   BLAST
     Reference version
**************************************************************** */
void
print_header1(FILE *fd, const char *argv_line,
	      const struct mngmsg *m_msp, const struct pstruct *ppst) {
  int i;

#ifdef PGM_DOC
  if (!(m_msp->markx & (MX_M8OUT+MX_MBLAST2))) fprintf(fd, "#%s\n",argv_line);
#endif

  if (m_msp->markx & MX_M11OUT) {
    fprintf(fd, "#:lav\n\nd {\n   \"%s\"\n}\n",argv_line+1);
  }

  if (m_msp->markx & MX_HTML) {
#ifdef HTML_HEAD
    fprintf(fd,"<html>\n<head>\n<title>%s Results</title>\n</head>\n<body>\n",prog_func);
#endif
    fprintf(fd,"<pre>\n");
  }

  if (m_msp->std_output) {
    fprintf(fd,"%s\n",iprompt0);
    if (refstr != NULL && refstr[0] != '\0') {
      fprintf(fd," version %s%s\nPlease cite:\n %s\n",verstr,mp_verstr,refstr);
    }
    else {
      fprintf(fd," version %s%s\n",verstr,mp_verstr);
    }
  }

  if (m_msp->markx & MX_MBLAST2) {
    if (refstr != NULL && refstr[0] != '\0') {
      fprintf(fd,"%s %s%s\n\nReference: %s\n\n", prog_func, verstr, mp_verstr, refstr);
    }
    else {
      fprintf(fd,"%s %s%s\n\n", prog_func, verstr, mp_verstr);
    }
  }

  fflush(fd);
}

/* ****************************************************************
   MX_HTML: <pre>
   Query:
     1>>>accession description # aa
   Annotation:
   Library:
**************************************************************** */
void
print_header2(FILE *fd, int qlib, char *info_qlabel, unsigned char **aa0,
	      const struct mngmsg *m_msp, const struct pstruct *ppst,
	      const char * info_lib_range_p) {
  int j;
  char tmp_str[MAX_STR];
  double db_tt;

  /* if (m_msp->markx & MX_HTML) fputs("<pre>\n",fd); */

  if (m_msp->std_output) {
    if (qlib==1) {
      fprintf(fd,"Query: %s\n", m_msp->tname);
    }

    if (m_msp->qdnaseq == SEQT_DNA || m_msp->qdnaseq == SEQT_RNA) {
      strncpy(tmp_str,(m_msp->qframe==1)? " (forward-only)" : "\0",sizeof(tmp_str));
      tmp_str[sizeof(tmp_str)-1]='\0';
    }
    else tmp_str[0]='\0';

    fprintf(fd,"%3d>>>%s%s\n", qlib,
	   m_msp->qtitle,
	   (m_msp->revcomp ? " (reverse complement)" : tmp_str));

    /* check for annotation */
    if (m_msp->ann_flg && m_msp->aa0a != NULL) {
      fprintf(fd,"Annotation: ");
      for (j=0; j<m_msp->n0; j++) {
	if (m_msp->aa0a[j] && m_msp->ann_arr[m_msp->aa0a[j]] != ' ' ) {
	  fprintf(fd,"|%ld:%c%c",
		  j+m_msp->q_off,m_msp->ann_arr[m_msp->aa0a[j]],ppst->sq[aa0[0][j]]);
	}
      }
    fprintf(fd,"\n");
    }

    fprintf(fd,"Library: %s%s\n", m_msp->ltitle,info_lib_range_p);

    if (m_msp->db.carry==0) {
      fprintf(fd, "  %7ld residues in %5ld sequences\n", m_msp->db.length, m_msp->db.entries);
    }
    else {
      db_tt = (double)m_msp->db.carry*(double)LONG_MAX + (double)m_msp->db.length;
      fprintf(fd, "  %.0f residues in %5ld library sequences\n", db_tt, m_msp->db.entries);
    }

  }
  else {
    if ((m_msp->markx & (MX_M8OUT + MX_M8COMMENT)) == (MX_M8OUT+MX_M8COMMENT)) {
      fprintf(fd,"# %s %s%s\n",prog_func,verstr,mp_verstr);
      fprintf(fd,"# Query: %s\n",m_msp->qtitle);
      fprintf(fd,"# Database: %s\n",m_msp->ltitle);
    }
  }
  if (m_msp->markx & MX_HTML) fputs("</pre>\n",fd);
  fflush(fd);
}

/* **************************************************************** */
/*   before showbest                                                */
/* **************************************************************** */
void print_header3(FILE *fd, int qlib, struct mngmsg *m_msp, struct pstruct *ppst) {

    if (m_msp->markx & MX_MBLAST2) {
      if (qlib == 1) {
	fprintf(fd, "\nDatabase: %s\n     %12ld sequences; %ld total letters\n\n\n",
		m_msp->ltitle, m_msp->db.entries, m_msp->db.length);
      }
      fprintf(fd, "\nQuery= %s\nLength=%d\n", m_msp->qtitle, m_msp->n0);
    }
}


/* **************************************************************** */
/* alignment tranistion                                             */
/* **************************************************************** */
void print_header4(FILE *fd, char *info_qlabel, char *argv_line, char *info_gstring3, char *info_hstring_p[2],
		   struct mngmsg *m_msp, struct pstruct *ppst) {

	if (m_msp->std_output && (m_msp->markx & (MX_AMAP+ MX_HTML + MX_M9SUMM)) && !(m_msp->markx & MX_M10FORM)) {
	  fprintf(fd,"\n>>>%s%s, %d %s vs %s library\n",
		  info_qlabel,(m_msp->revcomp ? "_rev":"\0"), m_msp->n0,
		  m_msp->sqnam,m_msp->lname);
	}

	if (m_msp->markx & MX_M10FORM) {
	  fprintf(fd,"\n>>>%s%s, %d %s vs %s library\n",
		  info_qlabel,(m_msp->revcomp ? "-":"\0"), m_msp->n0, m_msp->sqnam,
		  m_msp->lname);
	  fprintf(fd,"; pg_name: %s\n",m_msp->pgm_name);
	  fprintf(fd,"; pg_ver: %s%s\n",verstr,mp_verstr);
	  fprintf(fd,"; pg_argv: %s",argv_line);
	  fputs(info_gstring3,fd);
	  fputs(info_hstring_p[0],fd);
	  fputs(info_hstring_p[1],fd);
	}
}

void print_header4a(FILE *outfd, struct mngmsg *m_msp) {
  if (!(m_msp->markx & MX_M8OUT) && (m_msp->markx & (MX_M10FORM+MX_M9SUMM)) && m_msp->show_code != SHOW_CODE_ID) {
    fprintf(outfd,">>><<<\n");
  }
}

void print_header5(FILE *fd, int qlib, struct db_str *qtt,
		   struct mngmsg *m_msp, struct pstruct *ppst,
		   int in_mem, long tot_memK) {

  /* for MX_MBLAST2, show some statistics results */
  if (m_msp->markx & MX_MBLAST2) {
    fprintf(fd,"\n\nLambda      K     H\n");
    fprintf(fd," %6.3f  %6.3f  %6.3f\n\n",ppst->pLambda,ppst->pK,ppst->pH);
    fprintf(fd,"\nGapped\nLambda\n");
    fprintf(fd," %6.3f  %6.3f  %6.3f\n",ppst->pLambda,ppst->pK,ppst->pH);
    fprintf(fd,"\nEffective search space used: %ld\n\n",m_msp->db.entries);
  }

  if (m_msp->markx & MX_M8COMMENT) {
    fprintf(fd, "# %s processed %d queries\n",prog_func,qlib);
  }

  if ( !((m_msp->markx & MX_M8OUT) || (m_msp->markx & MX_HTML))
       && (m_msp->markx & (MX_M10FORM+MX_M9SUMM))) {
    fprintf(fd,">>>///\n");
  }

  if ( m_msp->markx & MX_HTML) fputs("<pre>",fd);
  if (m_msp->std_output) {
    print_sum(fd, qtt, &m_msp->db, in_mem, tot_memK);}
  if ( m_msp->markx & MX_HTML) fputs("</pre>\n",fd);
#ifdef HTML_HEAD
  if (m_msp->markx & MX_HTML) fprintf(fd,"</body>\n</html>\n");
#endif

  if (m_msp->markx & MX_MBLAST2) {
      fprintf(fd,"\n  Database: %s\n",m_msp->ltitle);
      fprintf(fd,"  Number of letters in database: %ld\n",m_msp->db.length);
      fprintf(fd,"  Number of sequences in database: %ld\n",m_msp->db.entries);
      fprintf(fd,"\n\n\nMatrix: %s\n",ppst->pam_name);
      fprintf(fd,"Gap Penalties: Existence: %d, Extension: %d\n",ppst->gdelval, ppst->ggapval);
  }
}

void
print_annot_header(FILE *fd, struct mngmsg *m_msp) {
  int i;

  if (m_msp->ann_arr_def[1]) {
    if (m_msp->markx & MX_HTML) {fprintf(fd,"<pre>");}
    fprintf(fd, "Annotation symbols:\n");
    for (i=1; m_msp->ann_arr[i]; i++) {
      if (m_msp->ann_arr_def[i]) {
	fprintf(fd, " %c : %s\n",m_msp->ann_arr[i], m_msp->ann_arr_def[i]);
      }
    }
    if (m_msp->markx & MX_HTML) {fputs("</pre><hr />\n",fd);}
  }
}

extern int fa_max_workers;

void
print_sum(FILE *fd, struct db_str *qtt, struct db_str *ntt, int in_mem, long tot_memK)
{
  double db_tt;
  char tstr1[26], tstr2[26];
  char memstr[256];

  strncpy(tstr1,ctime(&tdstart),sizeof(tstr1));
  strncpy(tstr2,ctime(&tddone),sizeof(tstr1));
  tstr1[24]=tstr2[24]='\0';

  /* Print timing to output file as well */

  fprintf(fd, "\n%ld residues in %ld query   sequences\n", qtt->length, qtt->entries);
  if (ntt->carry == 0)
    fprintf(fd, "%ld residues in %ld library sequences\n", ntt->length, ntt->entries);
  else {
    db_tt = (double)ntt->carry*(double)LONG_MAX + (double)ntt->length;
    fprintf(fd, "%.0f residues in %ld library sequences\n", db_tt, ntt->entries);
  }

  memstr[0]='\0';
  if (tot_memK && in_mem != 0) {
    sprintf(memstr," in memory [%ldG]",(tot_memK >> 20));
  }

#if !defined(COMP_THR) && !defined(PCOMPLIB)
  fprintf(fd," Scomplib [%s%s]\n start: %s done: %s\n",verstr,mp_verstr,tstr1,tstr2);
#endif
#if defined(COMP_THR)
  fprintf(fd," Tcomplib [%s%s] (%d proc%s)\n start: %s done: %s\n", verstr, mp_verstr,
	  fa_max_workers, memstr, tstr1,tstr2);
#endif
#if defined(PCOMPLIB)
  fprintf(fd," Pcomplib [%s%s] (%d proc%s)\n start: %s done: %s\n", verstr, mp_verstr,
	  fa_max_workers, memstr, tstr1,tstr2);
#endif
#ifndef COMP_MLIB
  fprintf(fd," Scan time: ");
  ptime(fd, tscan - tprev);
  fprintf (fd," Display time: ");
  ptime (fd, tdone - tscan);
#else
  fprintf(fd," Total Scan time: ");
  ptime(fd, ttscan);
  fprintf (fd," Total Display time: ");
  ptime (fd, ttdisp);
#endif
  fprintf (fd,"\n");
  fprintf (fd, "\nFunction used was %s [%s%s]\n", prog_func,verstr,mp_verstr);
}

extern double zs_to_Ec(double zs, long entries);
extern double zs_to_bit(double zs, int n0, int n1);
extern double zs_to_p(double zs);

#include "aln_structs.h"

void
prhist(FILE *fd, const struct mngmsg *m_msp,
       struct pstruct *ppst,
       struct hist_str hist,
       int nstats, int sstats,
       struct db_str ntt,
       char *stat_info2,
       char *lib_range,
       char **info_gstring2,
       char **info_hstring,
       long tscan)
{
  int i,j,hl,hll, el, ell, ev;
  char hline[80], pch, *bp;
  int mh1, mht;
  int maxval, maxvalt, dotsiz, ddotsiz,doinset;
  double cur_e, prev_e, f_int;
  double max_dev, x_tmp;
  double db_tt;
  int n_chi_sq, cum_hl=0, max_i=0, max_dev_i;
  double zs10_off;


  if (m_msp->markx & MX_HTML) fputs("<pre>\n",fd);
  else {fprintf(fd,"\n");}

  if (ppst->zsflag_f < 0) {
    if (!m_msp->nohist) {
      fprintf(fd, "  %7ld residues in %5ld sequences", ntt.length,ntt.entries);
      fprintf(fd, "%s\n",lib_range);
    }
    fprintf(fd,"Algorithm: %s\nParameters: %s\n",info_gstring2[0],info_gstring2[1]);
    return;
  }

  if (nstats > 20) {
    zs10_off = ppst->zs_off * 10.0;

    max_dev = 0.0;
    mh1 = hist.maxh-1;			/* max value for histogram */
    mht = (3*hist.maxh-3)/4 - 1;	/* x-coordinate for expansion */
    n_chi_sq = 0;

    if (!m_msp->nohist && mh1 > 0) {
      for (i=0,maxval=0,maxvalt=0; i<hist.maxh; i++) {
	if (hist.hist_a[i] > maxval) maxval = hist.hist_a[i];
	if (i >= mht &&  hist.hist_a[i]>maxvalt) maxvalt = hist.hist_a[i];
      }
      cum_hl = -hist.hist_a[0];
      dotsiz = (maxval-1)/60+1;
      ddotsiz = (maxvalt-1)/50+1;
      doinset = (ddotsiz < dotsiz && dotsiz > 2);

      if (ppst->zsflag_f>=0)
	fprintf(fd,"       opt      E()\n");
      else
	fprintf(fd,"     opt\n");

      prev_e =  zs_to_Ec((double)(hist.min_hist-hist.histint/2)-zs10_off,hist.entries);
      for (i=0; i<=mh1; i++) {
	pch = (i==mh1) ? '>' : ' ';
	pch = (i==0) ? '<' : pch;
	hll = hl = hist.hist_a[i];
	if (ppst->zsflag_f>=0) {
	  cum_hl += hl;
	  f_int = (double)(i*hist.histint+hist.min_hist)+(double)hist.histint/2.0;
	  cur_e = zs_to_Ec(f_int-zs10_off,hist.entries);
	  ev = el = ell = (int)(cur_e - prev_e + 0.5);
	  if (hl > 0  && i > 5 && i < (90-hist.min_hist)/hist.histint) {
	    x_tmp  = fabs(cum_hl - cur_e);
	    if ( x_tmp > max_dev) {
	      max_dev = x_tmp;
	      max_i = i;
	    }
	    n_chi_sq++;
	  }
	  if ((el=(el+dotsiz-1)/dotsiz) > 60) el = 60;
	  if ((ell=(ell+ddotsiz-1)/ddotsiz) > 40) ell = 40;
	  fprintf(fd,"%c%3d %5d %5d:",
		  pch,(i<mh1)?(i)*hist.histint+hist.min_hist :
		  mh1*hist.histint+hist.min_hist,hl,ev);
	}
	else fprintf(fd,"%c%3d %5d :",
		     pch,(i<mh1)?(i)*hist.histint+hist.min_hist :
		     mh1*hist.histint+hist.min_hist,hl);

	if ((hl=(hl+dotsiz-1)/dotsiz) > 60) hl = 60;
	if ((hll=(hll+ddotsiz-1)/ddotsiz) > 40) hll = 40;
	for (j=0; j<hl; j++) hline[j]='=';
	if (ppst->zsflag_f>=0) {
	  if (el <= hl ) {
	    if (el > 0) hline[el-1]='*';
	    hline[hl]='\0';
	  }
	  else {
	    for (j = hl; j < el; j++) hline[j]=' ';
	    hline[el-1]='*';
	    hline[hl=el]='\0';
	  }
	}
	else hline[hl] = 0;
	if (i==1) {
	  for (j=hl; j<10; j++) hline[j]=' ';
	  sprintf(&hline[10]," one = represents %d library sequences",dotsiz);
	}
	if (doinset && i == mht-2) {
	  for (j = hl; j < 10; j++) hline[j]=' ';
	  sprintf(&hline[10]," inset = represents %d library sequences",ddotsiz);
	}
	if (i >= mht&& doinset ) {
	  for (j = hl; j < 10; j++) hline[j]=' ';
	  hline[10]=':';
	  for (j = 11; j<11+hll; j++) hline[j]='=';
	  hline[11+hll]='\0';
	  if (ppst->zsflag_f>=0) {
	    if (ell <= hll) hline[10+ell]='*';
	    else {
	      for (j = 11+hll; j < 10+ell; j++) hline[j]=' ';
	      hline[10+ell] = '*';
	      hline[11+ell] = '\0';
	    }
	  }
	}

	fprintf(fd,"%s\n",hline);
	prev_e = cur_e;
      }
    }
    max_dev_i = max_i*hist.histint+hist.min_hist;
  }
  else {
    max_dev = 0.0;
    n_chi_sq = 0;
    max_i = 0;
    max_dev_i = 0;
  }

  if (ppst->zsflag_f >=0 ) {
    if (!m_msp->nohist) {
      if (ntt.carry==0) {
	fprintf(fd, "  %7ld residues in %5ld sequences", ntt.length, ntt.entries);
      }
      else {
	db_tt = (double)ntt.carry*(double)LONG_MAX + (double)ntt.length;
	fprintf(fd, "  %.0f residues in %5ld library sequences", db_tt, ntt.entries);
      }
      fprintf(fd, "%s\n",lib_range);
    }
    fprintf(fd,"Statistics: %s\n",hist.stat_info);
    if (stat_info2) {
      fprintf(fd," Statistics E2: %s\n",stat_info2);
    }

#ifdef SAMP_STATS
    fprintf(fd," statistics sampled from %ld (%d) to %ld sequences\n",
	    (hist.entries > nstats ? nstats : hist.entries),sstats, hist.entries);
#else
    fprintf(fd," statistics extrapolated from %ld to %ld sequences\n",
	    (hist.entries > nstats ? nstats : hist.entries),hist.entries);
#endif

    if (!m_msp->nohist && cum_hl > 0) {
      fprintf(fd," Kolmogorov-Smirnov  statistic: %6.4f (N=%d) at %3d\n",
	      max_dev/(float)cum_hl, n_chi_sq,max_dev_i);
    }
    if (m_msp->markx & MX_M10FORM) {
      while ((bp=strchr(hist.stat_info,'\n'))!=NULL) *bp=' ';
      if (cum_hl <= 0) cum_hl = -1;
      sprintf(info_hstring[0],"; mp_extrap: %d %ld\n; mp_stats: %s\n; mp_KS: %6.4f (N=%d) at %3d\n",
	      MAX_STATS,hist.entries,hist.stat_info,max_dev/(float)cum_hl,
	      n_chi_sq,max_dev_i);
    }
  }

  if (m_msp->markx & MX_M10FORM) {
    if ((bp = strchr(info_gstring2[1],'\n'))!=NULL) *bp = ' ';
    sprintf(info_hstring[1],"; mp_Algorithm: %s\n; mp_Parameters: %s\n",info_gstring2[0],info_gstring2[1]);
    if (bp != NULL ) *bp = '\n';
  }

  if (ppst->other_info != NULL) {
    fputs(ppst->other_info, fd);
  }

  fprintf(fd,"Algorithm: %s\nParameters: %s\n",info_gstring2[0],info_gstring2[1]);

  fprintf (fd," Scan time: ");
  ptime(fd,tscan);
  fprintf(fd,"\n");
  if (!m_msp->annot1_sname[0] &&  m_msp->markx & MX_HTML) {
    fputs("</pre>\n<hr />\n",fd);
  }

  fflush(fd);
}

extern char prog_name[], *verstr;

#ifdef PCOMPLIB
#include "mpi.h"
#endif

void s_abort (char *p,  char *p1)
{
  int i;

  fprintf (stderr, "\n***[%s] %s%s***\n", prog_name, p, p1);
#ifdef PCOMPLIB
  MPI_Abort(MPI_COMM_WORLD,1);
  MPI_Finalize();
#endif
  exit (1);
}

void w_abort (char *p, char *p1)
{
  fprintf (stderr, "\n***[%s] %s%s***\n\n", prog_name, p, p1);
  exit (1);
}

extern struct a_res_str *
build_ares_code(unsigned char *aa0, int n0,
		unsigned char *aa1, struct seq_record *seq,
		int frame, int *have_ares, int repeat_thresh,
		const struct mngmsg *m_msp, struct pstruct *ppst,
		void *f_str
		);

extern struct lmf_str *
re_openlib(struct lmf_str *, int outtty);

#define MAX_BLINE 2048
#define RANLIB (m_fptr->ranlib)

extern int
re_getlib(unsigned char *, struct annot_str **,
	  int, int, int, int, int, long *, long *,
	  struct lmf_str *m_fptr);

/*
   pre_load_best loads a set of sequences using re_getlib

   it should be used for getting sequences for shuffling, and for showbest() if m_msg->quiet

   it both opens the m_file_p buffer, gets the bline[] descriptions,
   and reads the actual sequences.  In reading the sequences, it
   should first allocate one large buffer so that individual buffers do not need to be freed.
*/

void
pre_load_best(unsigned char *aa1save, int maxn,
	      struct beststr **bbp_arr, int nbest,
	      struct mngmsg *m_msp, int debug)
{
  int i, n1, bl_len, tmp_bline_len, l_llen;
  int seq_buf_len;
  char bline[MAX_BLINE];
  unsigned char  *seq_buf_p;
  char *bline_buf_p;

  struct beststr *bbp;
  struct lmf_str *m_fptr;

  /*
     calculate how much room we need for sequences and blines
  */

  if (m_msp->pre_load_done) return;

  seq_buf_len = 1;
  for (i=0; i<nbest; i++) {
    /* we are not (currently) allocating more than n1+1, because alignment is not vectorized,
       if it were vectorized, we would need n+16
    */
#ifdef DEBUG
    if (bbp_arr[i]->n1 != bbp_arr[i]->seq->n1) {
      fprintf(stderr,"*** error [%s:%d] - n1 (%d) != seq->n1 (%d)\n",
	      __FILE__, __LINE__, bbp_arr[i]->n1, bbp_arr[i]->seq->n1);
    }
#endif

    if (bbp_arr[i]->seq->aa1b == NULL) {
      seq_buf_len += bbp_arr[i]->seq->n1 + 1;
    }
  }

  /* have required sequence space (seq_buf_len), allocate it */

  if ((m_msp->aa1save_buf_b=(unsigned char *)calloc(seq_buf_len, sizeof(char)))==NULL) {
    fprintf(stderr, "*** error [%s:%d] - cannot allocate space[%d] for sequence encoding\n",
	    __FILE__, __LINE__, seq_buf_len);
    exit(1);
  }
  else {
    seq_buf_p = m_msp->aa1save_buf_b+1;		/* ensure there is an initial '\0' */
  }

  /* adjust description line length */
  l_llen = m_msp->aln.llen;
  if ((m_msp->markx & MX_M9SUMM) && m_msp->show_code != SHOW_CODE_ID) {
    l_llen += 40;
    if (l_llen > 200) l_llen=200;
  }

  tmp_bline_len = sizeof(bline)-1;
  if (!(m_msp->markx & MX_M10FORM) && !m_msp->long_info) {tmp_bline_len = l_llen-5;}

  /* allocate more bline than we need for simplicity */
  if ((bline_buf_p=m_msp->bline_buf_b=(char *)calloc(nbest*tmp_bline_len, sizeof(char)))==NULL) {
    fprintf(stderr, "*** error [%s:%d] - cannot allocate space[%d] for bline descriptions\n",
	    __FILE__, __LINE__,  nbest*tmp_bline_len);
    exit(1);
  }

  for (i=0; i<nbest; i++) {
    bbp = bbp_arr[i];

    if ((m_fptr=re_openlib(bbp->mseq->m_file_p,!m_msp->quiet))==NULL) {
      fprintf(stderr,"*** error [%s:%d] - cannot re-open %s\n",
	      __FILE__, __LINE__, bbp->mseq->m_file_p->lb_name);
      exit(1);
    }
    RANLIB(bline,tmp_bline_len,bbp->mseq->lseek,bbp->mseq->libstr,m_fptr);
    bl_len = strlen(bline);
    bbp->mseq->bline = bline_buf_p;
    bbp->mseq->bline_max = m_msp->aln.llen;
    strncpy(bbp->mseq->bline, bline, bl_len);
    bline_buf_p += bl_len+1;

    /* make sure we get annotation if present, and sequence if necessary */
    if (bbp->seq->aa1b==NULL || (m_msp->ann_flg==1 && bbp->seq->annot_p==NULL)) {
      n1 = re_getlib(aa1save, (m_msp->ann_flg==1) ? &(bbp->seq->annot_p) : NULL,
		     maxn,m_msp->ldb_info.maxt3, m_msp->ldb_info.l_overlap,
		     bbp->mseq->cont,m_msp->ldb_info.term_code,
		     &bbp->seq->l_offset,&bbp->seq->l_off,bbp->mseq->m_file_p);
      if (n1 != bbp->seq->n1) {
	fprintf(stderr,"*** error [%s:%d] - n1[%d/%d] != n1[%d] from re_getlib() at %s [maxn:%d/maxt3:%d]\n",
		__FILE__, __LINE__,
		bbp->n1, bbp->seq->n1, n1, bbp->mseq->libstr, maxn, m_msp->ldb_info.maxt3);
      }

#ifdef DEBUG
      if (adler32(1L,aa1save,n1)!=bbp->adler32_crc) {
	fprintf(stderr,"*** error [%s:%d] - adler32_crc from re_getlib() at %d(%d):  %s\n",
		__FILE__, __LINE__,
		bbp->mseq->index,bbp->n1, bline);
      }
#endif

      /* if we don't have the sequence in the aa1b buffer, copy it from re_getlib */
      if (bbp->seq->aa1b == NULL)  {
	bbp->seq->aa1b = seq_buf_p;
	memcpy(bbp->seq->aa1b, aa1save, bbp->seq->n1+1);
	seq_buf_p += bbp->seq->n1+1;
      }
    }
  }

  /* here, we are getting query annots after all the bptr[]s have been processed */
  /* moved to comp_lib9.c */
  /*
  if (m_msp->annot0_sname[0]) {
    if (get_annot(m_msp->annot0_sname, m_msp, m_msp->qtitle, m_msp->q_offset+m_msp->q_off-1,m_msp->n0, &m_msp->annot_p, 0, debug) < 0) {
      fprintf(stderr,"*** error [%s:%d] - %s did not produce annotations\n",__FILE__, __LINE__, m_msp->annot0_sname);
      m_msp->annot0_sname[0] = '\0';
    }
    if (m_msp->annot_p && m_msp->annot_p->n_annot > 0) {
      m_msp->aa0a = m_msp->annot_p->aa1_ann;
    }
    if (!m_msp->ann_arr[0]) {m_msp->ann_arr[0] = ' '; m_msp->ann_arr[1] = '\0';}
  }
  */

  /* if we have an variant annotation script, execute it and capture the output */
  /* must do after bline is set */
  if (m_msp->annot1_sname[0]) {
    if (get_annot_list(m_msp->annot1_sname, m_msp, bbp_arr, nbest, 1, debug)< 0) {
      fprintf(stderr,"*** error [%s:%d] - %s did not produce annotations for %s\n",__FILE__, __LINE__, m_msp->annot1_sname,m_msp->qtitle);
      m_msp->annot1_sname[0] = '\0';
    };
    if (!m_msp->ann_arr[0]) {m_msp->ann_arr[0] = ' '; m_msp->ann_arr[1] = '\0';}
  }

  m_msp->pre_load_done = 1;
}

/*  merge_ares_chains()

    seeks to merge two ares chains, producing a single chain that is
    sorted by sw_score.

    Strategy -- choose the chain with the highest score, and go down
    it until the head of the other chain has higher score, then link
    the other chain to the main chain, breaking the first, and
    continue the process.

    The two pointers, max_next and alt_next, keep track of the best
    and the alternate chain
 */


#undef SHOW_MERGE_CHAIN

struct a_res_str *
merge_ares_chains(struct a_res_str *cur_ares,
		  struct a_res_str *tmp_ares,
		  int score_ix,
		  const char *msg)
{
  struct a_res_str *max_next, *max_ares, *alt_ares, *prev_next;

  if (!tmp_ares) return cur_ares;

#ifdef SHOW_MERGE_CHAIN
  fprintf(stderr,"cur_ares->");
  for (max_next = cur_ares; max_next; max_next = max_next->next) {
    fprintf(stderr,"%d->",max_next->rst.score[score_ix]);
  }

  fprintf(stderr,"||\n");
  fprintf(stderr,"tmp_ares->");
  for (max_next = tmp_ares; max_next; max_next = max_next->next) {
    fprintf(stderr,"%d->",max_next->rst.score[score_ix]);
  }
  fprintf(stderr,"||\n");
#endif

  /* start with the maximum score */

  if (cur_ares->rst.score[score_ix] >= tmp_ares->rst.score[score_ix]) {
    max_ares = max_next = prev_next = cur_ares;
    alt_ares = tmp_ares;
  }
  else {
    max_ares = max_next = prev_next = tmp_ares;
    alt_ares = cur_ares;
  }

  while (max_next && alt_ares) {
    /* this is guaranteed true for the first iteration */
    while (max_next && max_next->rst.score[score_ix] >= alt_ares->rst.score[score_ix]) {
      prev_next = max_next;
      max_next = max_next->next;
    }
    if (max_next==NULL) break;
    else {	/* max_next->rst.score[score_ix] no longer greater, switch
		   pointers */
      prev_next->next = alt_ares;
      alt_ares = max_next;
      max_next = prev_next->next;
    }
  }

  /* we quit whenever max_next or alt_ares == NULL; if
     (max_next==NULL), then continue adding the rest of alt_ares */

  if (max_next==NULL) {
    prev_next->next = alt_ares;
  }


#ifdef SHOW_MERGE_CHAIN
  fprintf(stderr,"[%s] merge_ares->",msg);
  for (max_next = max_ares; max_next; max_next = max_next->next) {
    fprintf(stderr,"%d->",max_next->rst.score[score_ix]);
  }
  fprintf(stderr,"||\n\n");
#endif

  return max_ares;
}

/* copies from from to to shuffling */

extern int my_nrand(int, void *);

void
shuffle(unsigned char *from, unsigned char *to, int n, void *rand_state)
{
  int i,j; unsigned char tmp;

  if (from != to) memcpy((void *)to,(void *)from,n);

  for (i=n; i>0; i--) {
    j = my_nrand(i, rand_state);
    tmp = to[j];
    to[j] = to[i-1];
    to[i-1] = tmp;
  }
  to[n] = 0;
}

/* shuffles DNA sequences as codons */
void
shuffle3(unsigned char *from, unsigned char *to, int n, void *rand_state)
{
  int i, j, i3,j3; unsigned char tmp;
  int n3;

  if (from != to) memcpy((void *)to,(void *)from,n);

  n3 = n/3;

  for (i=n3; i>0; i--) {
    j = my_nrand(i, rand_state);
    i3 = i*3;
    j3 = j*3;
    tmp = to[j3];
    to[j3] = to[i3-1];
    to[i3-1] = tmp;
    tmp = to[j3+1];
    to[j3+1] = to[i3];
    to[i3] = tmp;
    tmp = to[j3+2];
    to[j3+2] = to[i3+1];
    to[i3+1] = tmp;
  }
  to[n] = 0;
}

/* "shuffles" by reversing the sequence */
void
rshuffle(unsigned char *from, unsigned char *to, int n)
{
  unsigned char *ptr = from + n;

  while (n-- > 0) {
    *to++ = *ptr--;
  }
  *to = '\0';
}

static int ieven = 0;
/* copies from from to from shuffling, ieven changed for threads */
void
wshuffle(unsigned char *from, unsigned char *to, int n, int wsiz, void *rand_state)
{
  int i,j, k, mm;
  unsigned char tmp, *top;

  memcpy((void *)to,(void *)from,n);

  mm = n%wsiz;

  if (ieven) {
    for (k=0; k<(n-wsiz); k += wsiz) {
      top = &to[k];
      for (i=wsiz; i>0; i--) {
	j = my_nrand(i, rand_state);
	tmp = top[j];
	top[j] = top[i-1];
	top[i-1] = tmp;
      }
    }
    top = &to[n-mm];
    for (i=mm; i>0; i--) {
      j = my_nrand(i, rand_state);
      tmp = top[j];
      top[j] = top[i-1];
      top[i-1] = tmp;
    }
    ieven = 0;
  }
  else {
    for (k=n; k>=wsiz; k -= wsiz) {
      top = &to[k-wsiz];
      for (i=wsiz; i>0; i--) {
	j = my_nrand(i, rand_state);
	tmp = top[j];
	top[j] = top[i-1];
	top[i-1] = tmp;
      }
    }
    top = &to[0];
    for (i=mm; i>0; i--) {
      j = my_nrand(i, rand_state);
      tmp = top[j];
      top[j] = top[i-1];
      top[i-1] = tmp;
    }
    ieven = 1;
  }
  to[n] = 0;
}

int
sfn_cmp(int *q, int *s)
{
  if (*q == *s) return *q;
  while (*q && *s) {
    if (*q == *s) return *q;
    else if (*q < *s) q++;
    else if (*q > *s) s++;
  }
  return 0;
}

#ifndef MPI_SRC
void
revcomp(unsigned char *seq, int n, int *c_nt)
{
  unsigned char tmp;
  int i, ni;

  for (i=0, ni = n-1; i< n/2; i++,ni--) {
    tmp = c_nt[seq[i]];
    seq[i] = c_nt[seq[ni]];
    seq[ni] = tmp;
  }
  if ((n%2)==1) {
    i = n/2;
    seq[i] = c_nt[seq[i]];
  }
  seq[n]=0;
}
#endif

/* check to see whether this score (or a shuff score) should
   be included in statistics */
int samp_stats_idx (int *pre_nstats, int nstats, void *rand_state) {
  int jstats = -1;

  /* this code works when every score can be used for statistics
     estimates, but fails for fasta/[t]fast[xy] where only a fraction
     of scores are used */

  if (*pre_nstats < MAX_STATS) {
    jstats = (*pre_nstats)++;
  }

  /* here, the problem is that while we may have pre_nstats
     possible samplings, in some cases (-M subsets, fasta,
     [t]fast[xy] we don't have MAX_STATS samples yet.  Until we
     have MAX_STATS, we want more.  But the stats_idx strategy
     means that there may be additional samples in the buffers
     that are not reflected in nstats.
  */

  else {
#ifdef SAMP_STATS_LESS
    /* now we have MAX_STATS samples
       we want to sample 1/2 of 60K - 120K, 1/3 of 120K - 180K, etc */
    /* check every 15K to see if we have gone past the next threshold */

    /* pre_nstats cannot be incremented before the % to ensure
       that stats_inc is incremented exactly at 60000, 120000, etc.
       use ">=" in case increment comes later
       tests suggest the first 60K are sampled about 25% more
       than the rest
    */
    if (nstats < MAX_STATS) {
      jstats = MAX_STATS - my_nrand(MAX_STATS - nstats, rand_state)-1;
    }
    else if (((*pre_nstats)++ % (MAX_STATS/4)) == 0 &&
	     *pre_nstats >= stats_inc * MAX_STATS) {
      stats_inc = (*pre_nstats / MAX_STATS) + 1;
    }
    if ((*pre_nstats % stats_inc) == 0) {
      jstats = my_nrand(MAX_STATS, rand_state);
    }
#else
    /* this sampling strategy calls my_nrand() for every
       sequence > 60K, but provides a very uniform sampling */
    jstats = my_nrand(++(*pre_nstats), rand_state);
    if (jstats >= MAX_STATS) { jstats = -1;}
#endif
  }
  return jstats;
}

/* **************************************************************** */
/* build_link_data -- produces fasta file from m_msp->
   (1) generate a temporary file name
   (2) write out accessions \t expects to the temporary file
   (3) run script against temporary file, producing fasta_file_expansion_file
   (4) return fasta expansion filename for standard fasta openlib().

   returns: the expansion library file name
            **link_link_file_p is the name of the file with the data
              that will be removed.
*/
/* **************************************************************** */
char *
build_link_data(char **link_lib_file_p,
		struct mngmsg *m_msp, struct beststr **bestp_arr,
		int debug) {
  int i, status;
  char tmp_line[MAX_SSTR];
  char link_acc_file[MAX_STR];
  int link_acc_fd;
  char *link_lib_file;
  char *link_lib_str;
  char link_script[MAX_LSTR];
  int link_lib_type;
  char *bp, *link_bp;
  FILE *link_fd=NULL;		/* file for link accessions */

#ifndef UNIX
  return NULL;
#else
  /* get two tmpfiles, one for accessions, one for library */
  link_acc_file[0] = '\0';

  if ((link_lib_file=(char *)calloc(MAX_STR,sizeof(char)))==NULL) {
    fprintf(stderr,"*** error [%s:%d] - [build_link_data] Cannot allocate link_lib_file",
	    __FILE__, __LINE__);
  }
  link_lib_file[0] = '\0';

  if ((bp=getenv("TMP_DIR"))!=NULL) {
    strncpy(link_acc_file,bp,sizeof(link_acc_file));
    link_acc_file[sizeof(link_acc_file)-1] = '\0';
    SAFE_STRNCAT(link_acc_file,"/",sizeof(link_acc_file));
  }

  SAFE_STRNCAT(link_acc_file,"link_acc_XXXXXX",sizeof(link_acc_file));
  link_acc_fd = mkstemp(link_acc_file);
  strncpy(link_lib_file,link_acc_file,MAX_STR);
  link_acc_file[sizeof(link_acc_file)-1] = '\0';
  SAFE_STRNCAT(link_lib_file,".lib",MAX_STR);

  /* write out accessions to link_acc_file */
  if ((link_fd =fdopen(link_acc_fd,"w"))==NULL) {
    fprintf(stderr,"*** error [%s:%d] - Cannot open link_acc_file: %s\n",
	    __FILE__, __LINE__, link_acc_file);
    goto no_links;
  }

  for (i=0; i<m_msp->nskip + m_msp->nshow; i++) {
    if ((bp=strchr(bestp_arr[i]->mseq->bline,' '))!=NULL) {
      *bp = '\0';
    }
    fprintf(link_fd,"%s\t%.3g\n",bestp_arr[i]->mseq->bline,bestp_arr[i]->rst.escore);
    if (bp != NULL) *bp=' ';
  }
  fclose(link_fd);

  /* build link_script link_acc_file > link_lib_file */
  /* check for indirect */
  link_bp = &m_msp->link_lname[0];
  if (*link_bp == '!') {
    link_bp++;
  }
  if (*link_bp == '@') {
    link_bp++;
  }
  
  /* remove library type */
  if ((bp=strchr(link_bp,' '))!=NULL) {
    *bp = '\0';
    sscanf(bp+1,"%d",&link_lib_type);
  }
  else {
    link_lib_type = 0;
  }

  strncpy(link_script,link_bp,sizeof(link_script));
  link_script[sizeof(link_script)-1] = '\0';
  SAFE_STRNCAT(link_script," ",sizeof(link_script));
  SAFE_STRNCAT(link_script,link_acc_file,sizeof(link_script));
  SAFE_STRNCAT(link_script," >",sizeof(link_script));
  SAFE_STRNCAT(link_script,link_lib_file,sizeof(link_script));

  /* un-edit m_msp->link_lname */
  if (bp != NULL) *bp = ' ';

  /* run link_script link_acc_file > link_lib_file */
  status = system(link_script);
  if (!debug) {
#ifdef UNIX    
    unlink(link_acc_file);
#else
    _unlink(link_acc_file);
#endif
  }

  if (status == NO_FILE_EXIT) {	/* my specific return for no links */
    goto no_links;
  }

  if (status < 0 || status == 127) {
    fprintf(stderr,"*** error [%s:%d] - script: %s failed\n",
	    __FILE__, __LINE__,link_script);
    goto no_links;
  }

  if ((link_fd=fopen(link_lib_file,"r"))==NULL) {
    goto no_links;
  }
  else fclose(link_fd);

  if ((link_lib_str=(char *)calloc(MAX_STR,sizeof(char)))==NULL) {
    fprintf(stderr,"*** error [%s:%d] - [build_link_data] Cannot allocate link_lib_str",
	    __FILE__, __LINE__);
  }

  /* build the file string (possibly @link_lib_file libtype) */
  link_lib_str[0]='\0';
  if (m_msp->link_lname[0] == '@') {
    SAFE_STRNCAT(link_lib_str,"@",MAX_STR);
  }
  SAFE_STRNCAT(link_lib_str,link_lib_file,MAX_STR);
  if (link_lib_type > 0) {
    sprintf(tmp_line," %d",link_lib_type);
    SAFE_STRNCAT(link_lib_str,tmp_line,MAX_STR);
  }

  *link_lib_file_p = link_lib_file;
  return link_lib_str;

 no_links:
  free(link_lib_file);
  *link_lib_file_p = NULL;
  return NULL;
#endif
}

/* **************************************************************** */
/* build_lib_db -- produces fasta file from script
   (1) generate a temporary file name lib_db_file
   (2) run script producing data in lib_db_file

   returns: the expansion library file name
            **db_str_file_p is the name of the file with the data
              that will be removed.
*/
/* **************************************************************** */
char *
build_lib_db(char *script_file) {
  int i, status;
  char tmp_line[MAX_SSTR];
  char *lib_db_file, *lib_db_str;
  char lib_db_script[MAX_LSTR];
  int lib_db_indirect;
  int lib_db_type;
  int lib_db_str_len;
  char *bp, *lib_bp;

  if ((lib_db_file=(char *)calloc(MAX_STR,sizeof(char)))==NULL) {
    fprintf(stderr,"*** error [%s:%d] - [build_lib_db] Cannot allocate lib_db_file",
	    __FILE__, __LINE__);
    goto no_lib;
  }

  if ((bp=getenv("TMP_DIR"))!=NULL) {
    strncpy(lib_db_file,bp,MAX_STR);
    lib_db_file[sizeof(lib_db_file)-1] = '\0';
    SAFE_STRNCAT(lib_db_file,"/",sizeof(lib_db_file));
  }

  SAFE_STRNCAT(lib_db_file,"lib_db_XXXXXX",MAX_STR);
  mktemp(lib_db_file);
  lib_db_str_len = strlen(lib_db_file)+1;

  /* check for indirect */
  lib_bp = script_file;
  if (*lib_bp == '@') {
    lib_bp++;
    lib_db_str_len++;
  }
  /* remove library type */
  if ((bp=strchr(lib_bp,' '))!=NULL) {
    *bp = '\0';
    sscanf(bp+1,"%d",&lib_db_type);
    lib_db_str_len += (strlen(bp+1)+1);
  }
  else {
    lib_db_type = 0;
  }

  strncpy(lib_db_script,lib_bp,sizeof(lib_db_script));
  lib_db_script[sizeof(lib_db_script)-1] = '\0';
  SAFE_STRNCAT(lib_db_script," >",sizeof(lib_db_script));
  SAFE_STRNCAT(lib_db_script,lib_db_file,sizeof(lib_db_script));

  if (bp != NULL) *bp = ' ';

  /* run lib_db_script link_acc_file > lib_db_file */
  status = system(lib_db_script);

  if (status == NO_FILE_EXIT) {	/* my specific return for no links */
    goto no_lib;
  }

  if (status < 0 || status == 127) {
    fprintf(stderr,"*** error [%s:%d] - [build_lib_db] script: %s failed\n",
	    __FILE__, __LINE__, lib_db_script);
    goto no_lib;
  }

  /* build the file string (possibly @lib_db_str libtype) */
  if ((lib_db_str=calloc(lib_db_str_len+1,sizeof(char)))==NULL) {
    fprintf(stderr,"*** error [%s:%d] - [build_lib_db] cannot allocate lib_db_str[%d]\n",
	    __FILE__, __LINE__, lib_db_str_len+1);
    goto no_lib;
  }
  lib_db_str[0]='\0';
  if (*script_file == '@') {
    SAFE_STRNCAT(lib_db_str,"@",MAX_STR);
  }
  SAFE_STRNCAT(lib_db_str,lib_db_file,MAX_STR);
  if (lib_db_type > 0) {
    sprintf(tmp_line," %d",lib_db_type);
    SAFE_STRNCAT(lib_db_str,tmp_line,MAX_STR);
  }

  return lib_db_str;

 no_lib:
  return NULL;
}

/* used to temporarily allocate annotation array in next_annot_entry()*/
struct annot_mstr {
  int max_annot;
  struct annot_entry *tmp_arr_p;
};

/* init_tmp_annot_mstr(size) intializes the structure used to track annots  */
int
init_tmp_annot(struct annot_mstr *this, int size) {
  struct annot_entry *tmp_ann_astr;

  /* only reset if array is NULL */
  if (this->tmp_arr_p == NULL || this->max_annot <= 0) {
    this->max_annot = 32;
    if ((this->tmp_arr_p=(struct annot_entry *)calloc(this->max_annot, sizeof(struct annot_entry)))==NULL) {
      fprintf(stderr,"*** error [%s:%d] - cannot allocate annot_entry[%d]\n",
		__FILE__,__LINE__,this->max_annot);
	return 0;
    }
  }
  return 1;
}

int
update_tmp_annot(struct annot_mstr *this) {

  this->max_annot += (this->max_annot/2);
  if ((this->tmp_arr_p= (struct annot_entry *)realloc(this->tmp_arr_p, this->max_annot*sizeof(struct annot_entry)))==NULL) {
    fprintf(stderr,"[*** error [%s:%d] - cannot reallocate tmp_ann_astr[%d]\n",
	    __FILE__, __LINE__, this->max_annot);
    return 0;
  }
  return 1;
}

struct annot_str *
next_annot_entry(FILE *annot_fd, char *tmp_line, int n_tmp_line,
		 struct annot_str *annot_p,
		 struct annot_mstr *mtmp_annot_p,
		 struct mngmsg *m_msp, int target);

/* **************************************************************** */
/* get_annot_list -- produces fasta file from sname
   if sname[0]=='<', read the file directly, goto (4)
   if sname[0]=='!', run a script
   (1) generate a temporary file name
   (2) write out list of blines
   (3) run m_msp->annot1_sname[] script against temporary file, producing table of annotations

   (4) read in the annotations and merge them into beststr
   (5) return number of annotations
*/
/* **************************************************************** */

int
get_annot_list(char *sname, struct mngmsg *m_msp, struct beststr **bestp_arr, int nbest,
	       int target, int debug) {
  int i, status;
  long l_offset;
  char tmp_line[MAX_STR];
  char annot_bline_file[MAX_STR];
  int annot_bline_fd;
  char *annot_descr_file;
  char annot_script[MAX_LSTR];
  struct annot_str *annot_p;
  char *bp;
  int annot_seq_cnt;
  FILE *annot_fd=NULL;		/* file for annot accessions */
  struct annot_mstr mtmp_annot;

#ifndef UNIX
  return 0;
#else

  if (sname[0] == '!') {

    /* get two tmpfiles, one for bline, one for returned annotations
       (but it would make more sense to use popen() to get the
       annotations back
    */

    annot_bline_file[0] = '\0';

    if ((annot_descr_file=(char *)calloc(MAX_STR,sizeof(char)))==NULL) {
      fprintf(stderr,"*** error [%s:%d] - [get_annot_list] Cannot allocate annot_file",
	      __FILE__, __LINE__);
    }
    annot_descr_file[0] = '\0';

    if ((bp=getenv("TMP_DIR"))!=NULL) {
      strncpy(annot_bline_file,bp,sizeof(annot_bline_file));
      annot_bline_file[sizeof(annot_bline_file)-1] = '\0';
      SAFE_STRNCAT(annot_bline_file,"/",sizeof(annot_bline_file));
    }

    SAFE_STRNCAT(annot_bline_file,"annot_bline_XXXXXX",sizeof(annot_bline_file));
    annot_bline_fd = mkstemp(annot_bline_file);
    strncpy(annot_descr_file,annot_bline_file,MAX_STR);
    annot_bline_file[sizeof(annot_bline_file)-1] = '\0';
    SAFE_STRNCAT(annot_descr_file,".annot",MAX_STR);

    /* write out accessions to annot_bline_file */
    if ((annot_fd =fdopen(annot_bline_fd,"w"))==NULL) {
      fprintf(stderr,"*** error [%s:%d] - Cannot open annot_bline_file: %s\n",__FILE__, __LINE__, annot_bline_file);
      goto no_annots;
    }

    for (i=0; i<nbest; i++) {
      if (bestp_arr[i]->mseq->annot_req_flag) {	continue; }
      if ((strlen(bestp_arr[i]->mseq->bline) > DESCR_OFFSET) &&
	  (bp=strchr(bestp_arr[i]->mseq->bline+DESCR_OFFSET,' '))!=NULL) {*bp = '\0';}
      else {bp = NULL;}
      /* provide sequence length with offset, but only if offset is positive */
      l_offset = bestp_arr[i]->seq->l_offset+bestp_arr[i]->seq->l_off -1;
      if (l_offset < 0) { l_offset = 0;}
      fprintf(annot_fd,"%s\t%ld\n",bestp_arr[i]->mseq->bline,
	      l_offset + bestp_arr[i]->seq->n1);
      if (bp != NULL) *bp=' ';
      bestp_arr[i]->mseq->annot_req_flag = 1;
    }
    fclose(annot_fd);

    subs_env(annot_script, sname+1, sizeof(annot_script));
    annot_script[sizeof(annot_script)-1] = '\0';
    SAFE_STRNCAT(annot_script," ",sizeof(annot_script));
    SAFE_STRNCAT(annot_script,annot_bline_file,sizeof(annot_script));
    SAFE_STRNCAT(annot_script," >",sizeof(annot_script));
    SAFE_STRNCAT(annot_script,annot_descr_file,sizeof(annot_script));

    /* run annot_script annot_bline_file > annot_descr_file */
    status = system(annot_script);
    if (!debug) {
#ifdef UNIX
      unlink(annot_bline_file);
#else
      _unlink(annot_bline_file);
#endif
    }

    if (status == NO_FILE_EXIT) {	/* my specific return for no annots */
      goto no_annots;
    }

    if (status < 0 || status == 127) {
      fprintf(stderr,"*** error [%s:%d] - script: %s failed\n",
	      __FILE__, __LINE__, annot_script);
      goto no_annots;
    }
  }
  else if (sname[0] == '<') {
    annot_descr_file = sname+1;
  }
  else {
    fprintf(stderr,"*** error [%s:%d] - %s not script (!) or file (<)\n",__FILE__, __LINE__, sname);
    goto no_annots;
  }

  /* read annotation file */

  if ((annot_fd=fopen(annot_descr_file,"r"))==NULL) {
    goto no_annots;
  }

  /* be sure to ask for annotation once */
  for (i=0; i<nbest; i++) {
    bestp_arr[i]->mseq->annot_req_flag = 1;
  }
  /* we have some annotations */
  /* the annotation script MUST return the annotations ordered as in annot_descr_file,
     in "fasta" form:

     >bline_descr
     pos<tab>label<tab>value?<tab>comment (which is not read in this version)
     1 *
     11 V N
  */

  /* now read the annotation/variant file */

  /* read #comments, =annot_defs at beginning of file */
  tmp_line[0] = '#';
  while (tmp_line[0] == '#' || tmp_line[0] == '=') {
    if (tmp_line[0] == '=') add_annot_def(m_msp, tmp_line+1,1);
    if (fgets(tmp_line, sizeof(tmp_line), annot_fd)==NULL) {
      fprintf(stderr,"*** error [%s:%d] - premature annotation file end (%s)\n",
	      __FILE__,__LINE__, annot_descr_file);
      goto no_annots;
    }
  }

  /* set mtmp_annot to be initialized */
  mtmp_annot.tmp_arr_p = NULL;
  mtmp_annot.max_annot = 0;

  annot_seq_cnt = 0;

  /* now read in the annotations, but only the first time if asked multiple times */
  for (i=0; i<nbest; i++) {
    if (!bestp_arr[i]->mseq->annot_req_flag) {
      continue;
    }
    bestp_arr[i]->mseq->annot_req_flag = 0;

    if ((bp=strchr(tmp_line,'\n'))!=NULL) *bp = '\0';
    if ((bp=strchr(tmp_line,'\t'))!=NULL) *bp = '\0';
    if (tmp_line[0] != '>' || strncmp(&tmp_line[1], bestp_arr[i]->mseq->bline, strlen(&tmp_line[1])) != 0) {
      fprintf(stderr,"*** error [%s:%d] - %s description mismatch (%s:%s)\n",
	      __FILE__,__LINE__,annot_descr_file, tmp_line, bestp_arr[i]->mseq->bline);
      goto no_annots;
    }

    annot_p = next_annot_entry(annot_fd, tmp_line, sizeof(tmp_line), bestp_arr[i]->seq->annot_p, &mtmp_annot, m_msp, target);

    if (annot_p) {
      bestp_arr[i]->seq->annot_p = annot_p;
      s_annot_to_aa1a(bestp_arr[i]->seq->l_offset + bestp_arr[i]->seq->l_off - 1,
		      bestp_arr[i]->seq->n1, annot_p,m_msp->ann_arr, bestp_arr[i]->mseq->libstr);
      annot_seq_cnt++;
      mtmp_annot.tmp_arr_p = NULL;
    }
    else {
      if (bestp_arr[i]->seq->annot_p) {
	bestp_arr[i]->seq->annot_p->n_annot = 0;
      }
    }
  }

  if (mtmp_annot.tmp_arr_p) free(mtmp_annot.tmp_arr_p);

  fclose(annot_fd);
  if (sname[0]=='!') {
    if (!debug) {
#ifdef UNIX
      unlink(annot_descr_file);
#else
      _unlink(annot_descr_file);
#endif
    }
    free(annot_descr_file);
  }
  return annot_seq_cnt;

 no_annots:
  for (i=0; i<nbest; i++) {
    if (bestp_arr[i]->seq->annot_p) {
      if (bestp_arr[i]->seq->annot_p->n_annot > 0) {
	bestp_arr[i]->seq->annot_p->n_annot = 0;
      }
    }
  }
  if (sname[0] == '!') free(annot_descr_file);
  return -1;
#endif
}

void sort_annots(struct annot_entry **s_annot, int n_annot)
{
  int gap, i, j, k;
  struct annot_entry *tmp;
  int v;
  int incs[6] = { 112, 48, 21, 7, 3, 1 };

  for ( k = 0; k < 6; k++) {
    gap = incs[k];
    for (i=gap; i < n_annot; i++) {
      tmp = s_annot[i];
      j = i;
      v = s_annot[i]->pos;
      while ( j >= gap && s_annot[j-gap]->pos > v) {
	s_annot[j] = s_annot[j - gap];
	j -= gap;
      }
      s_annot[j] = tmp;
    }
  }
}

struct annot_str *
next_annot_entry(FILE *annot_fd, char *tmp_line, int n_tmp_line, struct annot_str *annot_p,
		 struct annot_mstr *mtmp_annot_p, struct mngmsg *m_msp, int target) {

  char ctmp_label, ctmp_value, tmp_comment[MAX_STR], annot_acc[MAX_STR];
  char *bp;
  int f_pos, f_end;
  int i_ann, l_doms, r_doms;
  int n_annot = 0;

  struct annot_entry *tmp_ann_entry_arr, **s_tmp_ann_entry_arr;

  SAFE_STRNCPY(annot_acc, tmp_line, sizeof(annot_acc));

  if (init_tmp_annot(mtmp_annot_p, 32)==0) return NULL;

  tmp_ann_entry_arr = mtmp_annot_p->tmp_arr_p;

  /* read through each annotation in file */
  while (fgets(tmp_line, n_tmp_line, annot_fd)!=NULL ) {
    if (tmp_line[0] == '>') goto next_bline;	/* start of new annotation */
    if (tmp_line[0] == '#') continue;		/* ignore comments */
    if (tmp_line[0] == '=') {	/* symbol definition */
      add_annot_def(m_msp, tmp_line+1,1);
      continue;
    }

    if (n_annot >= mtmp_annot_p->max_annot - 1) {
      /* try to expand annotation array */
      if (update_tmp_annot(mtmp_annot_p)==0) {
	return NULL;
      }
      tmp_ann_entry_arr = mtmp_annot_p->tmp_arr_p;
    }

    /* sscanf cannot give strings with blanks */
    /* sscanf(tmp_line,"%d %c %c %s", &f_pos, &ctmp_label, &ctmp_value, tmp_comment); */
    tmp_comment[0] = '\0';
    if ((bp=strchr(tmp_line,'\r')) || (bp=strchr(tmp_line,'\n'))) *bp='\0';	/* clean up EOL */
    if ((bp=strchr(tmp_line,'\t'))!=NULL) {	/* fields MUST be tab delimited */
      f_pos=atoi(tmp_line) - 1;	/* get first field -- f_pos, converted to 0-offset  */
      ctmp_label = bp[1];	/* get second field -- ctmp_label */
      if ((bp=strchr(bp+1,'\t'))!=NULL) {	/* next field could be f_end or ctmp_value */
	if (ctmp_label == '-') { f_end = atoi(bp+1) -1; ctmp_value = '\0';}
	else {ctmp_value = bp[1]; f_end = f_pos;} 	/* have variant, not coordinate */
	if ((bp=strchr(bp+1,'\t'))!=NULL) {		/* if last <tab>, get comment */
	  strncpy(tmp_comment,bp+1,sizeof(tmp_comment));
	}
      }
    }
    else {	/* no tab */
      continue;
    }

    tmp_ann_entry_arr[n_annot].pos = f_pos;
    tmp_ann_entry_arr[n_annot].end = f_end;
    tmp_ann_entry_arr[n_annot].label=ctmp_label;
    tmp_ann_entry_arr[n_annot].value=ctmp_value;
    tmp_ann_entry_arr[n_annot].comment = NULL;
    tmp_ann_entry_arr[n_annot].target = target;
    tmp_ann_entry_arr[n_annot].start = NULL;

    if (tmp_comment[0]) {
      if ((tmp_ann_entry_arr[n_annot].comment=(char *)calloc(strlen(tmp_comment)+1,sizeof(char)))!=NULL) {
	strncpy(tmp_ann_entry_arr[n_annot].comment,tmp_comment,strlen(tmp_comment));
      }
    }
    if (ctmp_label== 'V') {
      /* map the .value from ascii to encoded residue */
      /* this must be lascii, not qascii, because the script
	 describes the library, not the query, and for FASTX/TFASTX,
	 those are different */
      tmp_ann_entry_arr[n_annot].value = lascii[tmp_ann_entry_arr[n_annot].value];
    }
    else if (ctmp_label == '[') {
      i_ann = add_annot_char(m_msp->ann_arr, ctmp_label);
      if (i_ann > 0) {
	qascii[ctmp_label] = NANN + i_ann;
	m_msp->ann_arr_def[i_ann] = NULL;
      }
    }
    else if (ctmp_label == ']') {
      i_ann = add_annot_char(m_msp->ann_arr, ctmp_label);
      if (i_ann > 0) {
	qascii[ctmp_label] = NANN + i_ann;
	m_msp->ann_arr_def[i_ann] = NULL;
      }
    }
    else if (ctmp_label == '-') {	/* if ctmp_label == '-', have f_end, which must be added with ']' */
      /* make sure start/stop characters are in annotation alphabet */
      i_ann = add_annot_char(m_msp->ann_arr, '[');
      if (i_ann > 0) {
	qascii['['] = NANN + i_ann;
	m_msp->ann_arr_def[i_ann] = NULL;
      }

      tmp_ann_entry_arr[n_annot].label='[';
      n_annot++;

      if (n_annot >= mtmp_annot_p->max_annot - 1) {
	if (update_tmp_annot(mtmp_annot_p)==0) {
	  return NULL;
	}
	tmp_ann_entry_arr = mtmp_annot_p->tmp_arr_p;
      }

      /* only required for start - stop annotations; requires sort
	 later (which is not currently used */
      tmp_ann_entry_arr[n_annot].pos = f_end;
      tmp_ann_entry_arr[n_annot].end = f_end;
      tmp_ann_entry_arr[n_annot].label=']';
      tmp_ann_entry_arr[n_annot].value=ctmp_value;
      tmp_ann_entry_arr[n_annot].comment = NULL;
      tmp_ann_entry_arr[n_annot].target = target;
      tmp_ann_entry_arr[n_annot].start = &tmp_ann_entry_arr[n_annot-1];

      i_ann = add_annot_char(m_msp->ann_arr, ']');
      if (i_ann > 0) {
	qascii[']'] = NANN + i_ann;
	m_msp->ann_arr_def[i_ann] = NULL;
      }
    }
    else if ((i_ann = add_annot_char(m_msp->ann_arr, ctmp_label)) > 0) {
      m_msp->ann_arr_def[i_ann] = NULL;
      qascii[ctmp_label] = NANN + i_ann;
    }
    n_annot++;
  }

 next_bline:
  if (n_annot) {  /* if we have annotations, save them and set tmp_ann_entry_arr = NULL */
    tmp_ann_entry_arr = (struct annot_entry *)realloc(tmp_ann_entry_arr, (n_annot+1)*sizeof(struct annot_entry));

    if ((s_tmp_ann_entry_arr = (struct annot_entry **)calloc((n_annot+1),sizeof(struct annot_entry *)))==NULL) {
      fprintf(stderr,"*** error [%s:%d] - cannot alloc s_tmp_ann_entry_arr[%d]",
	      __FILE__,__LINE__, n_annot+1);
      exit(1);
    }

    /* pair every domain start/stop */
    /* (1) count number of domains/check for consistency */
    l_doms = r_doms = 0;
    for (i_ann=0; i_ann < n_annot; i_ann++) {
      if (tmp_ann_entry_arr[i_ann].label == '[') l_doms++;
      if (tmp_ann_entry_arr[i_ann].label == ']') r_doms++;
    }
    if (l_doms != r_doms) {
      fprintf(stderr,"*** error [%s:%d] - unpaired domains: %s %d != %d\n",
	      __FILE__,__LINE__, annot_acc, l_doms, r_doms);
#ifdef DEBUG      
      for (i_ann=0; i_ann < n_annot; i_ann++) {
	if (tmp_ann_entry_arr[i_ann].label == '[') 
	  fprintf(stderr, "%ld %c %s\n",tmp_ann_entry_arr[i_ann].pos,tmp_ann_entry_arr[i_ann].label,tmp_ann_entry_arr[i_ann].comment);
	if (tmp_ann_entry_arr[i_ann].label == ']')
	  fprintf(stderr, "%ld %c\n",tmp_ann_entry_arr[i_ann].pos,tmp_ann_entry_arr[i_ann].label);
      }
#endif
    }
    else if (l_doms > 0) {
      if ((tmp_domain_entry_arr = (struct annot_entry *)calloc((l_doms+1),sizeof(struct annot_entry)))==NULL) {
	fprintf(stderr,"*** error [%s:%d] - cannot alloc s_tmp_ann_entry_arr[%d]",
		__FILE__,__LINE__, l_doms+1);
      }
      else {
	l_doms = 0;
	for (i_ann=0; i_ann < n_annot+1; i_ann++) {
	  if (tmp_ann_entry_arr[i_ann].label == '[') {
	    tmp_domain_entry_arr[l_doms].pos = tmp_ann_entry_arr[i_ann].pos;
	    tmp_domain_entry_arr[l_doms].label = '-';
	    tmp_domain_entry_arr[l_doms].comment = tmp_ann_entry_arr[i_ann].comment;
	  }
	  else if (tmp_ann_entry_arr[i_ann].label == ']') {
	    tmp_domain_entry_arr[l_doms].end = tmp_ann_entry_arr[i_ann].pos;
	    l_doms++;
	  }
	}      
      }
    }

    for (i_ann=0; i_ann < n_annot+1; i_ann++) {
      s_tmp_ann_entry_arr[i_ann] = &tmp_ann_entry_arr[i_ann];
    }

    sort_annots(s_tmp_ann_entry_arr,n_annot);

    /* now allocate an annot_p if necessary, and link tmp_ann_entry_arr to it */
    if (annot_p || (annot_p = calloc(1,sizeof(struct annot_str)))!=NULL) {
      annot_p->annot_arr_p = tmp_ann_entry_arr;
      annot_p->s_annot_arr_p = s_tmp_ann_entry_arr;
      annot_p->n_annot = n_annot;
      annot_p->n_domains = l_doms;
      /* set to NULL to re-initialize */
    }
  }
  else {
    annot_p = NULL;
  }
  return annot_p;
}


/* **************************************************************** */
/* add_annot_char(ann_arr, ctmp_label) --

   (1) add annotation character to ann_arr if not present
   (2) return i_ann if added
*/
/* **************************************************************** */

int
add_annot_char(unsigned char *ann_arr, char ctmp_label) {
  int i_ann;

  if (ann_arr[0] == '\0') {
    ann_arr[0] = ' '; ann_arr[1] = '\0';
  }

  /* check to see if already there? */
  if (strchr((char *)ann_arr,ctmp_label)==NULL) {
    /* check for room for another character */
    if (strlen((char *)ann_arr) >= MAX_FN) {
      fprintf(stderr,"*** error [%s:%d] -  too many annotation characters: len(%s) + %c > %d\n",
	      __FILE__, __LINE__, ann_arr, ctmp_label, MAX_FN-1);
      return 0;
    }
    else {
      ann_arr[i_ann=strlen((char *)ann_arr)] = ctmp_label;      /* add the character */
      ann_arr[i_ann+1] = '\0';	      /* guarantee null termination */
      return i_ann;
    }
  }
  else {
    return 0;
  }
}

/* **************************************************************** */
/* get_annot -- produces fasta file from m_msp->sname script
   (modified 20-Sept-2012 to not use intermediate file)

   # (1) generate a temporary file name
   # (2) write out one bline (or portion that include accession)
   # (3) run sname[] script against temporary file, producing table of annotations
   (1) run script bline_id
   (4) read in the annotations and put them in struct annot_entry;
   (5) modify *annot_p to point to structure
   (6) return number of annotations
*/
/* **************************************************************** */

int
get_annot(char *sname, struct mngmsg *m_msp, char *bline, long offset, int n1, struct annot_str **annot_p,
	  int target, int debug) {

  char tmp_line[MAX_STR];
  FILE *annot_data_fd;
  char bline_descr[MAX_STR];
  char annot_data_file[MAX_LSTR];
  char annot_script[MAX_LSTR];
  long q_offset;

  char *bp;
  FILE *annot_fd=NULL;		/* file for annot accessions */
  struct annot_mstr mtmp_annot;

#ifndef UNIX
  return 0;
#else

  if (sname[0] == '!') {
    /* popen implementation */

    annot_data_file[0] = '\0';

    if (bline[0] == '>') {
      SAFE_STRNCPY(bline_descr, bline+1,sizeof(bline_descr));
    }
    else {
      SAFE_STRNCPY(bline_descr, bline,sizeof(bline_descr));
    }
    if ((strlen(bline_descr) > DESCR_OFFSET) && 
	(bp=strchr(bline_descr+DESCR_OFFSET,' '))!=NULL) {*bp = '\0';}
    else {bp = NULL;}

    q_offset = m_msp->q_offset + m_msp->q_off - 1;
    if (q_offset < 0) { q_offset = 0;}
    sprintf(annot_script,"%s \"%s\" %ld",sname+1, bline_descr,q_offset+m_msp->n0);
    annot_script[sizeof(annot_script)-1] = '\0';

    annot_fd = popen(annot_script,"r");
  }
  else if (sname[0] == '<') {
    SAFE_STRNCPY(annot_data_file,sname+1,sizeof(annot_data_file));
    annot_fd=fopen(annot_data_file,"r");
  }
  else {
    fprintf(stderr,"*** error [%s:%d] - %s not script (!) or file (<)\n",__FILE__, __LINE__, sname);
    goto no_annots;
  }

  if (!annot_fd) {
    goto no_annots;
  }
  else {	/* read the annotations into the array */

    /* read #comments, =annot_defs at beginning of file */
    tmp_line[0] = '#';
    while (tmp_line[0] == '#' || tmp_line[0] == '=') {
      if (tmp_line[0] == '=') add_annot_def(m_msp, tmp_line+1,1);
      if (fgets(tmp_line, sizeof(tmp_line), annot_fd)==NULL) {
	fprintf(stderr,"*** error [%s:%d] - premature annotation file end (%s)\n",
		__FILE__,__LINE__, annot_data_file);
	goto no_annots;
      }
    }

    /* set mtmp_annot to be initialized */
    mtmp_annot.tmp_arr_p = NULL;
    mtmp_annot.max_annot = 0;

    /* strlen(&tmp_line[1])-1 to remove '>' and beginning and '\n' at end */
    if (tmp_line[0] != '>') {
      fprintf(stderr,"*** error [%s:%d] - no %s description: [%s]\n",
	      __FILE__,__LINE__,annot_data_file, tmp_line);
      goto no_annots;
    }

    *annot_p = next_annot_entry(annot_fd, tmp_line, sizeof(tmp_line), *annot_p, &mtmp_annot, m_msp, target);

    if (sname[0] == '!') {
      pclose(annot_fd);
    }
    else {
      fclose(annot_fd);
    }

    /* now allocate an annot_p if necessary, and link tmp_ann_entry_arr to it */
    if (*annot_p) {
      s_annot_to_aa1a(offset, n1, (*annot_p),m_msp->ann_arr,"get_annot");
      return (*annot_p)->n_annot;
    }
    else {
      if (mtmp_annot.tmp_arr_p) free(mtmp_annot.tmp_arr_p);
      return 0;
    }
  }

 no_annots:
  return -1;
#endif
}

/* s_annot_to_aa1a -- takes an annot_entry[] and converts it to an *aa1_ann
 */
void
s_annot_to_aa1a(long offset, int n1, struct annot_str *annot_p, unsigned char *ann_arr, char *tmp_line) {
  unsigned char *aa1a_tmp;
  int i, ic, n_annot;
  struct annot_entry *this_annot;
  char *bp;

  if ((aa1a_tmp = (unsigned char *)calloc(n1+2,sizeof(char)))==NULL) {
    fprintf(stderr,"*** error [%s:%d] - cannot allocate aa1a_ann[%d] array\n",
	    __FILE__, __LINE__, n1);
    return;
  }

  if (offset < 0) offset++;

  for (i=0; i < annot_p->n_annot; i++) {
    this_annot = &annot_p->annot_arr_p[i];
    /* skip VAR labels */
    if (this_annot->label == 'V') { continue; }
    if (this_annot->label == '-') {
      aa1a_tmp[this_annot->pos]=qascii['['] - NANN;
      aa1a_tmp[this_annot->end]=qascii[']'] - NANN;
      continue;
    }
    if (strchr((char *)ann_arr, this_annot->label)==NULL) {continue;}
    if (this_annot->pos - offset < n1) {
      if (this_annot->pos >= offset) {	/* not an error, but annotation must be in range */
	aa1a_tmp[this_annot->pos - offset]=qascii[this_annot->label] - NANN;
      }
    }
    else {
      fprintf(stderr, "*** error [%s:%d] - this_annot->pos:[%ld - %ld] out of range: %d : %s\n",
	      __FILE__, __LINE__, this_annot->pos,offset, n1, tmp_line);
    }
  }
  annot_p->aa1_ann = aa1a_tmp;
}

/* save_best captures much of the complexity of saving the best scores
   and appropriately sampling the scores for statistical analysis. It
   does the following:

   (1) update s_info counts for functions like fasta/x/y that don't
       optimize every score

   (2) for every result in the buffer:
       (a) decide if it should be used for statistical sampling
       (b) if the number of samples > MAX_STATS, then run
           process_hist() and update all the zscores
       (c) reset everything for next sequence

   (3) must ensure that -BIGNUM are never in best[]

*/

#include "thr_buf_structs.h"
#ifndef PCOMPLIB
#define RESULTS_BUF reader_buf
#define XTERNAL
#include "thr_bufs2.h"
#else
#define RESULTS_BUF worker_buf
#include "pcomp_bufs.h"
#endif

extern char *prog_func;		/* function label */
extern int fa_max_workers;
extern struct buf_head *lib_buf2_list;
#ifdef DEBUG
void check_rbuf(struct buf_head *cur_buf);
#endif
extern void get_rbuf(struct buf_head **lib_buf, int max_work_buf);
extern void put_rbuf(struct buf_head *lib_buf, int max_work_buf);
extern void wait_rbuf(int max_work_buf);
extern void rbuf_done(int nthreads);
extern void put_rbuf_done(int nthreads, struct buf_head *lib_buf,
			  int max_work_buf);
extern int
process_hist(struct stat_str *sptr, int nstats,
	     const struct mngmsg *m_msg,
	     struct pstruct *ppst,
	     struct hist_str *hist, void **pstat_void, struct score_count_s *s_info, int do_hist);

extern void addhistz(double, struct hist_str *); /* scaleswn.c */
void selectbestz(struct beststr **, int, int );
extern double find_z(int score, double escore, int length, double comp, void *);
extern double zs_to_E(double zs,int n1, int dnaseq, long entries, struct db_str db);
extern struct beststr **bestp_arr;	/* array of pointers */
extern int nbest;
extern int nstats, nqstats, nrstats, pre_nstats, kstats, shuff_tot, sstats;
extern double zbestcut;	/* cut off for best z-score */
extern int bestfull;	/* index for selectbest() */
extern int stats_done;	/* flag for z-value processing */
extern void *rand_state;
extern struct stat_str *stats; /* array of scores for statistics from real
			   (or shuffled) sequences*/
extern struct stat_str *qstats;	/* array of scores for shuffled query stats */
extern struct stat_str *rstats;	/* array of scores from shuffled library */

/* in the current version (fasta_35_01) save_best is used by both
   threaded and unthreaded versions */

#define COPY_RST_P(d,s) 		\
{ d->rst.score[0] = s->rst.score[0];	\
  d->rst.score[1] = s->rst.score[1];	\
  d->rst.score[2] = s->rst.score[2];	\
  d->rst.valid_stat = s->rst.valid_stat; \
  d->rst.comp = s->rst.comp;		\
  d->rst.H = s->rst.H;			\
  d->rst.escore = s->rst.escore;	\
  d->rst.segnum = s->rst.segnum;	\
  d->rst.seglen = s->rst.seglen;	\
}

void
save_best(struct buf_head *lib_bhead_p,
	  const struct mngmsg *m_msp, struct pstruct *ppst,
	  struct db_str *ldb, FILE *fdata,
	  struct hist_str *histp, void **pstat_voidp,
	  struct score_count_s *s_info)
{
  double zscore;
  int i_score;
  struct beststr *bbp;
  struct buf2_data_s *rbuf_dp, *lib_buf2_dp;
  struct buf2_res_s *rbuf_rp, *lib_buf2_rp;
  int i, t_best, t_rbest, t_qrbest, tm_best, t_n1, sc_ix;
  int t_valid_stat, tr_valid_stat, use_shuff, zsflag_save;
  double e_score, tm_escore, t_rescore, t_qrescore;
  int buf2_cnt;

  if (!lib_bhead_p->hdr.have_results) return;
  if ((buf2_cnt=lib_bhead_p->hdr.buf2_cnt) <= 0) return;
  lib_buf2_dp = lib_bhead_p->buf2_data;
  lib_buf2_rp = lib_bhead_p->buf2_res;

  shuff_tot += lib_bhead_p->hdr.shuff_cnt;
  s_info->s_cnt[0] += lib_bhead_p->s_cnt_info.s_cnt[0];
  s_info->s_cnt[1] += lib_bhead_p->s_cnt_info.s_cnt[1];
  s_info->s_cnt[2] += lib_bhead_p->s_cnt_info.s_cnt[2];
  s_info->tot_scores += lib_bhead_p->s_cnt_info.tot_scores;;

  sc_ix = ppst->score_ix;

  t_best = t_rbest = t_qrbest = -BIGNUM;
  tm_escore = t_rescore = t_qrescore = FLT_MAX;
  t_valid_stat = tr_valid_stat = 0;
  if (ppst->zsflag >= 10 && ppst->zsflag < 20) { use_shuff = 1;}
  else { use_shuff = 0;}

#ifdef DEBUG
  if (fdata) {
    fprintf(fdata,">save_best: %d\n",buf2_cnt);
  }
#endif

  buf2_cnt = lib_bhead_p->hdr.buf2_cnt;
  while (buf2_cnt--) { /* count down the number of results */
    rbuf_rp = lib_buf2_rp++;	/* step through the results buffer */
    rbuf_dp = lib_buf2_dp++;	/* step through the data buffer */

    /* perhaps should use explicit flag to indicate no score */
    if (rbuf_rp->rst.score[0] == -BIGNUM) continue;

    /* i_score: current raw sorting score */
    i_score = rbuf_rp->rst.score[sc_ix];
    /* e_score, current escore */
    e_score = rbuf_rp->rst.escore;

    /* this should be done in the thread, and a sorted set of indexes
       should be produced by the thread, so we just go down the list
       to the zscore threshold */
    zscore = (double)i_score;
    if (stats_done) {
      zscore=find_z(i_score, e_score, rbuf_dp->seq->n1,(double)rbuf_rp->rst.comp,
			  *pstat_voidp);
    }

    /* we have complex logic to decide:
       (a) for multiframe results, which is the best
       (b) information about valid stats
       we should simply return a stats array where all this is figured
       out in the thread.
     */
    t_n1 = rbuf_dp->seq->n1;
    if (i_score > t_best) tm_best = t_best = i_score;
    if (e_score < tm_escore) tm_escore = e_score;
    if (rbuf_rp->rst.valid_stat > t_valid_stat) {
      t_valid_stat = 1;
    }

    /* this stuff happens only for fasts/fastm/fastf
       again, the t_qrbest stuff should be done in the thread
       rather than check for every hit, run through the loop
       only if necessary.
     */
    if (m_msp->qshuffle) {
      if (rbuf_rp->qr_score > t_qrbest)
	t_qrbest = rbuf_rp->qr_score;
      if (rbuf_rp->qr_escore < t_qrescore)
	t_qrescore = rbuf_rp->qr_escore;

      if (rbuf_dp->frame == m_msp->nitt1 && t_qrbest > 0 && nqstats < m_msp->shuff_max) {
	qstats[nqstats].n1 = rbuf_dp->seq->n1;	/* save the best score */
	qstats[nqstats].comp =  rbuf_rp->rst.comp;
	qstats[nqstats].H = rbuf_rp->rst.H;
	qstats[nqstats].escore = t_qrescore;
	qstats[nqstats++].score = t_qrbest;
	t_qrbest = -BIGNUM;	/* reset t_qrbest, t_qrescore */
	t_qrescore = FLT_MAX;
      }
    }	/* m_msp->qshuffle */

    if (use_shuff) {
      /* this check is required because some sequences scheduled to be
	 used for statistics may not in fact be returning a score (if
	 they are outside the -M range, for example.
       */
      if (rbuf_rp->r_rst.score[0] == -BIGNUM) { tr_valid_stat = 0; }
      if (rbuf_rp->r_rst.valid_stat > tr_valid_stat) {
	tr_valid_stat = 1;
      }
      if (rbuf_rp->r_rst.score[sc_ix] > t_rbest) {
	t_rbest = rbuf_rp->r_rst.score[sc_ix];
	t_rescore = rbuf_rp->r_rst.escore;
      }
    }

    /* need to look for frame 0 if TFASTA, then save stats at frame 6 */
    if (fdata) {
      fprintf(fdata,
	      "%-12s %6d %d %.5f %.5f %4d %4d %4d %2d %2d %4d %4d %4d %2d %2d %5d %8lld\n",
	      rbuf_dp->mseq->libstr, rbuf_dp->seq->n1,rbuf_dp->frame,rbuf_rp->rst.comp,rbuf_rp->rst.H,
	      rbuf_rp->rst.score[0],rbuf_rp->rst.score[1],rbuf_rp->rst.score[2],
	      t_valid_stat, rbuf_rp->rst.alg_info,
	      (rbuf_rp->r_rst.score[0]<0 ? -1 : rbuf_rp->r_rst.score[0]),
	      (rbuf_rp->r_rst.score[1]<0 ? -1 : rbuf_rp->r_rst.score[1]),
	      (rbuf_rp->r_rst.score[2]<0 ? -1 : rbuf_rp->r_rst.score[2]),
	      tr_valid_stat, rbuf_rp->r_rst.alg_info,
	      rbuf_dp->stats_idx, rbuf_dp->mseq->lseek);
    }

    /* statistics done for best score of set */

    if (rbuf_dp->frame == m_msp->nitt1) {
      ldb->entries++;
      ldb->length += t_n1;
      if (ldb->length > LONG_MAX) {
	ldb->length -= LONG_MAX; ldb->carry++;
      }
    }

    if (rbuf_dp->frame == m_msp->nitt1 && ppst->zsflag >= 0) {
      /* if this sample should be used for statistics */
      if (use_shuff) t_valid_stat = tr_valid_stat;
      if (t_valid_stat) {
	/* we've got our initial MAX_STATS values */
	if (nstats >= MAX_STATS) {
	  if (!stats_done) {
	    zsflag_save = ppst->zsflag;
	    if (ppst->zsflag > 20) {
	      ppst->zsflag -= 20;
	    }
	    ppst->zsflag_f = process_hist(stats,nstats,m_msp, ppst,
					  histp, pstat_voidp,s_info, 0);
	    ppst->zsflag = zsflag_save;
	    kstats = nstats;
	    if (ppst->zsflag >= 0) {	/* this is redundant, but rare */
	      stats_done = 1;
	      for (i=0; i< nstats; i++) {
		bestp_arr[i]->zscore =
		  find_z(bestp_arr[i]->rst.score[ppst->score_ix],
			 bestp_arr[i]->rst.escore, bestp_arr[i]->seq->n1,
			 bestp_arr[i]->rst.comp, *pstat_voidp);
	      }
	    }
	  }
	}
	else {
	  /* this logic allows stats_idx to be over-ruled for searches
	     where every query does not generate a score */
	  rbuf_dp->stats_idx = nstats;
	  nstats++;
	}
      }

      if (rbuf_dp->stats_idx >= 0 && t_valid_stat) {
	if (rbuf_dp->stats_idx >= MAX_STATS || nstats > MAX_STATS) {
	  fprintf(stderr, "*** error [%s:%d] - nstats index [%d] out of range [%d,%d]\n",
		  __FILE__, __LINE__,
		  rbuf_dp->stats_idx, nstats,MAX_STATS);
	}
	else {  /* stats_idx is in range */
	  sstats++;
	  stats[rbuf_dp->stats_idx].n1 = t_n1;
	  stats[rbuf_dp->stats_idx].comp = rbuf_rp->rst.comp;
	  stats[rbuf_dp->stats_idx].H = rbuf_rp->rst.H;
	  if (use_shuff) { /* use shuffled score */
	    stats[rbuf_dp->stats_idx].escore  = t_rescore;
	    stats[rbuf_dp->stats_idx].score = t_rbest;
	  }
	  else { /* real score, not shuffled */
	    stats[rbuf_dp->stats_idx].escore  = tm_escore;
	    stats[rbuf_dp->stats_idx].score = tm_best;
	  }
	} /* end stats_idx in range */
      }	/* end have valid stats_idx */

      if (t_valid_stat && stats_done && histp) {
	addhistz(find_z(t_best, tm_escore, rbuf_dp->seq->n1, (double) rbuf_rp->rst.comp,
			    *pstat_voidp), histp);
      }
      /* reset best scores */
      t_best = t_rbest = -BIGNUM;
      tm_escore = t_rescore = FLT_MAX;
      t_valid_stat = tr_valid_stat = 0;
    }

    /*
    if (rbuf_rp->rst.score[ppst->score_ix] > 200) {
      fprintf(stderr, "high score[%d]: %s %d: %d\n", rbuf_dp->seq->index,
	      rbuf_dp->mseq->libstr, rbuf_dp->seq->n1, rbuf_rp->rst.score[ppst->score_ix]);
    }
    */

    if (zscore > zbestcut) {
      if (nbest >= MAX_BEST) {
	bestfull = nbest-MAX_BEST/4;
	selectbestz(bestp_arr,bestfull-1,nbest);
	zbestcut = bestp_arr[bestfull-1]->zscore;
	nbest = bestfull;
      }
      bbp = bestp_arr[nbest++];

      COPY_RST_P(bbp, rbuf_rp);

      bbp->seq = rbuf_dp->seq;
      bbp->mseq = rbuf_dp->mseq;
      bbp->n1 = rbuf_dp->seq->n1;
#ifdef DEBUG
      bbp->adler32_crc = rbuf_dp->seq->adler32_crc;
#endif
      /* rbuf_dp->best_save is set after a rbuf_dp is entered into best_str */
      if (rbuf_dp->best_save) {
	/* a previous rbuf_dp->seq is in best_str at best_save */
	if (rbuf_dp->best_save->seq == rbuf_dp->seq) {
	  /* the best_save->seq matches the rbuf_dp->seq */
	  bbp->bbp_link = rbuf_dp->best_save;
	  /* bbp_link tells where this ->seq can be found */
	}
	else {
	  bbp->bbp_link = NULL;
	}
      }
      rbuf_dp->best_save = bbp;
      lib_bhead_p->hdr.have_best_save = 1;
      bbp->zscore = zscore;
      bbp->frame = rbuf_dp->frame;
    }
  }
}

void
save_best2(struct buf_head *lib_bhead_p,
	   const struct mngmsg *m_msp, struct pstruct *ppst,
	   struct db_str *ldb, FILE *fdata,
	   struct hist_str *histp, void **pstat_voidp,
	   struct score_count_s *s_info)
{
  double zscore;
  int i_score;
  struct beststr *bbp;
  struct buf2_data_s *rbuf_dp, *lib_buf2_dp;
  struct buf2_res_s *rbuf_rp, *lib_buf2_rp;
  int i, sc_ix;
  int t_valid_stat, use_shuff, zsflag_save;
  double e_score;
  int buf2_cnt;

  if (!lib_bhead_p->hdr.have_results) return;
  if ((buf2_cnt = lib_bhead_p->hdr.buf2_cnt) <= 0) return;

  /*
#ifdef DEBUG
  fprintf(stderr," save_best2: lib_bhead_p->buf2_data[0]->mseq->index/lseek: %d,%lld\n",
	  lib_bhead_p->buf2_data[0].mseq->index,lib_bhead_p->buf2_data[0].mseq->lseek);
#endif
  */
  if (ppst->zsflag >= 10 && ppst->zsflag < 20) { use_shuff = 1;}
  else {use_shuff = 0;}

  shuff_tot += lib_bhead_p->hdr.shuff_cnt;
  s_info->s_cnt[0] += lib_bhead_p->s_cnt_info.s_cnt[0];
  s_info->s_cnt[1] += lib_bhead_p->s_cnt_info.s_cnt[1];
  s_info->s_cnt[2] += lib_bhead_p->s_cnt_info.s_cnt[2];
  s_info->tot_scores += lib_bhead_p->s_cnt_info.tot_scores;;
  sc_ix = ppst->score_ix;

  /* save the raw data if requested */
  if (fdata) {
#ifdef DEBUG
    fprintf(fdata,">save_best: %d\n",buf2_cnt);
#endif
    lib_buf2_dp = lib_bhead_p->buf2_data;
    lib_buf2_rp = lib_bhead_p->buf2_res;
    buf2_cnt = lib_bhead_p->hdr.buf2_cnt;
    while (buf2_cnt--) { /* count down the number of results */
      
      rbuf_rp = lib_buf2_rp++;	/* step through the results buffer */
      rbuf_dp = lib_buf2_dp++;	/* step through the data buffer */

      /* perhaps should use explicit flag to indicate no score */
      if (rbuf_rp->rst.score[0] == -BIGNUM) continue;

      fprintf(fdata,
	      "%-12s %6d %d %.5f %.5f %4d %4d %4d %2d %2d %4d %4d %4d %2d %2d %5d %8lld\n",
	      rbuf_dp->mseq->libstr, rbuf_dp->seq->n1,rbuf_dp->frame,rbuf_rp->rst.comp,rbuf_rp->rst.H,
	      rbuf_rp->rst.score[0],rbuf_rp->rst.score[1],rbuf_rp->rst.score[2],
	      rbuf_rp->is_valid_stat, rbuf_rp->rst.alg_info,
	      (rbuf_rp->r_rst.score[0]<0 ? -1 : rbuf_rp->r_rst.score[0]),
	      (rbuf_rp->r_rst.score[1]<0 ? -1 : rbuf_rp->r_rst.score[1]),
	      (rbuf_rp->r_rst.score[2]<0 ? -1 : rbuf_rp->r_rst.score[2]),
	      rbuf_rp->is_valid_stat, rbuf_rp->r_rst.alg_info,
	      rbuf_dp->stats_idx, rbuf_dp->mseq->lseek);
    }
  }

  /* save the high-scoring data */
  lib_buf2_rp = lib_bhead_p->buf2_res;
  lib_buf2_dp = lib_bhead_p->buf2_data;
  buf2_cnt = lib_bhead_p->hdr.buf2_cnt;
  while (buf2_cnt--) {
    rbuf_rp = lib_buf2_rp++;
    rbuf_dp = lib_buf2_dp++;

    /* perhaps should use explicit flag to indicate no score */
    if (rbuf_rp->rst.score[0] == -BIGNUM) continue;

    /* i_score: current raw sorting score */
    i_score = rbuf_rp->rst.score[sc_ix];
    /* e_score, current escore */
    e_score = rbuf_rp->rst.escore;
    /* this should be done in the thread, and a sorted set of indexes
       should be produced by the thread, so we just go down the list
       to the zscore threshold */
    zscore = (double)i_score;
    if (stats_done) {
      zscore=find_z(i_score, e_score, rbuf_dp->seq->n1,(double)rbuf_rp->rst.comp,
			  *pstat_voidp);
    }

    if (rbuf_dp->frame == m_msp->nitt1) {
      ldb->entries++;
      ldb->length += rbuf_dp->seq->n1;
      if (ldb->length > LONG_MAX) {
	ldb->length -= LONG_MAX; ldb->carry++;
      }
    }

    if (zscore > zbestcut) {
      if (nbest >= MAX_BEST) {
	bestfull = nbest-MAX_BEST/4;
	selectbestz(bestp_arr,bestfull-1,nbest);
	zbestcut = bestp_arr[bestfull-1]->zscore;
	nbest = bestfull;
      }
      bbp = bestp_arr[nbest++];

      COPY_RST_P(bbp, rbuf_rp);

      bbp->seq = rbuf_dp->seq;
      bbp->mseq = rbuf_dp->mseq;
      bbp->n1 = rbuf_dp->seq->n1;
#ifdef DEBUG
      bbp->adler32_crc = rbuf_dp->seq->adler32_crc;
#endif
      /* rbuf_dp->best_save is set after a rbuf_dp is entered into best_str */
      if (rbuf_dp->best_save) {
	/* a previous rbuf_dp->seq is in best_str at best_save */
	if (rbuf_dp->best_save->seq == rbuf_dp->seq) {
	  /* the best_save->seq matches the rbuf_dp->seq */
	  bbp->bbp_link = rbuf_dp->best_save;
	  /* bbp_link tells where this ->seq can be found */
	}
	else {
	  bbp->bbp_link = NULL;
	}
      }
      rbuf_dp->best_save = bbp;
      lib_bhead_p->hdr.have_best_save = 1;
      bbp->zscore = zscore;
      bbp->frame = rbuf_dp->frame;
    }
  }

  /* process results for statistics */
  lib_buf2_dp = lib_bhead_p->buf2_data;
  lib_buf2_rp = lib_bhead_p->buf2_res;
  buf2_cnt = lib_bhead_p->hdr.buf2_cnt;
  while (buf2_cnt--) { /* count down the number of results */
    rbuf_dp = lib_buf2_dp++;	/* step through the results buffer */
    rbuf_rp = lib_buf2_rp++;	/* step through the results buffer */

    if (!rbuf_rp->is_valid_stat) { continue;}


    if (use_shuff) {
      i_score = rbuf_rp->r_rst.score[sc_ix];
    e_score = rbuf_rp->r_rst.escore;
    }
    else {
      i_score = rbuf_rp->rst.score[sc_ix];
      e_score = rbuf_rp->rst.escore;
    }

    if (rbuf_dp->stats_idx >= MAX_STATS || nstats > MAX_STATS) {
      fprintf(stderr, "*** error [%s:%d] - nstats index [%d] out of range [%d,%d]\n",
	      __FILE__, __LINE__,
	      rbuf_dp->stats_idx, nstats,MAX_STATS);
      continue;
    }

    if (nstats < MAX_STATS) {
      /* this logic allows stats_idx to be over-ruled for searches
	 where every query does not generate a score */
      rbuf_dp->stats_idx = nstats;
      nstats++;
    }

    if (stats_done && histp) {
      addhistz(find_z(i_score, e_score, rbuf_dp->seq->n1, (double) rbuf_rp->rst.comp,
		      *pstat_voidp), histp);
    }

    if (rbuf_dp->stats_idx < 0) {
      continue;
    }

    sstats++;
    stats[rbuf_dp->stats_idx].n1 = rbuf_dp->seq->n1;
    stats[rbuf_dp->stats_idx].comp = rbuf_rp->rst.comp;
    stats[rbuf_dp->stats_idx].H = rbuf_rp->rst.H;
    stats[rbuf_dp->stats_idx].escore  = e_score;
    stats[rbuf_dp->stats_idx].score = i_score;
  }


  /* fill the qstats[] array if m_msp->qshuffle */
  if (m_msp->qshuffle && nqstats < m_msp->shuff_max) {
    lib_buf2_dp = lib_bhead_p->buf2_data;
    lib_buf2_rp = lib_bhead_p->buf2_res;
    buf2_cnt = lib_bhead_p->hdr.buf2_cnt;
    while (buf2_cnt--) {
      rbuf_rp = lib_buf2_rp++;	/* step through the results buffer */
      rbuf_dp = lib_buf2_dp++;

      if (rbuf_rp->is_valid_stat && rbuf_rp->qr_score > 0
	  && nqstats < m_msp->shuff_max) {
	qstats[nqstats].n1 = rbuf_dp->seq->n1;	/* save the best score */
	qstats[nqstats].comp =  rbuf_rp->rst.comp;
	qstats[nqstats].H = rbuf_rp->rst.H;
	qstats[nqstats].escore = rbuf_rp->qr_escore;
	qstats[nqstats++].score = rbuf_rp->qr_score;
      }
    }	/* m_msp->qshuffle */
  }

  /* check if we have enough data to do stats */
  if (!stats_done && nstats >= MAX_STATS) {
    zsflag_save = ppst->zsflag;
    if (ppst->zsflag > 20) {
      ppst->zsflag -= 20;
    }
    ppst->zsflag_f = process_hist(stats,nstats,m_msp, ppst,
				  histp, pstat_voidp,s_info, 0);
    ppst->zsflag = zsflag_save;
    kstats = nstats;
    stats_done = 1;
    for (i=0; i< nstats; i++) {
      bestp_arr[i]->zscore =
	find_z(bestp_arr[i]->rst.score[ppst->score_ix],
	       bestp_arr[i]->rst.escore, bestp_arr[i]->seq->n1,
	       bestp_arr[i]->rst.comp, *pstat_voidp);
    }
  }

}

void
save_shuf(struct buf_head *lib_bhead_p, int nitt1, int shuff_max, int sc_ix,
	  struct score_count_s *s_info)
{
  struct buf2_data_s *rbuf_dp, *lib_buf2_dp;
  struct buf2_res_s *rbuf_rp, *lib_buf2_rp;
  int t_valid_stat;
  int t_rbest;
  double t_rescore;
  int buf2_cnt, jstats;
  static int kstats=0;


  lib_buf2_dp = lib_bhead_p->buf2_data;
  lib_buf2_rp = lib_bhead_p->buf2_res;
  buf2_cnt = lib_bhead_p->hdr.buf2_cnt;

  s_info->s_cnt[0] += lib_bhead_p->s_cnt_info.s_cnt[0];
  s_info->s_cnt[1] += lib_bhead_p->s_cnt_info.s_cnt[1];
  s_info->s_cnt[2] += lib_bhead_p->s_cnt_info.s_cnt[2];

  s_info->tot_scores += lib_bhead_p->s_cnt_info.tot_scores;
  /* this is done because we are not using r_rst->valid_stat to limit selection of scores */
  /*   s_info->s_cnt[sc_ix] = s_info->tot_scores; */

  t_rbest = -BIGNUM;
  t_valid_stat = 0;

  while (buf2_cnt--) { /* count down the number of results */
    rbuf_dp = lib_buf2_dp++;	/* step through the results buffer */
    rbuf_rp = lib_buf2_rp++;	/* step through the results buffer */

    /* perhaps should use explicit flag to indicate no score */
    if (rbuf_rp->r_rst.score[0] == -BIGNUM) continue;

    if (rbuf_rp->r_rst.score[sc_ix] > t_rbest) {
      t_rbest = rbuf_rp->r_rst.score[sc_ix];
      t_rescore = rbuf_rp->r_rst.escore;
    }

    if (rbuf_rp->r_rst.valid_stat > t_valid_stat) {
      t_valid_stat = 1;
    }

    /* statistics done for best score of set */
    /* currently no check for rst->valid_stat, which causes
       over-estimates of shuffles */

    if (rbuf_dp->frame == nitt1) {
      if (t_valid_stat) {
	if (nrstats < shuff_max ) { kstats = jstats = nrstats++; }
	else {	/* randomly replace */
	  jstats = my_nrand(++kstats,rand_state);
	  if (jstats >= shuff_max) goto done;
	}

	rstats[jstats].n1 = rbuf_dp->seq->n1;
	rstats[jstats].comp = rbuf_rp->r_rst.comp;
	rstats[jstats].H = rbuf_rp->r_rst.H;
	rstats[jstats].escore  = t_rescore;
	rstats[jstats].score = t_rbest;
      done:
	t_rbest = -BIGNUM;
      }
    }
  }
}

int
save_align(struct buf_head *lib_bhead_p, struct beststr **bestp_arr)
{
  struct buf2_ares_s *rbuf_ap, *lib_buf2_ap;
  int buf2_cnt;

  if (!lib_bhead_p->hdr.have_results || lib_bhead_p->hdr.buf2_cnt <= 0) return 0;

  lib_buf2_ap = lib_bhead_p->buf2_ares;
  buf2_cnt = lib_bhead_p->hdr.buf2_cnt;

  while (buf2_cnt-- > 0) { /* count down the number of results */
    rbuf_ap = lib_buf2_ap++;	/* step through the results buffer */
    if (bestp_arr[rbuf_ap->best_idx]->a_res == NULL) {
      bestp_arr[rbuf_ap->best_idx]->have_ares = rbuf_ap->have_ares;
      bestp_arr[rbuf_ap->best_idx]->a_res = rbuf_ap->a_res;
    }
#ifdef DEBUG
    else {
      fprintf(stderr,"*** error [%s:%d] - attempt to re-save a_res for [%d]: %s\n",
	      __FILE__, __LINE__, rbuf_ap->best_idx, bestp_arr[rbuf_ap->best_idx]->mseq->bline);
    }
#endif
  }

  buf2_cnt = lib_bhead_p->hdr.buf2_cnt;
  lib_bhead_p->hdr.have_results = 0;
  lib_bhead_p->hdr.buf2_cnt = 0;
  return buf2_cnt;
}

/* buf_do_work fills in the lib_bhead_p->buf2_res[] array with the
   do_work() results,

   inputs:  **aa0, n0 (query)
            lib_bhead_p->buf2_data lib_bhead_p->hdr.buf2_cnt library sequences
	    max_frame (used to set statistics info)
	    ppst,
	    void *f_struct prepared by init_work()

   results: lib_bhead_p->buf2_res[]

            included in buf2_res[] is use_stat, which captures the
            logic required to decide whether a value should be saved
            in the stats[] buffer.  This complexity mostly arises
            because there can be more scores than sequences, but there
            can only on statistics score per sequence (the best score).
*/
void
buf_do_work(unsigned char **aa0,  int n0,
	    struct buf_head *lib_bhead_p,
	    int max_frame,
	    struct pstruct *ppst, void **f_str) {

  int buf2_cnt;
  unsigned long atmp;
  struct buf2_data_s *lib_buf2_dp;
  struct buf2_res_s *lib_buf2_rp, *t_best_rp;
  int t_best, sc_ix;
  double t_escore;

  sc_ix = ppst->score_ix;

  lib_bhead_p->s_cnt_info.s_cnt[0] = lib_bhead_p->s_cnt_info.s_cnt[1] =
    lib_bhead_p->s_cnt_info.s_cnt[2] = lib_bhead_p->s_cnt_info.tot_scores = 0;

  buf2_cnt = lib_bhead_p->hdr.buf2_cnt;
  lib_buf2_dp = lib_bhead_p->buf2_data;
  lib_buf2_rp = lib_bhead_p->buf2_res;

  t_best_rp = NULL;
  t_best = -BIGNUM;
  t_escore = 1000.0;

  while (buf2_cnt-- > 0) {

    lib_buf2_rp->rst.score[0] =
      lib_buf2_rp->rst.score[1] =
      lib_buf2_rp->rst.score[2] = -BIGNUM;

    lib_buf2_rp->is_valid_stat = 0;

    if (lib_buf2_dp->seq->n1 < ppst->n1_low ||
	lib_buf2_dp->seq->n1 > ppst->n1_high ) {
      /* tells save_best() there is no stats score here -- not
	 necessary as -BIGNUM indicates no score */
      lib_buf2_dp->stats_idx = -1;
      goto next_seq;
    }

#ifdef DEBUG
    if (check_seq_range(lib_buf2_dp->seq->aa1b, lib_buf2_dp->seq->n1,
			ppst->nsqx, "buf_do_work()")) {
      fprintf(stderr, "*** error [%s:%d] - [%s/buf_do_work] range error at: %d/%d (n1:%d)\n",
	      __FILE__, __LINE__,
	      prog_func,lib_bhead_p->hdr.buf2_cnt - (buf2_cnt+1),
	      lib_bhead_p->hdr.buf2_cnt, lib_buf2_dp->seq->n1);
      goto next_seq;
    };

    /* also check for adler32_crc match */
    if (lib_buf2_dp->seq->adler32_crc != (atmp=adler32(1L,lib_buf2_dp->seq->aa1b, lib_buf2_dp->seq->n1))) {
      fprintf(stderr, "*** error [%s:%d] - [%s/buf_do_work] CRC error [%lu!=%lu] at: %d/%d (n1:%d/l_offset:%ld)\n",
	      __FILE__, __LINE__,
	      prog_func,lib_buf2_dp->seq->adler32_crc, atmp,
	      lib_bhead_p->hdr.buf2_cnt - (buf2_cnt+1),
	      lib_bhead_p->hdr.buf2_cnt,lib_buf2_dp->seq->n1,
	      lib_buf2_dp->seq->l_offset);
      goto next_seq;
    }
#endif

    do_work (aa0[lib_buf2_dp->frame], n0,
	     lib_buf2_dp->seq->aa1b, lib_buf2_dp->seq->n1,
	     lib_buf2_dp->frame, ppst, f_str[lib_buf2_dp->frame], 0, 0,
	     &(lib_buf2_rp->rst), &(lib_bhead_p->s_cnt_info));

    if (lib_buf2_rp->rst.valid_stat) {
      if (lib_buf2_rp->rst.escore < t_escore) {
	t_escore = lib_buf2_rp->rst.escore;
	t_best_rp = lib_buf2_rp;
      }
      if (lib_buf2_rp->rst.score[sc_ix] > t_best) {
	t_best = lib_buf2_rp->rst.score[sc_ix];
	t_best_rp = lib_buf2_rp;
      }
    }

    if (lib_buf2_dp->frame == max_frame) {
      if (t_best_rp!=NULL) {
	t_best_rp->is_valid_stat = 1;
	t_best_rp = NULL;
      }
      t_best = -BIGNUM;
      t_escore = 1000.0;
    }

  next_seq:
    lib_buf2_dp++;
    lib_buf2_rp++;
  }

  /* place to produce z_scores */
  /* place to produce sorted array */

  lib_bhead_p->hdr.have_results = 1;
}

void
buf_do_align(unsigned char **aa0,  int n0,
	     struct buf_head *lib_bhead_p,
	     struct pstruct *ppst, const struct mngmsg *m_msp,
	     void **f_str) {

  int buf2_cnt, i, nsq;
  struct buf2_data_s *lib_buf2_dp;
  struct buf2_res_s *lib_buf2_rp;
  struct buf2_ares_s *lib_buf2_ap;
  struct rstruct rst;

  buf2_cnt = lib_bhead_p->hdr.buf2_cnt;
  lib_buf2_dp = lib_bhead_p->buf2_data;
  lib_buf2_rp = lib_bhead_p->buf2_res;
  lib_buf2_ap = lib_bhead_p->buf2_ares;

  while (buf2_cnt-- > 0) {
    if ( m_msp->stages > 1) {
      /* this is not typically done unless m_msp->stages > 1 */
      do_opt (aa0[lib_buf2_dp->frame], n0, lib_buf2_dp->seq->aa1b, lib_buf2_dp->seq->n1,
	      lib_buf2_dp->frame, ppst, f_str[lib_buf2_dp->frame], &rst);
      lib_buf2_rp->rst.score[2]=rst.score[2];
    }

#ifdef DEBUG
    if (lib_buf2_dp->seq->aa1b == NULL) {
      fprintf(stderr,"*** error [%s:%d] - [buf_do_align] null aa1b\n",__FILE__, __LINE__);
      lib_buf2_ap->a_res = NULL;
      break;
    }
    if (check_seq_range(lib_buf2_dp->seq->aa1b, lib_buf2_dp->seq->n1,
			ppst->nsqx, "buf_do_align()")) {
      fprintf(stderr, "*** error [%s:%d] - [%s/buf_do_align] range error at: %d/%d (n1:%d)\n",
	      __FILE__, __LINE__,
	      prog_func,lib_bhead_p->hdr.buf2_cnt - (buf2_cnt+1),
	      lib_bhead_p->hdr.buf2_cnt, lib_buf2_dp->seq->n1);
    };

    /* also check for adler32_crc match */
    if (lib_buf2_dp->seq->adler32_crc != adler32(1L,lib_buf2_dp->seq->aa1b, lib_buf2_dp->seq->n1)) {
      fprintf(stderr, "*** error [%s:%d] - [%s/buf_do_align] CRC error at: %d/%d (n1:%d)\n",
	      __FILE__, __LINE__,
	      prog_func,lib_bhead_p->hdr.buf2_cnt - (buf2_cnt+1),
	      lib_bhead_p->hdr.buf2_cnt, lib_buf2_dp->seq->n1);
    }
#endif

    lib_buf2_ap->a_res = build_ares_code(aa0[lib_buf2_dp->frame], m_msp->n0,
					 lib_buf2_dp->seq->aa1b, lib_buf2_dp->seq,
					 lib_buf2_dp->frame, &lib_buf2_ap->have_ares,
					 lib_buf2_dp->repeat_thresh, m_msp, ppst, f_str[lib_buf2_dp->frame] );

    lib_buf2_dp++;
    lib_buf2_ap++;
    lib_buf2_rp++;
  }
  lib_bhead_p->hdr.have_results = 1;
}

void
buf_qshuf_work(unsigned char *aa0s,  int n0,
	       struct buf_head *lib_bhead_p,
	       int max_frame,
	       struct pstruct *ppst, void *qf_str,
	       int ix_score)
{
  int buf2_cnt;
  struct buf2_data_s *lib_buf2_dp;
  struct buf2_res_s *lib_buf2_rp, *tq_best_rp;
  struct rstruct rrst;
  struct score_count_s q_scnt_info;
  int tq_best;
  double tq_escore;

  buf2_cnt = lib_bhead_p->hdr.buf2_cnt;
  lib_buf2_dp = lib_bhead_p->buf2_data;
  lib_buf2_rp = lib_bhead_p->buf2_res;

  tq_best_rp = NULL;
  tq_best = -BIGNUM;
  tq_escore = 1000.0;

  while (buf2_cnt-- > 0) {
    rrst.score[0] = rrst.score[1] = rrst.score[2] = -BIGNUM;
    rrst.valid_stat = 0;

    if (lib_buf2_dp->seq->n1 < ppst->n1_low ||
	lib_buf2_dp->seq->n1 > ppst->n1_high ) {
      lib_buf2_dp++;
      lib_buf2_rp++;
      tq_best_rp = NULL;
      tq_best = -BIGNUM;
      tq_escore = 1000.0;
      continue;
    }

    do_work (aa0s, n0,
	     lib_buf2_dp->seq->aa1b, lib_buf2_dp->seq->n1,
	     lib_buf2_dp->frame, ppst, qf_str, 1, 0,
	     &rrst, &q_scnt_info);

    /* buf_qshuf_work() is always called after buf_do_work(), which
       sets rp->is_valid_stat */
    if (lib_buf2_rp->is_valid_stat) {
      tq_best_rp = lib_buf2_rp;
    }

    if (rrst.escore < tq_escore) {
      tq_escore = rrst.escore;
    }
    if (rrst.score[ix_score] > tq_best) {
      tq_best = rrst.score[ix_score];
    }

    if (lib_buf2_dp->frame == max_frame) {
      if (tq_best_rp!=NULL) {
	tq_best_rp->qr_score = tq_best;
	tq_best_rp->qr_escore = tq_escore;
	tq_best_rp = NULL;
      }
#ifdef DEBUG
      else {
	fprintf(stderr,"*** error [%s:%d] - tq_best_rp NULL at: %ld\n",
		__FILE__, __LINE__, lib_buf2_rp - lib_bhead_p->buf2_res);
      }
#endif
      tq_best = -BIGNUM;
      tq_escore = 1000.0;
    }
    lib_buf2_dp++;
    lib_buf2_rp++;
  }
}

void
buf_shuf_work(unsigned char **aa0,  int n0, unsigned char *aa1s, struct buf_head *lib_bhead_p,
	      int max_frame, struct pstruct *ppst, void **f_str,
	      int ix_score, void *rand_state)
{
  int buf2_cnt;
  int shuff_cnt;
  struct buf2_data_s *lib_buf2_dp;
  struct buf2_res_s *lib_buf2_rp, *tr_best_rp;
  int tr_best, sc_ix;
  double tr_escore;

  sc_ix = ppst->score_ix;

  lib_bhead_p->s_cnt_info.s_cnt[0] = lib_bhead_p->s_cnt_info.s_cnt[1] =
    lib_bhead_p->s_cnt_info.s_cnt[2] = lib_bhead_p->s_cnt_info.tot_scores = 0;

  shuff_cnt = 0;
  buf2_cnt = lib_bhead_p->hdr.buf2_cnt;
  lib_buf2_dp = lib_bhead_p->buf2_data;
  lib_buf2_rp = lib_bhead_p->buf2_res;

  tr_best_rp = NULL;
  tr_best = -BIGNUM;
  tr_escore = 1000.0;

  while (buf2_cnt-- > 0) {
    lib_buf2_rp->r_rst.score[0] = lib_buf2_rp->r_rst.score[1] =
      lib_buf2_rp->r_rst.score[2] = -BIGNUM;
    lib_buf2_rp->r_rst.valid_stat = lib_buf2_rp->is_valid_stat = 0;

    if ((lib_buf2_dp->stats_idx < 0) || lib_buf2_dp->seq->n1 < ppst->n1_low ||
	lib_buf2_dp->seq->n1 > ppst->n1_high ) {
      lib_buf2_dp++;
      lib_buf2_rp++;
      tr_best_rp = NULL;
      tr_best = -BIGNUM;
      tr_escore = 1000.0;
      continue;
    }

    shuff_cnt++;
    if (ppst->zs_win > 0) {
      wshuffle(lib_buf2_dp->seq->aa1b,aa1s,lib_buf2_dp->seq->n1,ppst->zs_win, rand_state);
    }
    else {
      if (ppst->shuffle_dna3) {shuffle3(lib_buf2_dp->seq->aa1b,aa1s,lib_buf2_dp->seq->n1, rand_state);}
      else {shuffle(lib_buf2_dp->seq->aa1b,aa1s,lib_buf2_dp->seq->n1, rand_state);}
    }

    /* rshuffle(lib_buf2_dp->seq->aa1b,aa1s,lib_buf2_dp->seq->n1); */

#ifdef DEBUG
    if (check_seq_range(aa1s, lib_buf2_dp->seq->n1,
			ppst->nsqx, "buf_do_align()")) {
      fprintf(stderr, "*** error [%s:%d] - [%s/buf_do_shuff] range error at: %d/%d (n1:%d)\n",
	      __FILE__, __LINE__,
	      prog_func,lib_bhead_p->hdr.buf2_cnt - (buf2_cnt+1),
	      lib_bhead_p->hdr.buf2_cnt, lib_buf2_dp->seq->n1);
    };
#endif

    do_work (aa0[lib_buf2_dp->frame], n0,
	     aa1s, lib_buf2_dp->seq->n1,
	     lib_buf2_dp->frame, ppst, f_str[lib_buf2_dp->frame], 0, 1,
	     &lib_buf2_rp->r_rst, &(lib_bhead_p->s_cnt_info));

    if (lib_buf2_rp->r_rst.valid_stat) {
      if (lib_buf2_rp->r_rst.escore < tr_escore) {
	tr_escore = lib_buf2_rp->r_rst.escore;
	tr_best_rp = lib_buf2_rp;
      }
      if (lib_buf2_rp->r_rst.score[sc_ix] > tr_best) {
	tr_best = lib_buf2_rp->r_rst.score[sc_ix];
	tr_best_rp = lib_buf2_rp;
      }
    }

    if (lib_buf2_dp->frame == max_frame) {
      if (tr_best_rp!=NULL) {
	tr_best_rp->is_valid_stat = 1;
	tr_best_rp = NULL;
      }
      tr_best = -BIGNUM;
      tr_escore = 1000.0;
    }

    lib_buf2_dp++;
    lib_buf2_rp++;
  }
  lib_bhead_p->hdr.shuff_cnt = shuff_cnt;
  lib_bhead_p->hdr.have_results = 1;
}

/* buf_shuf_seq is designed to:
   (1) take a list of sequences (specified by bptr[])
   (2) collect them from the database if they are not already available
   (3) send them to the threads or shuffle them directly and calculate scores
*/

void
buf_shuf_seq(unsigned char **aa0, int n0,
	     unsigned char **aa1shuff_b, unsigned char *aa1save, int maxn,
	     struct beststr **bestp_arr, int nbest,
	     struct pstruct *ppst, struct mngmsg *m_msp,
	     struct mng_thr *m_bufi_p
#if !defined(COMP_THR) && !defined(PCOMPLIB)
	     , void **f_str
#endif
	     , struct score_count_s *s_info)
{
  unsigned char *aa1shuff;
  struct beststr *bbp, **tmp_bestp;
  char l_bline[MAX_SSTR];
  int n1lib_req, shuff_mult;
  long loffset, l_off;
  int n1, itt;
  int max_do_cnt, ndiff, prev_index;
  int istats;
  int i, j;

  /* these variables track buffers of library sequences */
  int cur_buf_size, max_buf_size;
  struct buf_head *lib_bhead_p;
  struct buf2_data_s *lib_buf2_dp;
  struct buf2_res_s *lib_buf2_rp;

/* (1) get the sequences into a buffer - the sequence information is
   currently in the bestp_arr - find out how many we have, and how
   many we will need - the number to shuffle */

/* figure out how much space we need, first checking whether we have
   dups */
  if ((tmp_bestp = (struct beststr **)calloc(nbest, sizeof(struct beststr *)))==NULL) {
    fprintf(stderr,"*** error [%s:%d] - %s/buf_shuf_seq() *** cannot allocate tmp_bestp[%d]\n",
	    __FILE__, __LINE__, prog_name, nbest);
    exit(1);
  }
  for (i = 0; i < nbest; i++) {
    tmp_bestp[i] = bestp_arr[i];
  }

  /* sort tmp_bestp[] by sequence index, so duplicates are adjacent */
  sortbesti(tmp_bestp, nbest);

  /* count number of different sequence indices, get required space
     without dups */
  prev_index = -1;
  n1lib_req = ndiff = 0;
  for (i = 0; i < nbest; i++) {
    if (tmp_bestp[i]->seq->index > prev_index) {
      prev_index = tmp_bestp[i]->seq->index;
      n1lib_req += tmp_bestp[i]->n1+ 2;
      ndiff++;
    }
  }

#if !defined(COMP_THR) && !defined(PCOMPLIB)
      if (n1lib_req >= maxn) { /* we need new space, aa1shuff is too small */
	if ((*aa1shuff_b = aa1shuff =
	     (unsigned char *)realloc(*aa1shuff_b, n1lib_req*sizeof(char)))==NULL) {
	  fprintf(stderr,"*** error [%s:%d] - cannot realloc aa1shuff[%d]\n",
		  __FILE__, __LINE__, n1lib_req);
	  exit(1);
	}
      }
      else { aa1shuff = *aa1shuff_b;}
      *aa1shuff = '\0';
      aa1shuff++;

#else
      if (n1lib_req < 2) {
	fprintf(stderr,"*** error [%s:%d] - [%s/buf_shuf_seq] no residues to shuffle: %d (%d)\n",
		__FILE__, __LINE__,
		prog_func,n1lib_req,ndiff);
	exit(1);
      }

      if ((*aa1shuff_b = aa1shuff =
	   (unsigned char *)calloc(n1lib_req,sizeof(char)))==NULL) {
	fprintf(stderr,"*** error [%s:%d] - cannot calloc aa1shuff[%d]\n",
		__FILE__, __LINE__, n1lib_req);
	exit(1);
      }
      *aa1shuff = '\0';
      aa1shuff++;
#endif

      shuff_mult = (m_msp->shuff_max+1)/ndiff;
      istats = 0;

      /* setup lib_bhead buffers for shuffle comparisons */
#if defined(COMP_THR) || defined(PCOMPLIB)	/* threaded/parallel */
      /* max_do_cnt can be smaller than max_buf2_cnt, but not larger */
      max_do_cnt = min(m_bufi_p->max_buf2_res,
		       m_msp->shuff_max / (2 * fa_max_workers));
      /* we don't have a left over one, so we need one */
      get_rbuf(&lib_bhead_p,m_bufi_p->max_work_buf);
#else	/* not threaded */
      max_do_cnt = m_bufi_p->max_buf2_res;
      lib_bhead_p = lib_buf2_list;  /* equivalent to un-threaded get_rbuf() */
#endif
      max_buf_size = n1lib_req;
      cur_buf_size = 0;
      lib_bhead_p->hdr.buf2_cnt = 0;
      lib_bhead_p->hdr.have_results = 0;
      lib_bhead_p->hdr.stop_work = 0;
      lib_bhead_p->hdr.buf2_type=BUF2_DOSHUF;
      lib_bhead_p->hdr.seq_record_continuous = 0;
      lib_buf2_dp = lib_bhead_p->buf2_data;
      lib_buf2_rp = lib_bhead_p->buf2_res;

      /* read sequences into shuffle buffer */

      for (i = 0; i < ndiff; i++) {
	bbp = tmp_bestp[i];
	if (bbp->seq->aa1b == NULL) {
	  /* get the sequence */
	  (bbp->mseq->m_file_p->ranlib)(l_bline, sizeof(l_bline),
				       bbp->mseq->lseek,bbp->mseq->libstr,bbp->mseq->m_file_p);
	  n1 = re_getlib(aa1save,NULL, maxn,m_msp->ldb_info.maxt3,
			 m_msp->ldb_info.l_overlap,bbp->mseq->cont,m_msp->ldb_info.term_code,
			 &loffset,&l_off,bbp->mseq->m_file_p);

	  /* fprintf(stderr, " %d gets %d %d\n",i,tmp_bestp[i]->seq->n1,n1); */

	  memcpy(aa1shuff, aa1save, n1+1);
	  bbp->seq->aa1b = aa1shuff;
	  aa1shuff += n1 + 1;
	}

	/* lib_buf2_dp is used up by scores, the sequence is not sent multiple times */
	cur_buf_size += bbp->seq->n1+1;
	for (j = 0; j < shuff_mult; j++ ) {
	  for (itt = m_msp->revcomp; itt <= m_msp->nitt1; itt++) {
#ifdef PCOMPLIB
	    lib_buf2_dp->seq_dup = 0;	/* mark first ->seq as original, not duplicate */
#endif
	    lib_buf2_dp->seq = bbp->seq;
	    /* this invalidates lib_buf2_p->seq */
	    lib_buf2_dp->stats_idx = istats++;
	    lib_buf2_dp->frame = itt;
	    lib_buf2_dp++;		/* point to next buf2 */
	    lib_buf2_rp++;		/* point to next buf2 */
	    lib_bhead_p->hdr.buf2_cnt++;

	    if (lib_bhead_p->hdr.buf2_cnt >= max_do_cnt ||
		cur_buf_size >= max_buf_size) {
/* (2) send sequences for shuffling */
#if defined(COMP_THR) || defined(PCOMPLIB)	/* threaded - fill and empty buffers */
	      /* provide empty buffer to workers */
	      lib_bhead_p->hdr.aa1b_used = cur_buf_size;
	      lib_bhead_p->hdr.have_data = 1;
	      lib_bhead_p->hdr.seq_record_continuous = 0;
	      put_rbuf(lib_bhead_p,m_bufi_p->max_work_buf);
	      get_rbuf(&lib_bhead_p,m_bufi_p->max_work_buf);
#else		/* non-thread - just do the searches */
	      if (lib_bhead_p->hdr.buf2_type & BUF2_DOSHUF) {
		buf_shuf_work(aa0,m_msp->n0, aa1save, lib_bhead_p,
			      m_msp->nitt1, ppst, f_str, ppst->score_ix, rand_state);
	      }
#endif
/* (3) save results in the rstats structure */
	      if (lib_bhead_p->hdr.buf2_cnt > 0 && lib_bhead_p->hdr.have_results) {
		save_shuf(lib_bhead_p,m_msp->nitt1,m_msp->shuff_max,ppst->score_ix,s_info);
	      }

	      lib_bhead_p->s_cnt_info.s_cnt[0] = lib_bhead_p->s_cnt_info.s_cnt[1] =
		lib_bhead_p->s_cnt_info.s_cnt[2] = lib_bhead_p->s_cnt_info.tot_scores = 0;

	      lib_bhead_p->hdr.buf2_cnt = 0;
	      cur_buf_size = 0;
	      lib_bhead_p->hdr.have_results = 0;
	      lib_bhead_p->hdr.buf2_type=BUF2_DOSHUF;
	      lib_bhead_p->hdr.seq_record_continuous = 0; /* seq_records are coming from bestptr in any order */
	      lib_bhead_p->hdr.stop_work = 0;
	      lib_buf2_dp = lib_bhead_p->buf2_data;
	    }
	  } /* for (itt .. */
	}
      }			/* done with tmp_bestp[] */
      free(tmp_bestp);

#if defined(COMP_THR) || defined(PCOMPLIB)	/* if COMP_THR/PCOMPLIB - fill and empty buffers */
      /* check last buffers for any results */
      lib_bhead_p->hdr.seq_record_continuous = 0;
      put_rbuf(lib_bhead_p,m_bufi_p->max_work_buf);

      /* wait for the threads to finish */

      wait_rbuf(m_bufi_p->max_work_buf);
      /*
      fprintf(stderr, " num_reader[%d]-empty[%d]: %d\tnrstats: %d\n",
	      num_reader_bufs,empty_reader_bufs,
	      num_reader_bufs-empty_reader_bufs, nrstats);
      */

      for (i=0; i < num_reader_bufs; i++) {
	if (RESULTS_BUF[i]->hdr.buf2_cnt > 0 && RESULTS_BUF[i]->hdr.have_results) {
	  save_shuf(RESULTS_BUF[i],m_msp->nitt1, m_msp->shuff_max, ppst->score_ix, s_info);
	  RESULTS_BUF[i]->hdr.buf2_cnt = RESULTS_BUF[i]->hdr.have_results = 0;
	}
      }
#else	/* just do the searches */
      /* aa1save is used for shuffles, not aa1shuf, because aa1shuf
	 has library sequences */
      buf_shuf_work(aa0,m_msp->n0, aa1save, lib_bhead_p,
		    m_msp->nitt1, ppst, f_str, ppst->score_ix, rand_state);

      save_shuf(lib_bhead_p,m_msp->nitt1,m_msp->shuff_max, ppst->score_ix, s_info);
      lib_bhead_p->hdr.buf2_cnt = lib_bhead_p->hdr.have_results = 0;
#endif
}

/* buf_align_seq is structurally almost identical to buf_shuf_seq,
   except that the appropriate sequences are pre-loaded into bbp->seq
   (and ->bline), and it gets bbp->a_res, rather than scores */

void
buf_align_seq(unsigned char **aa0, int n0,
	      struct beststr **bestp_arr, int nbest,
	      struct pstruct *ppst, struct mngmsg *m_msp,
	      struct mng_thr *m_bufi_p
#if !defined(COMP_THR) && !defined(PCOMPLIB)
	      , void **f_str
#endif
	     )
{
  struct beststr *bbp;
  int max_align_cnt;
  int i, n_pre_align;
  int cur_buf_size, max_buf_size;
  struct buf_head *lib_bhead_p;
  struct buf2_data_s *lib_buf2_dp;
  struct buf2_ares_s *lib_buf2_ap;

  /* setup lib_bhead buffers for alignments */
#if defined(COMP_THR) || defined(PCOMPLIB)	/* threaded */
  /* max_do_cnt can be smaller than max_buf2_res, but not larger */
#ifdef COMP_THR
  max_align_cnt = min(m_bufi_p->max_buf2_res,
		      nbest / (4 * fa_max_workers));
#else
  max_align_cnt = min(m_bufi_p->max_buf2_res, nbest / fa_max_workers);
#endif
  if (max_align_cnt < 1) max_align_cnt = 1;

  get_rbuf(&lib_bhead_p,m_bufi_p->max_work_buf);
#else	/* not threaded */
  max_align_cnt = m_bufi_p->max_buf2_res;
  lib_bhead_p = lib_buf2_list;  /* equivalent to un-threaded get_rbuf() */
#endif

  max_buf_size = lib_bhead_p->hdr.aa1b_size;
  lib_bhead_p->hdr.buf2_cnt = 0;
  lib_bhead_p->hdr.have_results = 0;
  lib_bhead_p->hdr.stop_work = 0;
  lib_bhead_p->hdr.buf2_type=BUF2_DOALIGN;
  lib_buf2_dp = lib_bhead_p->buf2_data;
  lib_buf2_ap = lib_bhead_p->buf2_ares;

  /* read sequences into align buffer */

  n_pre_align = 0;
  cur_buf_size = 0;
  for (i = 0; i < nbest; i++) {
    bbp = bestp_arr[i];

    /* this invalidates lib_buf2_p->seq */
    lib_buf2_dp->seq = bbp->seq;
    cur_buf_size += bbp->seq->n1+1;
    lib_buf2_dp->frame = bbp->frame;
    lib_buf2_dp->repeat_thresh = bbp->repeat_thresh;
#ifdef PCOMPLIB
    lib_buf2_dp->seq_dup = 0;
#endif
    lib_buf2_ap->have_ares = 0;
    lib_buf2_ap->a_res = NULL;
    lib_buf2_ap->best_idx = i;
    lib_buf2_dp++;		/* point to next buf2_data */
    lib_buf2_ap++;		/* point to next buf2_ares */
    lib_bhead_p->hdr.buf2_cnt++;

    if (lib_bhead_p->hdr.buf2_cnt >= max_align_cnt ||
	cur_buf_size >= max_buf_size - m_msp->ldb_info.maxn) {
/* (2) send sequences for alignment */
#if defined(COMP_THR) || defined(PCOMPLIB)	/* threaded - fill and empty buffers */
      /* provide empty buffer to workers */
      lib_bhead_p->hdr.seqr_cnt = lib_bhead_p->hdr.buf2_cnt;	/* for alignments, they are the same */
      lib_bhead_p->hdr.have_data = 1;
      lib_bhead_p->hdr.aa1b_used = cur_buf_size;
      lib_bhead_p->hdr.seq_record_continuous = 0;
      put_rbuf(lib_bhead_p,m_bufi_p->max_work_buf);
      get_rbuf(&lib_bhead_p,m_bufi_p->max_work_buf);
#else		/* non-thread - just do the searches */
      buf_do_align(aa0, m_msp->n0, lib_bhead_p, ppst, m_msp, f_str);
#endif

/* (3) save alignments */
      if (lib_bhead_p->hdr.buf2_cnt > 0 && lib_bhead_p->hdr.have_results) {
	n_pre_align += save_align(lib_bhead_p,bestp_arr);
      }

      cur_buf_size = 0;
      max_buf_size = lib_bhead_p->hdr.aa1b_size;
      lib_bhead_p->hdr.buf2_cnt = 0;
      lib_bhead_p->hdr.have_results = 0;
      lib_bhead_p->hdr.buf2_type=BUF2_DOALIGN;
      lib_bhead_p->hdr.stop_work = 0;
      lib_buf2_dp = lib_bhead_p->buf2_data;
      lib_buf2_ap = lib_bhead_p->buf2_ares;
    }
  }			/* done with bestp_arr[] */

#if defined(COMP_THR) || defined(PCOMPLIB)	/* if COMP_THR - fill and empty buffers */
  /* check last buffers for any results */
  lib_bhead_p->hdr.seqr_cnt = lib_bhead_p->hdr.buf2_cnt;	/* for alignments, they are the same */
  lib_bhead_p->hdr.have_data = 1;
  lib_bhead_p->hdr.aa1b_used = cur_buf_size;
  lib_bhead_p->hdr.seq_record_continuous = 0;
  put_rbuf(lib_bhead_p,m_bufi_p->max_work_buf);

  /* wait for the threads to finish */

  wait_rbuf(m_bufi_p->max_work_buf);

  for (i=0; i < num_reader_bufs; i++) {
    if (RESULTS_BUF[i]->hdr.buf2_cnt > 0 && RESULTS_BUF[i]->hdr.have_results) {
      n_pre_align += save_align(RESULTS_BUF[i],bestp_arr);
      RESULTS_BUF[i]->hdr.buf2_cnt = RESULTS_BUF[i]->hdr.have_results = 0;
    }
  }
#else	/* just do the searches */
  buf_do_align(aa0, m_msp->n0, lib_bhead_p, ppst, m_msp, f_str);
  n_pre_align += save_align(lib_bhead_p,bestp_arr);
  lib_bhead_p->hdr.buf2_cnt = lib_bhead_p->hdr.have_results = 0;
#endif

  m_msp->align_done = 1;

  if (n_pre_align != nbest) {
    fprintf(stderr,"*** error [%s:%d] -  n_pre_align:%d != nbest: %d\n",
	    __FILE__, __LINE__, n_pre_align, nbest);
  }
  for (i=0; i < nbest; i++) {
    if (bestp_arr[i]->a_res == NULL) {
      fprintf(stderr, "*** error [%s:%d] - have NULL a_res: %d\n",
	      __FILE__, __LINE__, i);
    }
  }
}

int
check_seq_range(unsigned char *aa1b, int n1, int nsq, char *str) {
  int i, range_error;
  unsigned char *aa1p;

  range_error = 0;
  for (aa1p = aa1b, i=0; i < n1; i++, aa1p++) {
    if (*aa1p > nsq) {
      range_error = 1;
      /*      fprintf(stderr, "%s seq %d (%c) out of range at %d\n",
	      str, *aa1p, *aa1p,i);
      */
    }
  }
  return range_error;
}

struct stack_str {
  void **stack;
  int size;
  int inc;
  int top;
};

struct stack_str *init_stack(int size, int inc) {
  struct stack_str *stack;

  if ((stack=(struct stack_str *)calloc(1,sizeof(struct stack_str)))==NULL) {
    fprintf(stderr,"*** error [%s:%d] - cannot allocate stack\n",
	    __FILE__, __LINE__);
    return NULL;
  }

  if ((stack->stack=(void *)calloc(size,sizeof(void *)))==NULL) {
    fprintf(stderr,"*** error [%s:%d] - cannot allocate stack->stack[%d]\n",
	    __FILE__, __LINE__,size);
    free(stack);
    return NULL;
  }

  stack->size = size;
  stack->inc = inc;
  stack->top = 0;
  return stack;
}

void push_stack(struct stack_str *stack, void *value) {

  if (!stack) return;
  if (stack->top >= stack->size) {
    stack->size += stack->inc;
    if ((stack->stack = (void *)realloc(stack->stack, stack->size*sizeof(void *)))==NULL) {
      fprintf(stderr,"*** error [%s:%d] - cannot re-allocate stack to [%d]\n",
	      __FILE__, __LINE__, stack->size);
      return;
    }
  }
  stack->stack[stack->top++] = value;
}

void * pop_stack(struct stack_str *stack) {
  if (stack == NULL) {
#ifdef DEBUG
    fprintf(stderr," *** error [%s:%d] - pop_stack NULL stack\n",__FILE__, __LINE__);
#endif
    return NULL;
  }

  if (stack->top-- > 0) {
    return stack->stack[stack->top];
  }
  else {
    stack->top = 0;
    return NULL;
  }
}

void * free_stack(struct stack_str *stack) {
  if (stack==NULL) return NULL;
  if (stack->stack != NULL) free(stack->stack);
  free(stack);
  return NULL;
}

struct dyn_string_str *
init_dyn_string(int size, int inc) {
  struct dyn_string_str *dyn_string;

  if ((dyn_string=(struct dyn_string_str *)calloc(1,sizeof(struct dyn_string_str)))==NULL) {
    fprintf(stderr,"*** error [%s:%d] - cannot allocate dyn_string\n",
	    __FILE__, __LINE__);
    return NULL;
  }

  if ((dyn_string->string=(void *)calloc(size,sizeof(void *)))==NULL) {
    fprintf(stderr,"*** error [%s:%d] - cannot allocate dyn_string->string[%d]\n",
	    __FILE__, __LINE__,size);
    free(dyn_string);
    return NULL;
  }

  dyn_string->c_size = 0;
  dyn_string->inc = inc;
  dyn_string->mx_size = size;
  return dyn_string;
}

void
dyn_strcat(struct dyn_string_str *dyn_string, char *value) {
  size_t add_len;

  add_len = strlen(value);

  if (!dyn_string) return;
  if (add_len + dyn_string->c_size + 1 >= dyn_string->mx_size) {
    while (dyn_string->inc < add_len) { dyn_string->inc *= 2; }
    dyn_string->mx_size += dyn_string->inc;
    if ((dyn_string->string = (void *)realloc(dyn_string->string, dyn_string->mx_size))==NULL) {
      fprintf(stderr,"*** error [%s:%d] - cannot re-allocate dyn_string to [%d]\n",
	      __FILE__, __LINE__, dyn_string->mx_size);
      dyn_string->mx_size = 0;
      return;
    }
  }
  SAFE_STRNCAT(dyn_string->string,value,dyn_string->mx_size);
  dyn_string->c_size += add_len;
}

void dyn_strcpy(struct dyn_string_str *dyn_string, char *value) {
  size_t add_len;

  add_len = strlen(value);

  if (!dyn_string) return;
  if (add_len + 1>= dyn_string->mx_size) {
    while (dyn_string->inc < add_len) { dyn_string->inc *= 2; }
    dyn_string->mx_size += dyn_string->inc;
    if ((dyn_string->string = (void *)realloc(dyn_string->string, dyn_string->mx_size))==NULL) {
      fprintf(stderr,"*** error [%s:%d] - cannot re-allocate dyn_string to [%d]\n",
	      __FILE__, __LINE__, dyn_string->mx_size);
      dyn_string->mx_size = 0;
      return;
    }
  }
  SAFE_STRNCPY(dyn_string->string,value,dyn_string->mx_size);
}

void free_dyn_string(struct dyn_string_str *dyn_string) {
  if (dyn_string==NULL) return;
  if (dyn_string->string != NULL) free(dyn_string->string);
  free(dyn_string);
}

#include "a_mark.h"

/* *itmp has the current alignment score, if *annot_arr[i_annot].label='V',
     this can be increased (total increase in *v_delta)
   *pam2aa0v[.value] gives possibly better pam score for variant
   *ip is position in annotated sequence (&i0 for annot0_p)
   *ia is position in aligned sequence (&i1 for annot0_p)
   sp1 is the array for the (possibly modified) displayed sequence
   sp1a is the array for the associated annotation
   sq maps encoded residues to displayed characters
   i_annot -- current annotation index in annot0_p->annot_arr_p[i_annot]
   annot_arr = annot0/1_p->annot_arr_p
   annot_stack = save current annotation
   *have_push_features = set for annotations pushed in stack (not 'V')
   *v_delta = change in score from variant at this position
   **region_p = set for '[' region start
   init_score -- used to initialize tmp_region_p->score.

*/

int
next_annot_match(int *itmp, int *pam2aa0v, long ip, long ia, char *sp1, char *sp1a, const unsigned char *sq,
		 int i_annot, int n_annot, struct annot_entry **annot_arr, char **ann_comment,
		 void *annot_stack, int *have_push_features, int *v_delta,
		 struct annot_entry **region_p, struct annot_entry *tmp_region_p,
		 int init_score)  {
  int v_tmp;

  if (ann_comment) *ann_comment = NULL;

  /* count through the annotations at this position (long ip) */

  while (i_annot < n_annot && ip == annot_arr[i_annot]->pos) {
    if (annot_arr[i_annot]->label == 'V') { /* label == 'V' */
      v_tmp = pam2aa0v[annot_arr[i_annot]->value];
      if (v_tmp > *itmp) {
	*v_delta += (v_tmp- *itmp);
	*itmp = v_tmp;
	*sp1 = sq[annot_arr[i_annot]->value];
	if (sp1a) *sp1a = 'V';
	if (ann_comment) *ann_comment = annot_arr[i_annot]->comment;
      }
    }
    else if (annot_arr[i_annot]->label == '[') {
      /* region_p needs to point to a more sophisticated data
	 structure that keeps track of all the current regions being
	 updated

	 to start, region_p could include a linked list and a pointer to
	 the current left-most region, which would be used for ']'
	 detection

	 for efficiency, update the ->score only when a new
	 (overlapping) region is started or stopped

	 same for n_indent, n_aln
      */

      if (region_p) {
	memcpy(tmp_region_p, annot_arr[i_annot],sizeof(struct annot_entry));
        tmp_region_p->a_pos = ia;
	tmp_region_p->score = init_score;
	tmp_region_p->n_ident = tmp_region_p->n_aln = 0;
	*region_p = tmp_region_p;
      }
    }
    else if (annot_arr[i_annot]->label == ']') {
      if (have_push_features) *have_push_features = 1;
      push_stack(annot_stack, annot_arr[i_annot]);
    }
    else if (annot_stack) {
      if (have_push_features) *have_push_features = 1;
      push_stack(annot_stack, annot_arr[i_annot]);
    }
    i_annot++;
  }  /* everything at this alignment position is checked */
  return i_annot;
}

/* returns M_NEG, M_ZERO, M_POS, M_IDENT, M_DEL (a_mark.h)
   updates *aln->nsim, npos, nident, nmismatch

*/
int align_type(int score, char sp0, char sp1, int nt_align, struct a_struct *aln, int pam_x_id_sim) {
  int spa_val;

  if (score<0) {
    spa_val = M_NEG;
  }
  else if (score == 0) {
    spa_val = M_ZERO;
    if (aln) aln->nsim++;
  }
  else {
    spa_val = M_POS;
    if (aln) {aln->nsim++; aln->npos++;}
  }

  /* correct for score < 0 with 'N:N'/'X:X' */
  if (pam_x_id_sim > 0) {	/* > 0 -> identical, similar */
    if ((nt_align && toupper(sp0)=='N' && toupper(sp1)=='N') ||
	(!nt_align && toupper(sp0)=='X' && toupper(sp1)=='X')) {
      spa_val = M_POS;
      if (aln) {
	aln->nsim++;
      }
    }
  }

  if (aln) aln->nmismatch++;
  if (toupper(sp0) == toupper(sp1)) {
    spa_val = M_IDENT;
    if (aln) {
      aln->nident++;
      aln->nmismatch--;
    }
  }
  else if (nt_align) {
    if ((toupper(sp0) == 'T' && toupper(sp1) == 'U') ||
	(toupper(sp0)=='U' && toupper(sp1)=='T')) {
      spa_val = M_IDENT;
      if (aln) {
	aln->nident++;
	aln->nmismatch--;
      }
    }
    /* add to gap count for 'N' matches ?? */
    else if (aln && toupper(sp0) == 'N') aln->ngap_q++;
    else if (aln && toupper(sp1) == 'N') aln->ngap_l++;
  }

  /* correct nident, nmismatch for N:N / X:X */
  if (pam_x_id_sim < 0) {	/* > 0 -> identical, similar */
    if ((nt_align && toupper(sp0)=='N' && toupper(sp1)=='N') ||
	(!nt_align && toupper(sp0)=='X' && toupper(sp1)=='X')) {
      if (aln) {
	aln->nident--;
	aln->nmismatch++;
      }
    }
  }

  return spa_val;
}

/* seq_pos works with comment_var()/display_push_features()/do_url1() where
   i_offset = nn for reversed sequences
   off = 0 for 0 based offsets, 1 for 1-based offsets
 */
int
seq_pos(int pos, int rev, int off) {

  if (rev) {
    return -pos-1 + off;
  }
  else {
    return pos;
  }
}

/* target = 0 (aa0), 1 (aa1)

   d_type = display_type (annot_fmt in cal_cons.c):
            1 (long text),   d1_fmt = " Variant: %d%c%c%d%c : %c%d%c";
            2 (-m 9c code)   sprintf(tmp_str, "|%c%c:%ld%c%c%ld%c",

   i0_pos/i1_pos have already been converted to reverse coordinate if necessary
*/
void comment_var (long i0_pos, char sp0, long i1_pos, char sp1, char o_sp1,
		  char sim_char, const char *ann_comment,
		  struct dyn_string_str *annot_var_dyn, int target, int d_type)
{
  char tmp_str[MAX_LSTR], tc, ann_ch0, ann_ch1;
  char *d1_fmt;

  if (d_type == 1) {
    if (target ==1) {
      d1_fmt = " Variant: %d%c%c%d%c : %c%d%c";
      sprintf(tmp_str,d1_fmt,
	      i0_pos+1, sp0, sim_char, i1_pos+1,sp1, o_sp1,i1_pos+1,sp1);
    }
    else {
      d1_fmt = " qVariant: %d%c%c%d%c : %c%d%c";
      sprintf(tmp_str,d1_fmt,
	      i0_pos+1, sp0, sim_char, i1_pos+1,sp1, o_sp1,i1_pos+1,sp0);
    }

    /* SAFE_STRNCAT(annot_var_s,tmp_str,n_annot_var_s); */
    dyn_strcat(annot_var_dyn, tmp_str);

    if (ann_comment) {
      sprintf(tmp_str," : %s",ann_comment);
      /* SAFE_STRNCAT(annot_var_s,tmp_str,n_annot_var_s); */
      dyn_strcat(annot_var_dyn, tmp_str);
    }

    /* SAFE_STRNCAT(annot_var_s,"\n",n_annot_var_s); */
    dyn_strcat(annot_var_dyn, "\n");
  }
  else if (d_type == 2) {
    if (target == 1) {
      ann_ch0 = 'X';
      ann_ch1 = 'V';
    }
    else {
      ann_ch0 = 'V';
      ann_ch1 = 'X';
    }

    sprintf(tmp_str, "|%c%c:%ld%c%c%ld%c",
	    ann_ch0,ann_ch1,
	    i0_pos+1,sp0, sim_char,i1_pos+1,sp1);
    /* SAFE_STRNCAT(annot_var_s, tmp_str, n_annot_var_s); */
    dyn_strcat(annot_var_dyn, tmp_str);
  }
}

void
display_push_features(void *annot_stack, struct dyn_string_str *annot_var_dyn,
		      long i0_pos, char sp0, long i1_pos, char sp1, char sym,
		      struct annot_entry **region0_p,
		      struct annot_entry **region1_p,
		      int tot_score, double comp, int n0, int n1,
		      void *pstat_void, int d_type) {
  struct annot_entry *this_annot_p;
  double lbits, total_bits, zscore, lprob, lpercid;
  char *ann_comment, *bp;
  struct annot_entry *region_p;
  char tmp_lstr[MAX_LSTR], ctarget, tmp_sstr[MAX_SSTR];
  int q_min, q_max, l_min, l_max;
  char *dt1_fmt, *dt2_fmt;

  zscore = find_z(tot_score, 1.0, n1, comp, pstat_void);
  total_bits = zs_to_bit(zscore, n0, n1);

  while ((this_annot_p = (struct annot_entry *)pop_stack(annot_stack))!=NULL) {

    if (this_annot_p->label == ']') {
      if (this_annot_p->target == 1) {
	region_p = *region1_p;
	if (!region_p) {
	  fprintf(stderr,"*** error [%s:%d] *** -- target==1 but region1_p is null\n",__FILE__, __LINE__);
#ifdef DEBUG
	  fprintf(stderr,"*** qtitle: %s\n",ext_qtitle);
#endif
	  continue;
	}
	q_min = region_p->a_pos+1;
	l_min = region_p->pos+1;
	dt2_fmt = "|XR:%d-%d:%d-%d:s=%d;b=%.1f;I=%.3f;Q=%.1f";
      }
      else {
	region_p = *region0_p;
	if (!region_p) {
	  fprintf(stderr,"*** error [%s:%d] *** -- target==0 but region0_p is null\n",__FILE__, __LINE__);
#ifdef DEBUG
	  fprintf(stderr,"*** qtitle: %s\n",ext_qtitle);
#endif
	  continue;
	}
	q_min = region_p->pos+1;
	l_min = region_p->a_pos+1;
	dt2_fmt = "|RX:%d-%d:%d-%d:s=%d;b=%.1f;I=%.3f;Q=%.1f";
      }

      if (region_p->score < 0) {
	lbits = 0.0;
	lprob = 1.0;
      }
      else {
	lbits = total_bits * (double)region_p->score/tot_score;
	zscore = find_z(region_p->score, 1.0, n1, comp, pstat_void);
	lprob = zs_to_p(zscore);
      }

      if (lprob > 0.99) lprob = 0.0;
      else if (lprob < 1e-300) lprob = 3000.0;
      else lprob = -10.0*log(lprob)/log(10.0);

      if (region_p->n_aln > 0) {
	lpercid = ((double)region_p->n_ident)/(double)region_p->n_aln;
      }
      else lpercid = -1.0;

      if (d_type == 1) {
	if (this_annot_p->target == 0) {dt1_fmt = " qRegion: %d-%d:%d-%d : score=%d; bits=%.1f; Id=%.3f; Q=%.1f :  %s\n";}
	else {dt1_fmt = " Region: %d-%d:%d-%d : score=%d; bits=%.1f; Id=%.3f; Q=%.1f :  %s\n";}
	sprintf(tmp_lstr, dt1_fmt, q_min, i0_pos+1,
		l_min, i1_pos+1, region_p->score, lbits, lpercid, lprob,
		(region_p->comment) ? region_p->comment : '\0');

      }
      else if (d_type == 2) {
	sprintf(tmp_lstr,dt2_fmt,
		q_min, i0_pos+1,
		l_min, i1_pos+1, region_p->score, lbits,lpercid, lprob);

	if (region_p->comment) {
	  SAFE_STRNCPY(tmp_sstr,region_p->comment,sizeof(tmp_sstr));
	  if ((bp=strchr(tmp_sstr,' '))!=NULL) { *bp = '\0';}
	  SAFE_STRNCAT(tmp_lstr,";C=",sizeof(tmp_lstr));
	  SAFE_STRNCAT(tmp_lstr,tmp_sstr,sizeof(tmp_lstr));
	}
      }
      /* SAFE_STRNCAT(annot_var_s,tmp_lstr,n_annot_var_s); */
      dyn_strcat(annot_var_dyn, tmp_lstr);
      region_p->score = 0;
      region_p = NULL;
    }
    else if ((ann_comment = this_annot_p->comment)) {
      if (d_type == 1 ) {
	if (this_annot_p->target == 0) {dt1_fmt = " qSite:%c : %d%c%c%d%c : %s\n";}
	else {dt1_fmt = " Site:%c : %d%c%c%d%c : %s\n";}
	sprintf(tmp_lstr,dt1_fmt, this_annot_p->label,i0_pos+1, sp0,
		sym, i1_pos+1, sp1, ann_comment);
	/* SAFE_STRNCAT(annot_var_s,tmp_lstr,n_annot_var_s); */
	dyn_strcat(annot_var_dyn, tmp_lstr);
      }
    }
  }
}
