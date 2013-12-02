
/* copyright (c) 1996, 1997, 1998, 1999 William R. Pearson and the
   U. of Virginia */

/*  $Id: compacc.c 721 2011-05-06 17:45:58Z wrp $ */
/* $Revision: 721 $  */

/* Concurrent read version */

#include <stdio.h>
#include <stdlib.h>
#if defined(UNIX) || defined(WIN32)
#include <sys/types.h>
#endif

#include <limits.h>
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

extern time_t tdone, tstart;		/* Timing */
extern void abort ();
extern void ptime ();

#include "drop_func.h"	/* get init_work() */
void revcomp(unsigned char *seq, int n, int *c_nt);
extern void qshuffle(unsigned char *aa0, int n0, int nm0, void *);
#ifdef DEBUG
unsigned long adler32(unsigned long, const unsigned char *, unsigned int);
#endif


/* this function consolidates code in comp_lib4.c for non-threaded, and in
   work_thr2.c (threads) and work_comp2.c (worker nodes)
*/

void init_aa0(unsigned char **aa0, int n0, int nm0,
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
      fprintf(stderr," cannot allocate aa01[%d]\n", n0);
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
      fprintf(stderr,"cannot allocate aa0s[%d]\n",n0+2);
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
    fprintf(stderr,"unable to allocate shuffled library sequence [%d]\n", max_tot);
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

    if (m_msp->n0 > m_msp->max_tot/3) {
      fprintf(stderr," query sequence is too long %d > %d %s\n",
	      m_msp->n0,
	      m_msp->max_tot/3,
	      m_msp->sqnam);
      exit(1);
    }

    m_msp->ldb_info.l_overlap = over_len;
    m_msp->ldb_info.maxt3 = maxn-m_msp->ldb_info.l_overlap;
  }
  else {	/* is TFAST */
    if (m_msp->n0 > MAXTST) {
      fprintf(stderr," query sequence is too long %d %s\n",m_msp->n0,m_msp->sqnam);
      exit(1);
    }

    if (m_msp->n0*3 > maxn ) {	/* n0*3 for the three frames - this
				   will only happen if maxn has been
				   set low manually */

      if (m_msp->n0*4+2 < m_msp->max_tot) { /* m_msg0*3 + m_msg0 */
	fprintf(stderr,
		" query sequence too long for library segment: %d - resetting to %d\n",
	      maxn,m_msp->n0*3);
	maxn = m_msp->ldb_info.maxn = m_msp->n0*3;
      }
      else {
	fprintf(stderr," query sequence too long for translated search: %d * 4 > %d %s\n",
	      m_msp->n0,maxn, m_msp->sqnam);
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
  for (i=0; (size_t)i < strlen(str); i++) aaray[qascii[str[i]]]=1;
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
	  fprintf(stderr, "Unterminated ENV: %s\n",src);
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


void sortbeste(struct beststr **bptr, int nbest)
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

extern double zs_to_Ec(double zs, long entries);

#include "aln_structs.h"

/*
extern double ks_dev;
extern int ks_df; */

void
prhist(FILE *fd, const struct mngmsg *m_msp,
       struct pstruct *ppst, 
       struct hist_str hist,
       int nstats, int sstats,
       struct db_str ntt,
       char *stat_info2,
       char *lib_range,
       char **info_gstring2,
       char **info_hstring)
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


  fprintf(fd,"\n");
  
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
		struct mngmsg *m_msp, struct pstruct *ppst,
		void *f_str
		);

extern struct lmf_str *
re_openlib(struct lmf_str *, int outtty);

#define MAX_BLINE 2048
#define RANLIB (m_fptr->ranlib)

extern int
re_getlib(unsigned char *, unsigned char **, 
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
	      struct mngmsg *m_msp)
{
  int i, n1, bl_len, tmp_bline_len, l_llen;
  int seq_buf_len;
  char bline[MAX_BLINE];
  unsigned char  *seq_buf_p;
  char *bline_buf_p;

  struct beststr *bbp;
  struct rstruct rst;
  struct lmf_str *m_fptr;

  /* 
     calculate how much room we need for sequences and blines 
  */
  
  seq_buf_len = 1;
  for (i=0; i<nbest; i++) {
    /* we are not (currently) allocating more than n1+1, because alignment is not vectorized,
       if it were vectorized, we would need n+16
    */
    seq_buf_len += bbp_arr[i]->seq->n1 + 1;
  }

  if ((m_msp->aa1save_buf_b=(unsigned char *)calloc(seq_buf_len, sizeof(char)))==NULL) {
    fprintf(stderr, "*** error - cannot allocate space[%d] for sequence encoding\n",seq_buf_len);
    exit(1);
  }
  else {
    seq_buf_p = m_msp->aa1save_buf_b+1;		/* ensure there is an initial '\0' */
  }
  
  l_llen = m_msp->aln.llen;
  if ((m_msp->markx & MX_M9SUMM) && m_msp->show_code != SHOW_CODE_ID) {
    l_llen += 40;
    if (l_llen > 200) l_llen=200;
  }

  tmp_bline_len = sizeof(bline)-1;
  if (!(m_msp->markx & MX_M10FORM) && !m_msp->long_info) {tmp_bline_len = l_llen-5;}

  /* allocate more bline than we need for simplicity */
  if ((bline_buf_p=m_msp->bline_buf_b=(char *)calloc(nbest*tmp_bline_len, sizeof(char)))==NULL) {
    fprintf(stderr, "*** error - cannot allocate space[%d] for bline descriptions\n",nbest*tmp_bline_len);
    exit(1);
  }

  for (i=0; i<nbest; i++) {
    bbp = bbp_arr[i];

    if ((m_fptr=re_openlib(bbp->mseq->m_file_p,!m_msp->quiet))==NULL) {
      fprintf(stderr,"*** cannot re-open %s\n",bbp->mseq->m_file_p->lb_name);
      exit(1);
    }
    RANLIB(bline,tmp_bline_len,bbp->mseq->lseek,bbp->mseq->libstr,m_fptr);
    bl_len = strlen(bline);
    bbp->mseq->bline = bline_buf_p;
    bbp->mseq->bline_max = m_msp->aln.llen;
    strncpy(bbp->mseq->bline, bline, bl_len);
    bline_buf_p += bl_len+1;


    /* make sure we get annotation if present, and sequence if necessary */
    if (bbp->seq->aa1b==NULL || (m_msp->ann_flg && bbp->seq->aa1_ann==NULL)) {
      n1 = re_getlib(aa1save, m_msp->ann_flg ? &(bbp->seq->aa1_ann) : NULL, 
		     maxn,m_msp->ldb_info.maxt3, m_msp->ldb_info.l_overlap,
		     bbp->mseq->cont,m_msp->ldb_info.term_code,
		     &bbp->seq->l_offset,&bbp->seq->l_off,bbp->mseq->m_file_p);
      if (n1 != bbp->seq->n1) {
	fprintf(stderr," *** error n1[%d] != n1[%d] from re_getlib()\n", 
		bbp->n1, n1);
      }

#ifdef DEBUG
      if (adler32(1L,aa1save,n1)!=bbp->adler32_crc) {
	fprintf(stderr," *** error adler32_crc from re_getlib()\n");
      }
#endif

      if (bbp->seq->aa1b == NULL)  {
	bbp->seq->aa1b = seq_buf_p;
	memcpy(bbp->seq->aa1b, aa1save, bbp->seq->n1+1);
	seq_buf_p += bbp->seq->n1+1;
      }
    }
  }
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

#define ESS 49

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
