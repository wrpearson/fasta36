/*  $Id: work_thr2.c $ */

/* copyright (c) 1996, 1997, 1998, 1999, 2014 by William R. Pearson
   and The The Rector & Visitors of the University of Virginia */

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

/* work_thr.c - threaded worker */

/* modified 21-Oct-1998 to work with reverse complement for DNA */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <sys/types.h>
#include <signal.h>

#include "defs.h"		/* various constants */
#include "best_stats.h"			/* defines beststr */
#include "structs.h"
#include "param.h"		/* pstruct rstruct */
#include "thr_buf_structs.h"

/***************************************/
/* thread global variable declarations */
/***************************************/

#ifndef PCOMPLIB
#define XTERNAL
#include "thr_bufs2.h"
#undef XTERNAL
#else
#include "msg.h"
#define XTERNAL
#include "uascii.h"
#undef XTERNAL
#ifdef MPI_SRC
#include "mpi.h"
#endif
#endif

void alloc_pam (int, int, struct pstruct *);
int **alloc_pam2p(int **,int, int);
void revcomp(unsigned char *seq, int n, int *c_nt);

#if defined(WIN32) || !defined(THR_EXIT)
void pthread_exit(void *);
#define THR_EXIT pthread_exit
#else
void THR_EXIT(void *);
#endif

#ifdef DEBUG
extern struct buf_head *lib_buf2_list;
#endif

/* functions getting/sending buffers to threads (thr_sub.c) */
extern void wait_thr(void);
extern int get_wbuf(struct buf_head **cur_buf, int max_work_buf);
extern void put_wbuf(struct buf_head *cur_buf, int max_work_buf);

/* dropxx.c functions */
#include "drop_func.h"

extern void *my_srand();
extern unsigned int my_nrand(int, void *);
extern void qshuffle(unsigned char *aa0, int n0, int nm0, void *);
extern void free_pam2p(int **);

void init_aa0(unsigned char **aa0, int n0, int nm0,
	      unsigned char **aa0s, unsigned char **aa1s, 
	      int qframe, int qshuffle_flg, int max_tot,
	      struct pstruct *ppst, void **f_str, void **qf_str,
	      void *my_rand_state);

extern void
buf_do_work(unsigned char **aa0, int n0, struct buf_head *lib_bhead_p,
	    int max_frame, struct pstruct *ppst, void **f_str);
extern void
buf_qshuf_work(unsigned char *aa0s, int n0, struct buf_head *lib_bhead_p,
	       int max_frame, struct pstruct *ppst, void *qf_str, int score_ix);
extern void
buf_shuf_work(unsigned char **aa0, int n0, unsigned char *aa1s,
	      struct buf_head *lib_bhead_p, int max_frame, struct pstruct *ppst,
	      void **f_str, int score_ix, void *);

void
buf_do_align(unsigned char **aa0,  int n0,
	     struct buf_head *lib_bhead_p, 
	     struct pstruct *ppst, const struct mngmsg *my_msp,
	     void **f_str);

#ifndef PCOMPLIB
#define FIRSTNODE 0
void
work_thread (struct thr_str *work_info)
#else
#if defined(TFAST)
extern void aainit(int tr_type, int debug);
#endif

int g_worker;

void work_comp(int my_worker)
#endif
{
  struct buf_head *cur_buf, *my_cur_buf;
  char info_lib_range[MAX_FN];
  unsigned char *aa1s=NULL;
#ifndef PCOMPLIB

  const struct mngmsg *my_msp;
  int my_worker;
#else
#ifdef MPI_SRC
  struct mngmsg *my_msp;
  MPI_Status mpi_status;
  int buf_alloc_flag = 0;
#endif
  struct mngmsg my_msg;
  int int_msg_b[4];
  struct buf2_data_s *my_buf2_data;
  struct buf2_res_s *my_buf2_res;
  struct buf2_ares_s *my_buf2_ares;
  struct seq_record *my_seq_buf;
  unsigned char *my_aa1b_buf;
#endif
  int i, j, npam, n0, nm0;
  int max_work_buf, max_buf2_res, max_chain_seqs, seq_buf_size;
  void *my_rand_state;

  struct pstruct my_pst, *my_ppst;
  unsigned char *aa0[6], *aa0s;
  void *f_str[6], *qf_str;

  my_rand_state=my_srand();

#ifndef PCOMPLIB
  my_worker = work_info->worker;
  max_work_buf = work_info->max_work_buf;
  wait_thr();	/* wait for start_thread predicate to drop to  0 */

  my_msp = work_info->m_msp;
#else 	/* PCOMPLIB */

#ifdef DEBUG
/*  fprintf(stderr,"%d: work_comp started\n",my_worker); */
#endif
  g_worker = my_worker;
  my_msp = &my_msg;

#ifdef MPI_SRC
 pcomp_loop:

  MPI_Recv(int_msg_b,4,MPI_INT,0, STARTTYPE0,MPI_COMM_WORLD,
	   &mpi_status);
  
  max_work_buf = int_msg_b[0];
  max_buf2_res = int_msg_b[1];
  max_chain_seqs = int_msg_b[2];
  seq_buf_size = int_msg_b[3];

  /* quit the main loop with a message of 0 max_work_buf */
  if (max_work_buf == 0) { goto pcomp_final;}

  MPI_Recv((void *)my_msp,sizeof(struct mngmsg),MPI_BYTE,0,STARTTYPE1,MPI_COMM_WORLD,
	   &mpi_status);

  MPI_Recv((void *)&my_pst,(int)sizeof(struct pstruct),MPI_BYTE,0,STARTTYPE2,MPI_COMM_WORLD,
	   &mpi_status);
  my_ppst = &my_pst;

#endif	/* MPI_SRC */

  if (!buf_alloc_flag) {
    buf_alloc_flag = 1;
    /* must allocate buffers for data, sequences, results */
    if ((my_cur_buf = cur_buf = (struct buf_head *)calloc(1,sizeof(struct buf_head)))==NULL) {
      fprintf(stderr,"cannot allocate buf_head\n");
      exit(1);
    }

    /* allocate results array */
    if ((my_buf2_res = (struct buf2_res_s*)calloc(max_buf2_res+1,sizeof(struct buf2_res_s)))==NULL) {
      fprintf(stderr,"cannot allocate buf2_data[%d]\n",max_buf2_res);
      exit(1);
    }
    cur_buf->buf2_res = my_buf2_res;

    /* allocate buffers for ares alignment encodings */
    if ((my_buf2_ares = (struct buf2_ares_s*)calloc(max_buf2_res+1,sizeof(struct buf2_ares_s)))==NULL) {
      fprintf(stderr,"cannot allocate buf2_data[%d]\n",max_buf2_res);
      exit(1);
    }
    cur_buf->buf2_ares = my_buf2_ares;

    /* allocate buffers for data */
    if ((my_buf2_data = (struct buf2_data_s*)calloc(max_buf2_res+1,sizeof(struct buf2_data_s)))==NULL) {
      fprintf(stderr,"cannot allocate buf2_data[%d]\n",max_buf2_res);
      exit(1);
    }
    cur_buf->buf2_data = my_buf2_data;

    /* also must allocate seq_records */
    if ((my_seq_buf = 
	 (struct seq_record *)calloc((size_t)(max_buf2_res+1), sizeof(struct seq_record)))
        ==NULL) {
      fprintf(stderr,"%d: cannot allocate seq_record buffer[%d]\n",my_worker,max_buf2_res+1);
      exit(1);
    }
    cur_buf->buf2_data[0].seq = cur_buf->hdr.seq_b = my_seq_buf;

    if ((my_aa1b_buf = (unsigned char *)calloc((size_t)(seq_buf_size+1),sizeof(unsigned char)))
        ==NULL) {
      fprintf(stderr,"%d: cannot allocate sequence buffer[%d]\n",my_worker, seq_buf_size);
      exit(1);
    }
    else {	  /* now associate the my_aa1b_buf with cur_buf */
      my_aa1b_buf++;
      cur_buf->hdr.aa1b_start = cur_buf->buf2_data[0].seq->aa1b = my_aa1b_buf;
      cur_buf->hdr.aa1b_size = seq_buf_size;
   }
  }
  else {
    cur_buf = my_cur_buf;
    cur_buf->buf2_data = my_buf2_data;
    cur_buf->buf2_data[0].seq = cur_buf->hdr.seq_b = my_seq_buf;
    cur_buf->buf2_res = my_buf2_res;
    cur_buf->buf2_ares = my_buf2_ares;
    cur_buf->hdr.aa1b_start = cur_buf->buf2_data[0].seq->aa1b = my_aa1b_buf;
    cur_buf->hdr.aa1b_size = seq_buf_size;
  }

#if defined(TFAST)
    /* set up translation tables: faatran.c */
  aainit(my_ppst->tr_type,my_ppst->debug_lib);
#endif

#endif	/* PCOMPLIB */

  /* the pam allocation stuff is very different for threaded vs PCOMPLIB,
     so the code is separate */
#if !defined(PCOMPLIB)
  /* make certain that all but 0 have their own copy of pst */
  if (my_worker== 0) {
    my_ppst=work_info->ppst;
  }
  else {
    my_ppst = &my_pst;
    memcpy(my_ppst,work_info->ppst,sizeof(struct pstruct));
    /* #else we already have the stuff in my_pst from initialization */

    my_ppst->pam2p[0] = my_ppst->pam2p[1] = NULL;

    alloc_pam(MAXSQ, MAXSQ, my_ppst);

    npam = my_pst.nsqx;

    /* allocate local copy of pam2[][] */
    for (i=0; i<npam; i++) {
      for (j=0; j<npam; j++) {
	my_pst.pam2[0][i][j] = work_info->ppst->pam2[0][i][j];
	my_pst.pam2[1][i][j] = work_info->ppst->pam2[1][i][j];
      }
    }
  }
#endif
#if defined(PCOMPLIB)	/* PCOMPLIB */
  my_ppst = &my_pst;	/* for all workers */
  alloc_pam(my_msg.pamd1,my_msg.pamd2,my_ppst);
#ifdef MPI_SRC
  MPI_Recv(&my_pst.pam2[0][0][0],my_msg.pamd1*my_msg.pamd2,MPI_INT,0,
	   STARTTYPE3, MPI_COMM_WORLD,&mpi_status);

  MPI_Recv(&my_pst.pam2[1][0][0],my_msg.pamd1*my_msg.pamd2,MPI_INT,0,
	   STARTTYPE3, MPI_COMM_WORLD,&mpi_status);
  /* no code for profiles */

  /* get pascii (only for fasty/tfasty */
  pascii = aascii;
  MPI_Recv(pascii, sizeof(aascii), MPI_BYTE, 0, STARTTYPE4, MPI_COMM_WORLD, &mpi_status);
#endif
#endif

  /* fill in info_lib_range */
  if (my_worker == FIRSTNODE) {
    /* label library size limits */
    if (my_ppst->n1_low > 0 && my_ppst->n1_high < BIGNUM) {
      sprintf(info_lib_range," (range: %d-%d)",my_ppst->n1_low,my_ppst->n1_high);}
    else if (my_ppst->n1_low > 0) {
      sprintf(info_lib_range," (range: >%d)",my_ppst->n1_low);}
    else if (my_ppst->n1_high < BIGNUM) {
      sprintf(info_lib_range," (range: <%d)",my_ppst->n1_high);}
    else {
      info_lib_range[0]='\0';
    }
    info_lib_range[sizeof(info_lib_range)-1]='\0';
#ifndef PCOMPLIB
    strncpy(work_info->info_lib_range,info_lib_range,MAX_SSTR);
    /* this does not work on some architectures */
    work_info->f_str_ap = &f_str[0];
#endif
  }

#ifdef PCOMPLIB
#ifdef MPI_SRC
  /* send back sync message */
  int_msg_b[0]=my_worker;
  MPI_Send(int_msg_b,1,MPI_INT,0,MSEQTYPE0,MPI_COMM_WORLD);
  if (my_worker == FIRSTNODE) {
    MPI_Send(info_lib_range,MAX_FN,MPI_BYTE,0,MSEQTYPE0,MPI_COMM_WORLD);
  }
#endif
#endif

  /* do the aa0[] stuff after m_msg/my_pst are initialized, for later
     inclusion in a loop */

#ifdef PCOMPLIB
#ifdef MPI_SRC
  MPI_Recv(int_msg_b,2,MPI_INT,0,
	   QSEQTYPE0, MPI_COMM_WORLD, &mpi_status);

  n0 = int_msg_b[0];
  nm0 = int_msg_b[1];
#endif
#else	/* COMP_THR */
  n0 = my_msp->n0;
  nm0 = my_msp->nm0;
  if (my_worker != FIRSTNODE) {
    /* if this is a pssm search, allocate local copy of pam2p[][]*/
    if (work_info->ppst->pam_pssm && work_info->ppst->pam2p[0]) {
      my_ppst->pam2p[0] = alloc_pam2p(my_ppst->pam2p[0],n0,npam);
      my_ppst->pam2p[1] = alloc_pam2p(my_ppst->pam2p[1],n0,npam);

      for (i=0; i<n0; i++) {
	for (j=0; j < npam; j++) {
	  my_pst.pam2p[0][i][j] = work_info->ppst->pam2p[0][i][j];
	  my_pst.pam2p[1][i][j] = work_info->ppst->pam2p[1][i][j];
	}
      }
    }
  }
#endif

  if ((aa0[0]=(unsigned char *)calloc((size_t)n0+2+SEQ_PAD,sizeof(unsigned char)))
      ==NULL) {
    fprintf(stderr," cannot allocate aa00[%d] for worker %d\n",
	    n0, my_worker);
    exit(1);
  }
  *aa0[0]='\0';
  aa0[0]++;

#ifndef PCOMPLIB
  memcpy(aa0[0],work_info->aa0,n0+1);
#else
#ifdef MPI_SRC
  /* get aa0[0] from host */
  MPI_Recv(aa0[0],n0+1,MPI_BYTE,0,
	   QSEQTYPE1,MPI_COMM_WORLD, &mpi_status);

  /* also get annotation if available */
  if (my_msp->ann_flg && my_msp->aa0a != NULL) {
    if ((my_msp->aa0a = (unsigned char *)calloc(my_msp->n0+2,sizeof(char)))==NULL) {
      fprintf(stderr, "*** error -- cannot allocate annotation array\n");
      exit(1);
    }
    MPI_Recv(my_msp->aa0a, (my_msp->n0+2)*sizeof(char), MPI_BYTE, 0,
	     QSEQTYPE1, MPI_COMM_WORLD, &mpi_status);
  }
#endif
#endif

  init_aa0(aa0, n0, nm0, &aa0s, &aa1s,
	   my_msp->qframe, my_msp->qshuffle, my_msp->max_tot,
	   my_ppst,  &f_str[0], &qf_str, my_rand_state);

/* **************************************************************** */
/* main work loop */

  while (get_wbuf(&cur_buf,max_work_buf)) {

    if (cur_buf->hdr.stop_work) break;

    /* exit thread on specific command -- this option is not used 
       for threads - get_wbuf() stops when rbuf_done() sets reader_done==1
       but it is used for PCOMPLIB
    */

    if (cur_buf->hdr.buf2_cnt <= 0) {	/* buffers can be empty */
      cur_buf->hdr.have_results = 0;
      goto res_done;
    }

    if (cur_buf->hdr.buf2_type & BUF2_DOWORK) {

      buf_do_work(aa0, n0, cur_buf, my_msp->nitt1, my_ppst, f_str);

      if (my_msp->qshuffle) {
	buf_qshuf_work(aa0s, n0, cur_buf, my_msp->nitt1,
		       my_ppst, qf_str, my_ppst->score_ix);
      }
    }

    if (cur_buf->hdr.buf2_type & BUF2_DOSHUF) {
      buf_shuf_work(aa0, n0, aa1s,  cur_buf, my_msp->nitt1,
		    my_ppst, f_str, my_ppst->score_ix, my_rand_state);
    }

    /*
    if (cur_buf->hdr.buf2_type & BUF2_DOOPT) {
      buf_do_opt(aa0, n0, cur_buf, my_ppst, f_str);
    }
    */

    if (cur_buf->hdr.buf2_type & BUF2_DOALIGN) {
      buf_do_align(aa0, n0, cur_buf, my_ppst, my_msp, f_str);
    }
    cur_buf->hdr.have_results = 1;

  res_done:
    cur_buf->hdr.have_data = 0;

    put_wbuf(cur_buf,max_work_buf);

  } /* end main while */

/* **************************************************************** */
/* all done - clean-up */

  close_work(aa0[0], n0, my_ppst, &f_str[0]);
  free(aa0[0]-1);
  if (my_msp->qframe == 2) {
    close_work(aa0[1], n0, my_ppst, &f_str[1]);
    free(aa0[1]-1);
  }

  if (my_msp->qshuffle) {
    close_work(aa0s, n0, my_ppst, &qf_str);
    free(aa0s-1);
  }

  free(aa1s-1);

#ifdef PCOMPLIB
  if (my_msp->ann_flg && my_msp->aa0a) { free(my_msp->aa0a);}
#endif

  if (my_worker) {
    free(my_pst.pam2[1][0]);
    free(my_pst.pam2[0][0]);
    free(my_pst.pam2[1]);
    free(my_pst.pam2[0]);
  }

  if (my_worker && my_pst.pam_pssm) {
    free_pam2p(my_pst.pam2p[0]);
    free_pam2p(my_pst.pam2p[1]);
  }

/* **************************************************************** */
/* and exit */

#ifdef DEBUG
  /*   fprintf(stderr,"worker [%d] done\n",my_worker); */
#endif

#ifndef PCOMPLIB
  free(my_rand_state);
  THR_EXIT(&work_info->status);
#else
  /* the PCOMPLIB version loops after a search, waiting for another max_work_buf */
  /* max_work_buf==0 signals end of queries */
  goto pcomp_loop;

 pcomp_final:
  free(my_rand_state);
#endif
}  /* end work_thread */
