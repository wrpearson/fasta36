/* $Id: pcomp_subs2.c $ */

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

/* modified to do more initialization of work_info here, rather than
   in main() */

/* this file provides the same functions for PCOMPLIB as pthr_subs2.c does for COMP_THR */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <sys/types.h>
#ifdef UNIX
#include <unistd.h>
#endif
#include <signal.h>

#include "defs.h"
#include "structs.h"		/* mngmsg, libstruct */
#include "param.h"		/* pstruct, thr_str, buf_head, rstruct */
#include "thr_buf_structs.h"

#ifdef MPI_SRC
#include "mpi.h"
#endif

#include "msg.h"
#include "pcomp_bufs.h"
#define XTERNAL
#include "uascii.h"
#undef XTERNAL
#include "pthr_subs.h"

#ifdef DEBUG
unsigned long adler32(unsigned long, const unsigned char *, unsigned int);
#endif

static int next_worker_idx, num_workers_idle;
extern int g_worker;

/* used for debugging */
/*
int check_seq_range(unsigned char *aa1b, int n1, int nsq, char *);
*/

/* start the workers, nworkers == number of workers, not nodes */
void
init_thr(int nworkers, char *info_lib_range_p, const struct mngmsg *m_msp, struct pstruct *ppst,
	 unsigned char *aa0, struct mng_thr *m_bufi_p)
{
#ifdef MPI_SRC
  MPI_Status mpi_status;
  int int_msg_b[4];	/* general purpose buffer for integers */
#endif
  int node, snode;

  /* start the worker processes */

  if (work_q == NULL) {
    if ((work_q=(int *)calloc(nworkers, sizeof(int)))==NULL) {
      fprintf(stderr, " cannot allocate work_q[%d] structure\n",
	      nworkers);
      exit(1);
    }
    else {max_worker_q = nworkers;}
  }
  num_workers_idle = 0;

  /* setup thread buffer info */
  if (aa0 == NULL) {
    int_msg_b[0] = int_msg_b[1] = int_msg_b[2] = 0;
  }
  else {
    int_msg_b[0] = nworkers;
    int_msg_b[1] = m_bufi_p->max_buf2_res;
    int_msg_b[2] = m_bufi_p->max_chain_seqs;
    int_msg_b[3] = m_bufi_p->seq_buf_size;
  }

  /* send thread info */
  for (node=FIRSTNODE; node < nworkers+FIRSTNODE; node++) {

    MPI_Send(int_msg_b, 4, MPI_INT, node, STARTTYPE0, MPI_COMM_WORLD);

    if (aa0 == NULL) { continue;}

    /* send mngmsg */
    MPI_Send((void *)m_msp, sizeof(struct mngmsg), MPI_BYTE, node,
	     STARTTYPE1, MPI_COMM_WORLD);

    MPI_Send(ppst, sizeof(struct pstruct), MPI_BYTE, node,
	     STARTTYPE2, MPI_COMM_WORLD);

    /* send the rest of the pieces of pam[2] */
    MPI_Send(&ppst->pam2[0][0][0],m_msp->pamd1*m_msp->pamd2,MPI_INT,node,STARTTYPE3,
	     MPI_COMM_WORLD);
    MPI_Send(&ppst->pam2[1][0][0],m_msp->pamd1*m_msp->pamd2,MPI_INT,node,STARTTYPE3,
	     MPI_COMM_WORLD);

    /* send pascii (only for fasty/tfasty */
    MPI_Send(pascii, sizeof(aascii), MPI_BYTE, node, STARTTYPE4, MPI_COMM_WORLD);
  }

  if (aa0 == NULL) {
    /* all done */
    free(work_q);
    return;
  }

  /* wait for returned status results */
  while (num_workers_idle < max_worker_q) {
    MPI_Recv(&node, 1, MPI_INT, MPI_ANY_SOURCE,MSEQTYPE0,
	     MPI_COMM_WORLD, &mpi_status);
    snode= mpi_status.MPI_SOURCE;
    if (snode == FIRSTNODE) {
      MPI_Recv(info_lib_range_p, MAX_FN, MPI_BYTE, snode,MSEQTYPE0,
	       MPI_COMM_WORLD, &mpi_status);
    }

    if (snode != node) {
      fprintf(stderr, " initial node mismatch [%d!=%d]\n",node, snode);
    }
    worker_buf[snode-FIRSTNODE]->hdr.have_data = 0;
    worker_buf[snode-FIRSTNODE]->hdr.have_results = 0;
    worker_buf[snode-FIRSTNODE]->hdr.worker_idx = snode;
    work_q[num_workers_idle++] = snode;
  }
  next_worker_idx = 0;

  /* send query sequence info to workers */
  for (node=FIRSTNODE; node < nworkers+FIRSTNODE; node++) {
    /* send thread buffer info */
    int_msg_b[0] = m_msp->n0;
    int_msg_b[1] = m_msp->nm0;
    MPI_Send(int_msg_b, 2, MPI_INT, node, QSEQTYPE0, MPI_COMM_WORLD);
    MPI_Send(aa0, m_msp->n0+1, MPI_BYTE, node, QSEQTYPE1, MPI_COMM_WORLD);
    if (m_msp->ann_flg && m_msp->aa0a) {
      MPI_Send(m_msp->aa0a, m_msp->n0+2, MPI_BYTE, node, QSEQTYPE1, MPI_COMM_WORLD);
    }
  }
}

/* get_rbuf() provides buffers containing sequences to the main
   program. max_work_q buffers are available, with each
   buffer tied to a worker.

   As the main program runs, it calls get_rbuf() to get a worker
   buffer (reader buffers are not used with PCOMPLIB), fills it with
   sequences, and sends it to a worker with put_rbuf().

   At the same time, the worker programs call get_wbuf(), to get a
   filled worker buffer sent by put_rbuf(), takes the sequences from
   the buffer and does the comparisons, and sends the results back to
   the manager by calling put_wbuf().
*/

/* wait for results from any worker */
struct buf_head *
next_work_result(int *snode) {
  int this_node, buf2_cnt;
  int int_msg_b[4];	/* general purpose int buffer */
  int i;
  struct buf2_hdr_s buf2_head;
  struct buf_head *this_buf_p, tmp_buf_head;
  struct seq_record *seq_b_save;
  struct mseq_record *mseq_b_save;
  unsigned char *aa1b_start_save;
  struct a_res_str *new_ares_p, *prev_ares_p;
#ifdef MPI_SRC
  MPI_Status mpi_status;
#endif

  /* wait for a returned result */
  MPI_Recv(&tmp_buf_head, sizeof(struct buf_head), MPI_BYTE, MPI_ANY_SOURCE,RES_TYPE0,
	   MPI_COMM_WORLD, &mpi_status);
  this_node = mpi_status.MPI_SOURCE;
  buf2_cnt = tmp_buf_head.hdr.buf2_cnt;

#ifdef DEBUG
  /*
  fprintf(stderr," %d: %d results\n", this_node, buf2_cnt);
  */
#endif

  this_buf_p = worker_buf[this_node-FIRSTNODE];
  /* move things selectively to avoid over-writing pointers to res, a_res arrays */
  
  aa1b_start_save = this_buf_p->hdr.aa1b_start;
  seq_b_save = this_buf_p->hdr.seq_b;
  mseq_b_save = this_buf_p->hdr.mseq_b;

  memcpy(&this_buf_p->hdr,&tmp_buf_head.hdr,sizeof(struct buf2_hdr_s));

  this_buf_p->hdr.aa1b_start = aa1b_start_save;
  this_buf_p->hdr.seq_b = seq_b_save;
  this_buf_p->hdr.mseq_b =mseq_b_save;

  memcpy(&this_buf_p->s_cnt_info,&tmp_buf_head.s_cnt_info,sizeof(struct score_count_s));

  if (this_buf_p->hdr.have_results) {
    if (this_buf_p->hdr.buf2_type & (BUF2_DOWORK + BUF2_DOSHUF + BUF2_DOOPT)) {
      MPI_Recv(this_buf_p->buf2_res, sizeof(struct buf2_res_s)*buf2_cnt,
	       MPI_BYTE, this_node, RES_TYPE1, MPI_COMM_WORLD, &mpi_status);
      /*
      for (i=0; i < buf2_cnt; i++) {
	if (this_buf_p->buf2_res[i].rst.score[2] > 200) {
	  fprintf(stderr, "HS[%d:%d,%d]: %d (%d:%d)\n",i,this_node, tmp_buf_head.hdr.worker_idx, this_buf_p->buf2_res[i].rst.score[2],
		  this_buf_p->buf2_data[i].seq->index,this_buf_p->buf2_data[i].seq->n1);
	}
      }
      */
    }

    if (this_buf_p->hdr.buf2_type & BUF2_DOALIGN) {
      /* (1) get a message that has "have_ares"
	 (2) allocate space for each a_res and receive it individually
	 (3) reset the ->next pointers for the a_res chain
      */

      for (i = 0; i < buf2_cnt; i++) {
	MPI_Recv(int_msg_b, 1, MPI_INT, this_node, ALN_TYPE0, MPI_COMM_WORLD, &mpi_status);
	this_buf_p->buf2_ares[i].have_ares = int_msg_b[0];
	this_buf_p->buf2_ares[i].a_res = NULL;	/* pre-initialize */

	if (this_buf_p->buf2_ares[i].have_ares) {
	  /* allocate space to receive it */
	  if ((new_ares_p = (struct a_res_str *)calloc(1,sizeof(struct a_res_str)))==NULL) {
	    fprintf(stderr, "cannot allocate a_res from %d\n",this_node);
	    exit(1);
	  }
	  /* save the head of the ares_chain */
	  this_buf_p->buf2_ares[i].a_res = new_ares_p;

	  /* get the first a_res */
	  MPI_Recv(new_ares_p, sizeof(struct a_res_str), MPI_BYTE, this_node,
		   ALN_TYPE1, MPI_COMM_WORLD, &mpi_status);
	  /* get the associated res[nres] */
	  if ((new_ares_p->res = (int *)calloc(new_ares_p->nres,sizeof(int)))==NULL) {
	    fprintf(stderr, "cannot allocate res for a_res from %d\n",this_node);
	    exit(1);
	  }
	  MPI_Recv(new_ares_p->res, new_ares_p->nres, MPI_INT, this_node,
		   ALN_TYPE2, MPI_COMM_WORLD, &mpi_status);

	  /* now get alignment encodings if available */
	  if (new_ares_p->aln_code) {
	    if ((new_ares_p->aln_code = (char *)calloc(new_ares_p->aln_code_n+1,sizeof(char)))==NULL) {
	      fprintf(stderr, "cannot allocate aln_code for a_res from %d\n",this_node);
	      exit(1);
	    }
	    MPI_Recv(new_ares_p->aln_code, new_ares_p->aln_code_n+1, MPI_BYTE, this_node,
		   ALN_TYPE3, MPI_COMM_WORLD, &mpi_status);
	  }
	  if (new_ares_p->ann_code) {
	    if ((new_ares_p->ann_code = (char *)calloc(new_ares_p->ann_code_n+1,sizeof(char)))==NULL) {
	      fprintf(stderr, "cannot allocate ann_code for a_res from %d\n",this_node);
	      exit(1);
	    }
	    MPI_Recv(new_ares_p->ann_code, new_ares_p->ann_code_n+1, MPI_BYTE, this_node,
		   ALN_TYPE3, MPI_COMM_WORLD, &mpi_status);

	  }

	  while (new_ares_p->next) {	/* while the chain continues */
	    prev_ares_p = new_ares_p;	/* save pointer to previous a_res to fix prev_ares->next */
	    if ((new_ares_p = (struct a_res_str *)calloc(1,sizeof(struct a_res_str)))==NULL) {
	      fprintf(stderr, "cannot allocate a_res from %d\n",this_node);
	      exit(1);
	    }
	    prev_ares_p->next = new_ares_p;
	    MPI_Recv(new_ares_p, sizeof(struct a_res_str), MPI_BYTE, this_node,
		     ALN_TYPE1, MPI_COMM_WORLD, &mpi_status);
	    if ((new_ares_p->res = (int *)calloc(new_ares_p->nres,sizeof(int)))==NULL) {
	      fprintf(stderr, "cannot allocate res for a_res from %d\n",this_node);
	      exit(1);
	    }
	    MPI_Recv(new_ares_p->res, new_ares_p->nres, MPI_INT, this_node,
		     ALN_TYPE2, MPI_COMM_WORLD, &mpi_status);
	    /* now get alignment encodings if available */
	    if (new_ares_p->aln_code) {
	      if ((new_ares_p->aln_code = (char *)calloc(new_ares_p->aln_code_n+1,sizeof(char)))==NULL) {
		fprintf(stderr, "cannot allocate aln_code for a_res from %d\n",this_node);
		exit(1);
	      }
	      MPI_Recv(new_ares_p->aln_code, new_ares_p->aln_code_n+1, MPI_BYTE, this_node,
		       ALN_TYPE3, MPI_COMM_WORLD, &mpi_status);
	    }
	    if (new_ares_p->ann_code) {
	      if ((new_ares_p->ann_code = (char *)calloc(new_ares_p->ann_code_n+1,sizeof(char)))==NULL) {
		fprintf(stderr, "cannot allocate ann_code for a_res from %d\n",this_node);
		exit(1);
	      }
	      MPI_Recv(new_ares_p->ann_code, new_ares_p->ann_code_n+1, MPI_BYTE, this_node,
		       ALN_TYPE3, MPI_COMM_WORLD, &mpi_status);

	    }
	  }	/* finished with the ares_chain */
	} /* done with have_ares */
	else {
#ifdef DEBUG
	  fprintf(stderr, " getting alignment with no have_ares[%d]: %d/%d",
		  this_buf_p->hdr.worker_idx,i,this_buf_p->buf2_ares[i].best_idx);
#endif
	}
      }	/* done with buf2_ares[buf2_cnt] */
    } /* done with BUF_DOALIGN */
  } /* done with have_results */
  *snode = this_node;
  return this_buf_p;
}

/* wait until a worker/buffer is available */
void get_rbuf(struct buf_head **cur_buf, int max_work_buf)
{
  int node, snode;
  int i_msg_b[2], nresults;
  struct buf_head *this_buf_p;

#ifdef MPI_SRC
  MPI_Status mpi_status;
#endif

  if (num_workers_idle == 0) {
    this_buf_p = next_work_result(&snode);

    work_q[next_worker_idx] = snode;
    num_workers_idle++;
  }
  else {
    this_buf_p = worker_buf[work_q[next_worker_idx]-FIRSTNODE];
  }

  *cur_buf = this_buf_p;

  /* update worker queue */
  next_worker_idx = (next_worker_idx+1)%(max_work_buf);
}

/* put_rbuf() takes a buffer filled with sequences to be compared 
   sends it to a worker */

void put_rbuf(struct buf_head *cur_buf, int max_work_buf)
{
#ifdef MPI_SRC
  MPI_Status mpi_status;
#endif
  struct seq_record *cur_seq_p, *tmp_seq_p;
  int i, j, snode, buf2_cnt, seqr_cnt;
  int cur_aa1b_size, max_aa1b_size;

  /* do not send msg if no data */
  if (!cur_buf->hdr.have_data || !(cur_buf->hdr.buf2_cnt > 0)) {return;}

  /* here, since we have a buffer, we have a worker, just send the info */
  snode = cur_buf->hdr.worker_idx;
  buf2_cnt = cur_buf->hdr.buf2_cnt;
  seqr_cnt = cur_buf->hdr.seqr_cnt;
  max_aa1b_size = cur_buf->hdr.aa1b_size;

#ifdef DEBUG
  /*   fprintf(stderr," sending %d/%d seqs to %d\n", buf2_cnt, seqr_cnt, snode); */
#endif
  /* send header */
  MPI_Send(&cur_buf->hdr, sizeof(struct buf2_hdr_s), MPI_BYTE, snode,
	   MSEQTYPE0, MPI_COMM_WORLD);

  /* send data */
  MPI_Send(cur_buf->buf2_data, sizeof(struct buf2_data_s)*buf2_cnt,
	   MPI_BYTE, snode, MSEQTYPE1, MPI_COMM_WORLD);

  /* before sending sequence records, we need to check to see if we
     need to transfer to a continuous location (or send lots of short
     records) */

#ifdef DEBUG
  cur_aa1b_size = 0;
  for (i=0; i < buf2_cnt; i++) {
    cur_seq_p = cur_buf->buf2_data[i].seq;
    if (!cur_buf->buf2_data[i].seq_dup) {
      cur_aa1b_size += cur_seq_p->n1+1;
    }
    if (check_seq_range(cur_seq_p->aa1b, cur_seq_p->n1, 50, "put_rbuf()")) {
      fprintf(stderr, "[put_rbuf] range error at: %d\n", i);
    }
  }

  if (cur_aa1b_size != cur_buf->hdr.aa1b_used) {
    fprintf(stderr,"[put_rbuf:%d] aa1b_used size mismatch: %d != %d\n",
	    snode, cur_aa1b_size, cur_buf->hdr.aa1b_used);
  }
#endif

  if (cur_buf->hdr.seq_record_continuous) {
    /* send sequence records associated with data in one message */
    MPI_Send(cur_buf->hdr.seq_b, sizeof(struct seq_record)*seqr_cnt,
	     MPI_BYTE, snode, MSEQTYPE2, MPI_COMM_WORLD);
    MPI_Send(cur_buf->hdr.aa1b_start, cur_buf->hdr.aa1b_used+1,
	     MPI_BYTE, snode, MSEQTYPE3, MPI_COMM_WORLD);
  }
  else {
    /* send individual sequence records */
    cur_aa1b_size = 0;
    for (i=0; i < buf2_cnt; i++) {
      cur_seq_p = cur_buf->buf2_data[i].seq;
      if (!cur_buf->buf2_data[i].seq_dup) {	/* don't send sequence if its a duplicate */
	MPI_Send(cur_seq_p, sizeof(struct seq_record),
		 MPI_BYTE, snode, MSEQTYPE4, MPI_COMM_WORLD);
	MPI_Send(cur_seq_p->aa1b, cur_seq_p->n1+1,
		 MPI_BYTE, snode, MSEQTYPE5, MPI_COMM_WORLD);
      }
    }
  }

  /* reduce the number of idle workers */
  num_workers_idle--;
}

/* wait_rbuf() -- wait for the worker threads to finish with the
   current sequence buffers.
*/
void wait_rbuf(int used_reader_bufs) {
  int snode;

  while (num_workers_idle < max_worker_q) {
    next_work_result(&snode);
    num_workers_idle++;
  }

  /* all workers are idle, re-initialize work_q */
  for (snode = 0; snode < max_worker_q; snode++) {
    work_q[snode] = snode + FIRSTNODE;
  }
}

void rbuf_done(int nthreads)
{
#ifdef MPI_SRC
  MPI_Status mpi_status;
#endif
  int status, i;

  /* use a dummy buf_head to send buf2_cnt=0, stop_work=1 */
  struct buf2_hdr_s tmp_buf2_hdr;

  tmp_buf2_hdr.stop_work = 1;
  tmp_buf2_hdr.buf2_cnt = 0;

  /* send a message to all the workers waiting for get_wbuf()
     to quit
   */
 
  for (i=FIRSTNODE; i < nthreads+FIRSTNODE; i++) {
    MPI_Send(&tmp_buf2_hdr, sizeof(struct buf2_hdr_s), MPI_BYTE, i,
	     MSEQTYPE0, MPI_COMM_WORLD);
  }
}

/* get_wbuf() -- called in workers
   get a buffer full of sequences to be compared from the main program

   this function should follow put_rbuf() message for message

   In the PCOMPLIB version, there is no queue of buffers to be read,
   but we must have space to put the messages in as we receive them,
   and we must fix the pointers in the seq_records
 */

int get_wbuf(struct buf_head **cur_buf, int max_work_buf)
{
#ifdef MPI_SRC
  MPI_Status mpi_status;
#endif

  /* we need to preserve some sequence pointer information so it is not
     over-written by the messages */

  struct seq_record *seq_base, *cur_seq_p, *prev_seq_p, *old_host_seq_p, *host_seq_p;
  struct buf2_data_s *cur_buf2_dp;
  unsigned char *aa1b_start_save, *old_aa1b_start, *cur_aa1b;
  unsigned char *host_aa1b, *old_host_aa1b;
  int buf2_cnt, i, j, cur_n1, seqr_cnt;
  int max_aa1b_size, aa1b_size_save;
  int cur_aa1b_size;
  int snode;

  snode = (*cur_buf)->hdr.worker_idx;
  seq_base = (*cur_buf)->hdr.seq_b;
  aa1b_start_save = (*cur_buf)->hdr.aa1b_start;
  max_aa1b_size = aa1b_size_save = (*cur_buf)->hdr.aa1b_size;

  /* put invalid bytes in aa1b to check for transmission errors */
  memset(aa1b_start_save, 127, aa1b_size_save);

  MPI_Recv(&(*cur_buf)->hdr, sizeof(struct buf2_hdr_s), MPI_BYTE, 0, 
	   MSEQTYPE0, MPI_COMM_WORLD, &mpi_status);

  buf2_cnt = (*cur_buf)->hdr.buf2_cnt;
  seqr_cnt = (*cur_buf)->hdr.seqr_cnt;

  if (buf2_cnt <= 0 || (*cur_buf)->hdr.stop_work) { return 0; }

  /* get the buf2_data array, which has seq_dup and ->seq records */
  MPI_Recv((*cur_buf)->buf2_data, sizeof(struct buf2_data_s)*buf2_cnt,
	   MPI_BYTE, 0, MSEQTYPE1, MPI_COMM_WORLD, &mpi_status);
  
#ifdef DEBUG
  /*   fprintf(stderr,"[%d/get_wbuf] receiving %d/%d sequences\n",snode, buf2_cnt, seqr_cnt); */
#endif

  /* get seq_records (but not mseq_records, don't need them) */
  if ((*cur_buf)->hdr.seq_record_continuous) {
    MPI_Recv(seq_base, sizeof(struct seq_record)*seqr_cnt,
	     MPI_BYTE, 0, MSEQTYPE2, MPI_COMM_WORLD, &mpi_status);

    /* now get the sequence data */
    MPI_Recv(aa1b_start_save, (*cur_buf)->hdr.aa1b_used+1,
	     MPI_BYTE, 0, MSEQTYPE3, MPI_COMM_WORLD, &mpi_status);
  
    /* map the seq records back into buf2_data */
    /* must check for duplicate sequence records, initialize buf2_data[i]->seq 
       AND seq.aa1b in the same pass */

    cur_buf2_dp = (*cur_buf)->buf2_data;
    cur_seq_p = prev_seq_p = seq_base;

    cur_aa1b = aa1b_start_save;
    cur_aa1b_size = 0;

    for (i=0; i < buf2_cnt; i++, cur_buf2_dp++) {
      if (!cur_buf2_dp->seq_dup) {	/* not a duplicate */
	cur_seq_p->aa1b = cur_aa1b;
	cur_aa1b += cur_seq_p->n1 + 1;
	cur_aa1b_size += cur_seq_p->n1 + 1;
	cur_buf2_dp->seq = cur_seq_p++;
      }
      else {		/* duplicate */
	cur_buf2_dp->seq = prev_seq_p;	/* point to the previous value */
	prev_seq_p = cur_seq_p;
      }
    }

    if (cur_aa1b_size != (*cur_buf)->hdr.aa1b_used) {
      fprintf(stderr, "[%d] incorrect cur_aa1b_size: %d != %d [%d]\n",
	      snode, cur_aa1b_size, (*cur_buf)->hdr.aa1b_used);
    }
  }
  else { /* not continuous, get seq_records one at a time */
    cur_seq_p = seq_base;
    cur_aa1b = aa1b_start_save;
    cur_buf2_dp = (*cur_buf)->buf2_data;
    cur_aa1b_size = 0;
    for (i=0; i < buf2_cnt; i++) {
      /* get a seq record */
      if (!(*cur_buf)->buf2_data[i].seq_dup) {	/* not a duplicate, so get it */
	MPI_Recv(cur_seq_p, sizeof(struct seq_record),
		 MPI_BYTE, 0, MSEQTYPE4, MPI_COMM_WORLD, &mpi_status);
	/* get the sequence itself */
	prev_seq_p = cur_seq_p;
	cur_n1 = cur_seq_p->n1;
	cur_aa1b_size += cur_n1+1;
	if (cur_aa1b_size >= max_aa1b_size) {
	  fprintf(stderr,"[get_wbuf:%d] -- receive buffer too small %d > %d\n",
		  (*cur_buf)->hdr.worker_idx, cur_aa1b_size, max_aa1b_size);
	  exit(1);
	}

	MPI_Recv(cur_aa1b, cur_n1+1, MPI_BYTE, 0, MSEQTYPE5, MPI_COMM_WORLD, &mpi_status);
	cur_seq_p->aa1b = cur_aa1b;
#ifdef DEBUG
	if (cur_seq_p->adler32_crc != adler32(1L,cur_aa1b,cur_n1)) {
	  fprintf(stderr," [get_wbuf:%d] -- adler32 mismatch; n1: %d\n",
		  (*cur_buf)->hdr.worker_idx, cur_n1);
	}
#endif

	cur_buf2_dp->seq = cur_seq_p++;
	cur_aa1b += cur_n1+1;
      }
      else {	/* its a duplicate, so point to the original version */
	cur_buf2_dp->seq = prev_seq_p;
      }
      cur_buf2_dp++;
    }
  }

  /* restore the seq_b, aa1b_start that were over-written */
  (*cur_buf)->hdr.seq_b = seq_base;
  (*cur_buf)->hdr.aa1b_start = aa1b_start_save;
  (*cur_buf)->hdr.aa1b_size = aa1b_size_save;

  /*
  for (i=0; i < buf2_cnt; i++) {
    cur_seq_p = (*cur_buf)->buf2_data[i].seq;
    if (check_seq_range(cur_seq_p->aa1b, cur_seq_p->n1, 50, "get_wbuf()")) {
      fprintf(stderr, "[%d] (get_wbuf) range error at: %d/%d (seqr_cnt: %d)\n", 
	      (*cur_buf)->hdr.worker_idx, i, buf2_cnt, seqr_cnt);
    }
  }
  */

  return 1;
}

/* put_wbuf() -- called in workers

   In the PCOMPLIB version, there is no queue of buffers to be read,
   so just send the buffer to the manager
 */
void put_wbuf(struct buf_head *cur_buf, int max_work_buf)
{
  int int_msg_b[4], i;
  struct buf2_ares_s *buf2_ares_p;
  struct a_res_str *cur_ares_p, *next_ares_p;
#ifdef MPI_SRC
  MPI_Status mpi_status;
#endif

  MPI_Send(&cur_buf->hdr, sizeof(struct buf_head), MPI_BYTE, 0,
	   RES_TYPE0, MPI_COMM_WORLD);

  if (!cur_buf->hdr.have_results) { return;}

  /* have buf2_res type results */
  if (cur_buf->hdr.buf2_type & (BUF2_DOWORK + BUF2_DOSHUF+BUF2_DOOPT)) {
    MPI_Send(cur_buf->buf2_res, sizeof(struct buf2_res_s)*cur_buf->hdr.buf2_cnt, MPI_BYTE, 0,
	     RES_TYPE1, MPI_COMM_WORLD);
  }

  /* have buf2_ares type results */
  if (cur_buf->hdr.buf2_type & BUF2_DOALIGN) {
    /* buf2_ares does not have much useful information, except have_ares and a chain of *a_res pointers.
       so we need to:
       (1) send have_ares
       (2) send each part of the a_res chain individually
    */
       
    buf2_ares_p = cur_buf->buf2_ares;
    for (i=0; i < cur_buf->hdr.buf2_cnt; i++) {
      int_msg_b[0] = buf2_ares_p->have_ares;
      MPI_Send(int_msg_b, 1, MPI_INT, 0, ALN_TYPE0, MPI_COMM_WORLD);
      if (buf2_ares_p->have_ares) {
	/* (a) send the first one */
	for (cur_ares_p = buf2_ares_p->a_res; cur_ares_p; cur_ares_p = cur_ares_p->next) {
	  MPI_Send(cur_ares_p, sizeof(struct a_res_str), MPI_BYTE, 0, ALN_TYPE1, MPI_COMM_WORLD);
	  MPI_Send(cur_ares_p->res, cur_ares_p->nres ,MPI_INT, 0, ALN_TYPE2, MPI_COMM_WORLD);
	  if (cur_ares_p->aln_code) {
	    MPI_Send(cur_ares_p->aln_code, cur_ares_p->aln_code_n+1 ,MPI_BYTE, 0, ALN_TYPE3, MPI_COMM_WORLD);
	  }
	  if (cur_ares_p->ann_code) {
	    MPI_Send(cur_ares_p->ann_code, cur_ares_p->ann_code_n+1 ,MPI_BYTE, 0, ALN_TYPE3, MPI_COMM_WORLD);
	  }
	} /* done with a_res chain */

	/* free the chain */
	cur_ares_p = buf2_ares_p->a_res;
	while (cur_ares_p) {
	  if (cur_ares_p->aln_code) free(cur_ares_p->aln_code);
	  if (cur_ares_p->ann_code) free(cur_ares_p->ann_code);
	  if ((buf2_ares_p->have_ares & 0x1) && cur_ares_p->res) free(cur_ares_p->res);
	  next_ares_p = cur_ares_p->next;
	  free(cur_ares_p);
	  cur_ares_p = next_ares_p;
	}
	buf2_ares_p->a_res = NULL;
      } /* done with have_ares */
      buf2_ares_p->have_ares = 0;	/* must be zero-ed out for later use */
      buf2_ares_p++;
    } /* done with buf2_ares[buf2_cnt] */
  } /* done with BUF2_DOALIGN */
} /* done with put_wbuf() */
