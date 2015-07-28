/* $Id: pthr_subs2.c 625 2011-03-23 17:21:38Z wrp $ */

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

/* modified to do more initialization of work_info here, rather than in main() */

/* this file isolates the pthreads calls from the main program */

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

#include <pthread.h>
#define XTERNAL
#include "thr_bufs2.h"
#undef XTERNAL
#include "pthr_subs.h"

extern void work_thread (struct thr_str *);

/* start the threads working */

void init_thr(int nthreads, struct thr_str *work_info,
	      const struct mngmsg *m_msp, struct pstruct *ppst,
	      unsigned char *aa0, struct mng_thr *m_bufi_p)
{
  int status, i;
  pthread_attr_t thread_attr;

  if (fa_threads == NULL && (fa_threads=(pthread_t *)calloc(nthreads, sizeof(pthread_t)))==NULL) {
      fprintf(stderr, "Cannot allocate %d pthread_t\n",nthreads);
      exit(1);
  }

  /* set up work_info[] structure, set parameters */

  for (i=0; i<nthreads; i++) {
    work_info[i].m_msp = m_msp;
    work_info[i].n0 = m_msp->n0;
    work_info[i].nm0 = m_msp->nm0;
    work_info[i].qframe = m_msp->qframe;
    work_info[i].qshuffle = m_msp->qshuffle;
    work_info[i].ppst = ppst;
    work_info[i].aa0 = aa0;
    work_info[i].max_work_buf=m_bufi_p->max_work_buf;
    work_info[i].worker=i;
    work_info[i].max_tot=m_msp->max_tot;
  }

  /* mutex and condition variable initialisation */

  status = pthread_mutex_init(&reader_mutex, NULL);
  check(status,"Reader_mutex init bad status\n");
   
  status = pthread_mutex_init(&worker_mutex, NULL);
  check(status,"Worker_mutex init bad status\n");

  status = pthread_cond_init(&reader_cond_var, NULL);
  check(status,"Reader_cond_var init bad status\n");

  status = pthread_cond_init(&worker_cond_var, NULL);
  check(status,"Worker_cond_var init bad status\n");

  status = pthread_mutex_init(&start_mutex, NULL);
  check(status,"Start_mutex init bad status\n");

  status = pthread_cond_init(&start_cond_var, NULL);
  check(status,"Start_cond_var init bad status\n");

  /* change stacksize on threads */    /***************************/

  status = pthread_attr_init( &thread_attr );
  check(status,"attribute create bad status\n");

#ifdef IRIX
  if (pthread_attr_setscope( &thread_attr, 2) != NULL) 
    status = pthread_attr_setscope( &thread_attr,PTHREAD_SCOPE_PROCESS);
  check(status,"set scope on IRIX bad status\n");
#endif

#ifdef FASTA_setscope
  status = pthread_attr_setscope( &thread_attr, PTHREAD_SCOPE_SYSTEM);
  check(status,"set scope bad status\n");
#endif

  /* start the worker threads */

  for (i=0; i < nthreads; i++) {
    /**********************/
    status=pthread_create(&fa_threads[i],&thread_attr,
			  (void *(*)(void *))&work_thread,&work_info[i]);
    check(status,"Pthread_create failed\n");
  }
}

/* start_mutex/start_cont_var provides exclusive access to 
   extern int start_thread */

void start_thr()
{
  int status;

  /* tell threads to proceed */

  status = pthread_mutex_lock(&start_mutex);
  check(status,"Start_mutex lock bad status in main\n");

  start_thread = 0;  /* lower predicate */

  status = pthread_cond_broadcast(&start_cond_var);
  status = pthread_mutex_unlock(&start_mutex);
  check(status,"Start_mutex unlock bad status in main\n");
}

/* get_rbuf() provides buffers containing sequences to the main program 
   initially, max_work_buf buffers are allocated and are available.
   As the main program runs, it calls get_rbuf() to get a reader
   buffer, fills it with sequences, and puts it on the queue with
   put_rbuf().

   At the same time, the worker programs call get_wbuf(), which gets a
   filled buffer put on the queue by put_rbuf(), takes the sequences
   from the buffer and does the comparisons, and puts the results back
   in that buffer, finally calling put_wbuf().

   locks on reader_mutex
   increments reader_buf_readp
*/

/* wait until all reader bufs are available */
void get_rbuf(struct buf_head **cur_buf, int max_work_buf)
{
  int status;

  status = pthread_mutex_lock(&reader_mutex);  /* lock reader_buf structure */

  reader_wait = 0;

  check(status,"Reader_mutex lock in master bad status\n");

  /* no reader bufs:  wait for signal to proceed */
  while (num_reader_bufs == 0) {
    pthread_cond_wait(&reader_cond_var,&reader_mutex);
  }

  *cur_buf = reader_buf[reader_buf_readp];  /* get the buffer address */
  reader_buf_readp = (reader_buf_readp+1)%(max_work_buf);  /* increment index */
  num_reader_bufs--;

  /* fprintf(stderr, " rb: %3d consumed by %lld\n",num_reader_bufs,pthread_self()); */

  status = pthread_mutex_unlock(&reader_mutex);  /* unlock structure */
  check(status,"Reader_mutex unlock in master bad status\n");
}

/* put_rbuf() takes a buffer filled with sequences to be compared and
   puts it in the queue.

   locks on worker_mutex
   increments worker_buf_readp;
 */

void put_rbuf(struct buf_head *cur_buf, int max_work_buf)
{
  int status;

  /* give the buffer to a thread, and wait for more */
  status = pthread_mutex_lock(&worker_mutex);  /* lock worker_buf_structure */
  check(status,"Worker_mutex lock in master bad status\n");

  /*  Put buffer onto available for workers list */
  worker_buf[worker_buf_readp] = cur_buf;
  worker_buf_readp = (worker_buf_readp+1)%(max_work_buf);
  num_worker_bufs++;   /* increment number of buffers available to workers */

  /*  Signal one worker to wake and start work */
  status = pthread_cond_signal(&worker_cond_var);

  status = pthread_mutex_unlock(&worker_mutex);
  check(status,"Worker_mutex unlock in master bad status\n"); 
}

/* this function is not currently used */
void put_rbuf_done(int nthreads, struct buf_head *cur_buf, int max_work_buf)
{
  int status, i;
  void *exit_value;

  /* give the buffer to a thread, and wait for more */
  status = pthread_mutex_lock(&worker_mutex);  /* lock worker_buf_structure */
  check(status,"Worker_mutex lock in master bad status\n");

  /*  Put buffer onto available for workers list */
  worker_buf[worker_buf_readp] = cur_buf;
  worker_buf_readp = (worker_buf_readp+1)%(max_work_buf);
  num_worker_bufs++;   /* increment number of buffers available to workers */

  /*  Signal one worker to wake and start work */

  reader_done = 1;	/* this causes the next get_wbuf() in the
			   thread to return 0 which causes the thread
			   to exit the main while() loop and quit */

  status = pthread_cond_broadcast(&worker_cond_var);

  status = pthread_mutex_unlock(&worker_mutex);
  check(status,"Worker_mutex unlock in master bad status\n"); 

  /* wait for all buffers available (means all do_workers are done) */
 
  for (i=0; i < nthreads; i++) {
    status = pthread_join( fa_threads[i], &exit_value);
    check(status,"Pthread_join bad status\n");
  } 
}

/* wait_rbuf() -- wait for the worker threads to finish with the
   current sequence buffers.
*/
void wait_rbuf(int used_reader_bufs) {
  int status;

  /* get a buffer to work on */
  status = pthread_mutex_lock(&reader_mutex);
  check(status,"wait_rbuf reader_mutex lock in worker bad status\n");

  reader_wait = 1;

  /*  worker_bufs still available:  wait for worker to consume them */
  while (num_reader_bufs < used_reader_bufs) {
    /* wait for a signal that no more worker_bufs are available */
    pthread_cond_wait(&reader_cond_var,&reader_mutex);
  }

  status = pthread_mutex_unlock(&reader_mutex);
  check(status,"wait_rbuf reader_mutex unlock in worker bad status\n");
}

/*
  Called once at the end of each search.  
 */

void rbuf_done(int nthreads)
{
  int status, i;
  void *exit_value;

  /* give the buffer to a thread, and wait for more */
  status = pthread_mutex_lock(&worker_mutex);  /* lock worker_buf_structure */
  check(status,"Worker_mutex lock in master bad status\n");

  /*  Signal one worker to wake and start work */

  reader_done = 1;	/* this causes the next get_wbuf() in the
			   thread to return 0 which causes the thread
			   to exit the main while() loop and quit */

  status = pthread_cond_broadcast(&worker_cond_var);

  status = pthread_mutex_unlock(&worker_mutex);
  check(status,"Worker_mutex unlock in master bad status\n"); 

  /* wait for all buffers available (means all do_workers are done) */
 
  for (i=0; i < nthreads; i++) {
    status = pthread_join( fa_threads[i], &exit_value);
    check(status,"Pthread_join bad status\n");
  } 
}

/* wait for extern int start_thread == 0 */

void wait_thr()
{
  int status;

  /* Wait on master to give start signal */
  status = pthread_mutex_lock(&start_mutex);
  check(status,"Start_mutex lock bad status in worker\n");

  while (start_thread) {
         status = pthread_cond_wait(&start_cond_var, &start_mutex);
         check(status,"Start_cond_wait bad status in worker\n");
  }

  status = pthread_mutex_unlock(&start_mutex);
  check(status,"Start_mutex unlock bad status in worker\n");
}

/* get_wbuf() -- used in worker threads
   get a buffer full of sequences to be compared from the main program

   locks on worker_mutex
   increments worker_buf_workp
 */

int get_wbuf(struct buf_head **cur_buf, int max_work_buf)
{
  int status;

  /* get a buffer to work on */
  status = pthread_mutex_lock(&worker_mutex);
  check(status,"First worker_mutex lock in worker bad status\n");

  /*  No worker_bufs available:  wait for reader to produce some */
  while (num_worker_bufs == 0) {
    /*  Exit if reader has finished */
    if (reader_done) {
      pthread_mutex_unlock(&worker_mutex);
      return 0;
    }
    pthread_cond_wait(&worker_cond_var,&worker_mutex);
  } /* end while */

  /*  Get the buffer from list */
  *cur_buf = worker_buf[worker_buf_workp];
  worker_buf_workp = (worker_buf_workp+1)%(max_work_buf);
  num_worker_bufs--;

  status = pthread_mutex_unlock(&worker_mutex);
  check(status,"First worker_mutex unlock in worker bad status\n");
  return 1;
}

/* put_wbuf() -- called in worker threads
   return a filled buffer of scores to the main program

   locks on reader_mutex
   increments reader_buf_workp
 */
void put_wbuf(struct buf_head *cur_buf, int max_work_buf)
{
  int status;

  /* put buffer back on list for reader */
  status = pthread_mutex_lock(&reader_mutex);
  check(status,"Reader_mutex lock in worker bad status\n");
    
  reader_buf[reader_buf_workp] = cur_buf;
  reader_buf_workp = (reader_buf_workp+1)%(max_work_buf);
  num_reader_bufs++;

  /*   fprintf(stderr, " rb: %3d produced\n",num_reader_bufs); */

  /* we used to do this only when num_reader_bufs==1 */
  if (num_reader_bufs == 1 || reader_wait == 1) {
    pthread_cond_signal(&reader_cond_var);
  }

  status = pthread_mutex_unlock(&reader_mutex);
  check(status,"Reader_mutex unlock in worker bad status\n");
}
