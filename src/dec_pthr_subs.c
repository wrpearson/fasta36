/* $Id: dec_pthr_subs.c 625 2011-03-23 17:21:38Z wrp $ */

/* copyright (c) 1996, 1997, 1998, 1999 by William R. Pearson and the
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


/* this file isolates the pthreads calls from the main program */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <sys/types.h>
#include <signal.h>

#include "param.h"

#include <pthread.h>
#define XTERNAL
#include "thr.h"
#undef XTERNAL
#include "pthr_subs.h"

extern void work_thread (struct thr_str *work_info);

/* start the threads working */

void init_thr(int nthreads, struct thr_str *work_info)
{
  int status, i;
  pthread_attr_t thread_attr;

  if (nthreads > MAX_WORKERS) {
    fprintf ( stderr," cannot start %d threads, max: %d\n",
	      nthreads, MAX_WORKERS);
    exit(1);
  }

  /* mutex and condition variable initialisation */

  status = pthread_mutex_init(&reader_mutex, pthread_mutexattr_default);
  check(status,"Reader_mutex init bad status\n");
   
  status = pthread_mutex_init(&worker_mutex, pthread_mutexattr_default);
  check(status,"Worker_mutex init bad status\n");

  status = pthread_cond_init(&reader_cond_var, pthread_condattr_default);
  check(status,"Reader_cond_var init bad status\n");

  status = pthread_cond_init(&worker_cond_var, pthread_condattr_default);
  check(status,"Worker_cond_var init bad status\n");

  status = pthread_mutex_init(&start_mutex, pthread_mutexattr_default);
  check(status,"Start_mutex init bad status\n");

  status = pthread_cond_init(&start_cond_var, pthread_condattr_default);
  check(status,"Start_cond_var init bad status\n");

  /* change stacksize on threads */    /***************************/

  status = pthread_attr_create( &thread_attr );
  check(status,"attribute create bad status\n");

  status = pthread_attr_setstacksize( &thread_attr, 1000000);
  check(status,"stacksize change bad status\n");

  /* start the worker threads */

  for (work_info->worker=0; work_info->worker < nthreads;
       work_info->worker++) {
    /**********************/
    status=pthread_create(&threads[work_info->worker],thread_attr,
			  (pthread_startroutine_t)&work_thread,
			  (pthread_addr_t)work_info);
    check(status,"Pthread_create failed\n");
  }
}

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

void get_rbuf(struct buf_head **cur_buf, int max_work_buf)
{
  int status;

  status = pthread_mutex_lock(&reader_mutex);  /* lock reader_buf structure */

  check(status,"Reader_mutex lock in master bad status\n");

  /* no reader bufs:  wait for signal to proceed */
  while (num_reader_bufs == 0) {
    pthread_cond_wait(&reader_cond_var,&reader_mutex);
  }

  *cur_buf = reader_buf[reader_buf_readp];  /* get the buffer address */
  reader_buf_readp = (reader_buf_readp+1)%(max_work_buf);  /* increment index */
  num_reader_bufs--;

  status = pthread_mutex_unlock(&reader_mutex);  /* unlock structure */
  check(status,"Reader_mutex unlock in master bad status\n");
}

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

  reader_done = 1;
  status = pthread_cond_broadcast(&worker_cond_var);

  status = pthread_mutex_unlock(&worker_mutex);
  check(status,"Worker_mutex unlock in master bad status\n"); 

  /* wait for all buffers available (means all do_workers are done) */
 
  for (i=0; i < nthreads; i++) {
    status = pthread_join( threads[i], &exit_value);
    check(status,"Pthread_join bad status\n");

    status = pthread_detach( &threads[i]);
    check(status,"Pthread_detach bad status\n");
  } 
}

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

void put_wbuf(struct buf_head *cur_buf, int max_work_buf)
{
  int status;

  /* put buffer back on list for reader */
  status = pthread_mutex_lock(&reader_mutex);
  check(status,"Reader_mutex lock in worker bad status\n");
    
  reader_buf[reader_buf_workp] = cur_buf;
  reader_buf_workp = (reader_buf_workp+1)%(max_work_buf);
  num_reader_bufs++;

     /* No reader_bufs available:  wake reader */
  if (num_reader_bufs == 1) {
    pthread_cond_signal(&reader_cond_var);
  }

  status = pthread_mutex_unlock(&reader_mutex);
  check(status,"Reader_mutex unlock in worker bad status\n");
}
