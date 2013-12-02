
/***************************************/
/* thread global variable declarations */
/***************************************/

/* $Id: pcomp_bufs.h 625 2011-03-23 17:21:38Z wrp $ */
/* $Revision: 625 $  */

/* this file serves the same purpose as thr_bufs2.h in the threaded version

   Unlike thr_bufs2.h, which has two separate lists of buffers,
   reader_buf[] and worker_buf[], so that the threads can be asynchronous.

   pcomp workers are asynchonous as well, but, initially, there is an
   array that has an entry for each worker (worker_buf[]), and a queue
   that captures which buffers are ready (work_q[]).
 */

#ifndef MAX_WORKERS
#define MAX_WORKERS 2
#endif

#ifndef XTERNAL
struct buf_head **worker_buf;  /* pointers buffers of sequences/results */
int *work_q;	/* next worker available */
int max_worker_q;
int num_reader_bufs;
#else
extern struct buf_head **worker_buf;
extern int *work_q;
extern int max_worker_q;
extern int num_reader_bufs;
#endif
