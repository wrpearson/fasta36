/***************************************/
/* thread global variable declarations */
/***************************************/

/* $Id: thr_bufs2.h 625 2011-03-23 17:21:38Z wrp $ */

#ifndef MAX_WORKERS
#define MAX_WORKERS 2
#endif

#ifndef XTERNAL
struct buf_head **worker_buf;  /* pointers to full buffers */
struct buf_head **reader_buf;  /* pointers to empty buffers */

/* protected by worker_mutex/worker_cond_var */
 /* indices into full-buffers ptrs */
int worker_buf_workp;	/* modified by get_wbuf() */
int worker_buf_readp;	/* modified by put_rbuf() */
int num_worker_bufs;
int reader_done;

/* protected by reader_mutex/reader_cond var */
 /* indices into empty-buffers ptrs */
int reader_buf_workp;	/* modified by put_wbuf() */
int reader_buf_readp;	/* modified by get_rbuf(), main()-- after rbuf_wait */
int num_reader_bufs;
int reader_wait;

/* protected by start_mutex/start_cont_var */
int start_thread=1;        /* start-up predicate, 0 starts */
#else
extern struct buf_head **worker_buf;
extern struct buf_head **reader_buf;
extern int num_worker_bufs, reader_done;
extern int num_reader_bufs, reader_wait;
extern int worker_buf_workp, worker_buf_readp;
extern int reader_buf_workp, reader_buf_readp;

extern int start_thread;
#endif

