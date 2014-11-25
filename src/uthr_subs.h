/***************************************/
/* thread global variable declarations */
/***************************************/

/* $Id: uthr_subs.h 625 2011-03-23 17:21:38Z wrp $ */

#ifndef MAX_WORKERS
#define MAX_WORKERS 2
#endif
#define NUM_WORK_BUF 2*MAX_WORKERS

#include <synch.h>
#include <thread.h>

#define check(status,string) \
     if (status == -1) perror(string)   /* error macro for thread calls */

#ifndef XTERNAL

thread_t threads[MAX_WORKERS];

/* mutex stuff */

mutex_t reader_mutex;      /* empty buffer pointer structure lock */
mutex_t worker_mutex;      /* full buffer pointer structure lock */

/* condition variable stuff */

cond_t reader_cond_var;    /* condition variable for reader */
cond_t worker_cond_var;    /* condition variable for workers */

mutex_t start_mutex;       /* start-up synchronisation lock */
cond_t start_cond_var;     /* start-up synchronisation condition variable */

#else

extern thread_t threads[];

/* mutex stuff */

extern mutex_t reader_mutex;
extern mutex_t worker_mutex;

/* condition variable stuff */

extern cond_t reader_cond_var;
extern cond_t worker_cond_var;

extern mutex_t start_mutex;
extern cond_t start_cond_var;

#endif
