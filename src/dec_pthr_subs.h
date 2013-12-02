
/* $Id: dec_pthr_subs.h 625 2011-03-23 17:21:38Z wrp $ */
/* $Revision: 625 $  */

#include <pthread.h>

#define check(status,string) \
     if (status == -1) perror(string)   /* error macro for thread calls */

#ifndef XTERNAL
pthread_t threads[MAX_WORKERS];

/* mutex stuff */

pthread_mutex_t reader_mutex;      /* empty buffer pointer structure lock */
pthread_mutex_t worker_mutex;      /* full buffer pointer structure lock */

/* condition variable stuff */

pthread_cond_t reader_cond_var;    /* condition variable for reader */
pthread_cond_t worker_cond_var;    /* condition variable for workers */

pthread_mutex_t start_mutex;       /* start-up synchronisation lock */
pthread_cond_t start_cond_var;     /* start-up synchronisation condition variable */

extern pthread_t threads[];

/* mutex stuff */

extern pthread_mutex_t reader_mutex;
extern pthread_mutex_t worker_mutex;

/* condition variable stuff */

extern pthread_cond_t reader_cond_var;
extern pthread_cond_t worker_cond_var;

extern pthread_mutex_t start_mutex;
extern pthread_cond_t start_cond_var;
extern int start_thread;

#endif
