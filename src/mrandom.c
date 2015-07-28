/* mrandom.c 28-Jan-2010 */

/*  $Id:  */
/* $Revision: 625 $  */

/* system versions of random/nrand48/nrand tend have thread contention
   issues.  This version uses a random number generator from Wikipedia
   that maintains state in a separate buffer, so that there is no contention.
*/

#include <stdio.h>
#include <stdlib.h>
#ifdef UNIX
#include <sys/time.h>
#else
#include <time.h>
#endif

/* minimal standard random number generator taken from:
   S. K. Park and K. W. Miller (1988) "Random number generators: Good
   ones are hard to find" Comm. ACM 31:1192-1201
*/
#define MIN_STD_RAND

struct m_rand_struct {
#ifndef MIN_STD_RAND
  unsigned int mw;
  unsigned int mz;
#else
  int seed;
#endif
};

#ifdef MIN_STD_RAND
#define A 16807
#define M 2147483647
#define Q 127773	/* M / A */
#define R 2836		/* M % A */
#endif

void *
my_srand(int set)	/* initialize random number generator */
{
#ifdef UNIX
  struct timeval t;
#endif
  int n;
  struct m_rand_struct *my_rand_state;

  if ((my_rand_state = (struct m_rand_struct *)calloc(1, sizeof(struct m_rand_struct)))==NULL) {
    fprintf(stderr," *** [my_srand] cannot allocate random state ***\n");
    exit(1);
  }
  
#ifdef UNIX
  gettimeofday(&t,NULL);
  n = t.tv_usec % 65535;
#else
  n = time(NULL);
#endif
  if ((n % 2)==0) n++;

#ifndef MIN_STD_RAND
  my_rand_state->mw = n;
  /* swap things around, since the next time will be close */
  n = ((n & 0xFFF) << 12) + ((n>>12) & 0xFFF);
  if ((n%2)==0) n++;
  my_rand_state->mz = n;
#else
  if (set > 0) {  my_rand_state->seed = set;}
  else {my_rand_state->seed = n;}
#endif
  return my_rand_state;
}

/* returns a random number between 0 and n-1 where n < 2^31) */
unsigned int
my_nrand(int n, struct m_rand_struct *my_rand_state)
{
  unsigned int rn;
#ifdef MIN_STD_RAND
  int lo, hi, test;

  hi = my_rand_state->seed / Q;
  lo = my_rand_state->seed % Q;
  test = A * lo - R * hi;
  if (test > 0) { my_rand_state->seed = test;}
  else {my_rand_state->seed = test + M;}
  rn = my_rand_state->seed;
#else 
  my_rand_state->mz = 36969 * (my_rand_state->mz & 65535) + (my_rand_state->mz >> 16);
  my_rand_state->mw = 18000 * (my_rand_state->mw & 65535) + (my_rand_state->mw >> 16);
  rn =  (my_rand_state->mz << 16) + my_rand_state->mw;  /* 32-bit result */
#endif

  return rn%n;
}
