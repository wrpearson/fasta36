#include <stdio.h>
#include <stdlib.h>
#include <sys/time.h>

void *my_srand(int);

main(argc, argv)
     int argc; char **argv;
{
  int i, n, s;
  struct timeval t;
  void *my_rand_state;

  if (argc < 2) n = 10;
  else n = atoi(argv[1]);

  /*
    gettimeofday(&t,NULL);
    printf(" seed: %d\n",t.tv_usec);
  */  

  my_rand_state = my_srand(0);

  for (i=0; i<n; i++) {
    s = my_nrand(n,my_rand_state);
    printf("%d\n",s);
  }

  /*
  for (i=0; i< 9999; i++) {

  }
  n = my_nrand(2147483648,my_rand_state);
  printf("number 10000: %d\n",n);
  n = my_nrand(2147483648,my_rand_state);
  printf("number 10001: %d\n",n);
  */
}
