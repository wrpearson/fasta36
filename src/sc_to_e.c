/* $Id: sc_to_e.c 625 2011-03-23 17:21:38Z wrp $ */

/* copyright (c) 1996, 1997, 1998, 1999, 2014 by William R. Pearson and the
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

/* sc_to_e  uses statistical parameters from search and
	    score, length, and database size to calculate E()
*/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

double mean_var, mu, rho;

main(argc, argv)
     int argc; char **argv;
{
  char line[128];
  int score, length, db_size;
  double z_val, s_to_zv(), zv_to_E();

  if (argc == 4) {
    sscanf(argv[1],"%lf",&rho);
    sscanf(argv[2],"%lf",&mu);
    sscanf(argv[3],"%lf",&mean_var);
  }
  else {
    fprintf(stderr," enter rho mu mean_var: ");
    fgets(line,sizeof(line),stdin);
    sscanf(line,"%lf %lf %lf",&rho, &mu, &mean_var);
  }

  while (1) {
    fprintf(stderr," enter score length db_size: ");
    if (fgets(line,sizeof(line),stdin)==NULL) exit(0);
    if (line[0]=='\n') exit(0);
    sscanf(line,"%d %d %d",&score, &length, &db_size);
    if (db_size < 1) db_size = 50000;

    z_val = s_to_zv(score, length);

    printf(" s: %d (%d) E(%d): %4.2g\n",score,length,db_size,zv_to_E(z_val,db_size));
  }
}

double s_to_zv(int score, int length)
{
  return ((double)score - rho * log((double)length) - mu)/sqrt(mean_var);
}

double zv_to_E(double zv, int db_size)
{
  double e;

  e = exp(-1.282554983 * zv - .577216);
  return (double)db_size * (e > .01 ? 1.0 - exp(-e) : e);
}
