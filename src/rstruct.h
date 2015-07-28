/* $Id: rstruct.h 625 2011-03-23 17:21:38Z wrp $ */

#ifndef RSTRUCT
#define RSTRUCT
struct rstruct
{
  int score[3];
  int valid_stat;
  int alg_info;
  double comp;
  double H;
  double escore;
  int segnum;
  int seglen;
};
#endif
