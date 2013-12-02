/* Concurrent read version */

/* $Id: mw.h 625 2011-03-23 17:21:38Z wrp $ */
/* $Revision: 625 $  */

#include "param.h"

#ifndef FSEEK_T_DEF
#ifndef USE_FSEEKO
typedef long fseek_t;
#else
typedef off_t fseek_t;
#endif
#endif

struct beststr {
  struct seq_record *seq;	/* sequence info */
  struct mseq_record *mseq;	/* sequence meta-info */
  struct beststr *bbp_link;	/* link to a previous beststr entry with the same sequence */
  struct rstruct rst;		/* results info */

  int n1;		/* duplicate of seq.n1, used for error checking/debugging */
#ifdef DEBUG
  long adler32_crc;	/* duplicate of seq.adler32_crc for error checking/debugging */
#endif
  int frame;		/* in buf2_str */
  int repeat_thresh;	/* threshold for additional alignments */
  double zscore;	/* the z-score mostly exists for sorting best scores */
  double zscore2;	/* z-score - from high-scoring shuffles  */
  double bit_score;	/* move to bit-scores for consistency */
  double bit_score2;	/* bit-score for second shuffle */

  int a_res_cnt;
  struct a_res_str *a_res;	/* need only a_res, not a_res[2], because different frames
				   for the same sequence are stored separately */
  int have_ares;
  float percent, gpercent;
};

struct stat_str {
  int score;
  int n1;
  double comp;
  double H;
  double escore;
  int segnum;
  int seglen;
};

