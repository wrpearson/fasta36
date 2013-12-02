/* $Id: dropgsw2.h 625 2011-03-23 17:21:38Z wrp $ */
/* $Revision: 625 $  */

/* global definitions shared by dropgsw.c and altivec.c */

/* definitions for SW */

struct f_struct {
  struct swstr *ss;
  int *waa_s, *waa_a;
  int **pam2p[2];
  int max_res;
  double aa0_f[MAXSQ];
  double *kar_p;
  double e_cut;
  int show_ident;
  int max_repeat;
#if defined(SW_ALTIVEC) || defined(SW_SSE2)
  unsigned char      bias;
  unsigned char      ceiling;
  unsigned short *   word_score;
  unsigned char *    byte_score;
  void *             workspace;
  int                alphabet_size;
  void *             word_score_memory;
  void *             byte_score_memory;
  void *             workspace_memory;
  int                try_8bit;
  int                done_8bit;
  int                done_16bit;
#endif
};

#ifdef LALIGN
void SIM(const unsigned char *A, /* seq1 indexed A[1..M] */
	 const unsigned char *B, /* seq2 indexed B[1..N] */
	 int M, int N,		 /* len seq1, seq2 */
	 struct pstruct *ppst,	/* parameters */
	 int nseq,		 /* nseq - number of different sequences */
	 int mini_score,	 /* cut-off score */
	 int max_count,		 /* number of alignments */
	 struct a_res_str *a_res);	/* alignment result structure */

int same_seq(const unsigned char *aa0, int n0,
	     const unsigned char *aa1, int n1);
#endif
