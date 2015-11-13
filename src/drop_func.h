/* drop_func.h */

/* $Id: drop_func.h 1196 2013-07-19 20:18:21Z wrp $ */
/* $Revision: 1196 $  */

/* copyright (c) 2005, 2014 by William R. Pearson and The Rector & Vistors
   of the University of Virginia */

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

/* functions provided by each of the drop files */

#ifdef DEBUG
unsigned long adler32(unsigned long, const unsigned char *, unsigned int);
#endif

void	/* initializes f_struct **f_arg */
init_work (unsigned char *aa0, int n0,
	   struct pstruct *ppst,
#ifndef DROP_INTERN
	   void **f_arg
#else
	   struct f_struct **f_arg
#endif
);


void	/* frees memory allocated in f_struct */
close_work (const unsigned char *aa0, int n0,
	    struct pstruct *ppst,
#ifndef DROP_INTERN
	   void **f_arg
#else
	   struct f_struct **f_arg
#endif
);

void	/* documents search function, parameters */
get_param (const struct pstruct *pstr, 
	   char **pstring1, char *pstring2,
	   struct score_count_s *);

void	/* calculates alignment score(s), returns them in rst */
do_work (const unsigned char *aa0, int n0,
	 const unsigned char *aa1, int n1,
	 int frame,
	 const struct pstruct *ppst,
#ifndef DROP_INTERN
	 void *f_arg,
#else
	 struct f_struct *f_arg,
#endif
	 int qr_flg, int shuff_flg, struct rstruct *rst,
	 struct score_count_s *);

void	/* calculates optimal alignment score */
do_opt (const unsigned char *aa0, int n0,
	const unsigned char *aa1, int n1,
	int frame,
	struct pstruct *ppst,
#ifndef DROP_INTERN
	void *f_arg,
#else
	struct f_struct *f_arg,
#endif
	struct rstruct *rst
	);

struct a_res_str *	/* produces encoding of alignment */
do_walign (const unsigned char *aa0, int n0,
	   const unsigned char *aa1, int n1,
	   int frame, int repeat_thresh,
	   struct pstruct *ppst, 
#ifndef DROP_INTERN
	   void *f_arg,
#else
	   struct f_struct *f_arg,
#endif
	   int *have_ares);

void
pre_cons(const unsigned char *aa, int n, int frame, 
#ifndef DROP_INTERN
	   void *f_arg
#else
	   struct f_struct *f_arg
#endif
	);

void 
aln_func_vals(int frame, struct a_struct *aln);

#include "dyn_string.h"

/* calc_cons_a - takes aa0, aa1, a_res, and produces seqc0, seqc1, 
 *             and seqc0a, seqc1a - the annotated sequences 
 */
int
calc_cons_a(const unsigned char *aa0, int n0,
	    const unsigned char *aa1, int n1,
	    int *nc,
	    struct a_struct *aln,
	    struct a_res_str *a_res,
	    struct pstruct *ppst,
	    char *seqc0, char *seqc1, char *seqca, int *seqc_score,
	    const unsigned char *ann_arr,
	    const unsigned char *aa0a, const struct annot_str *annot0_p, char *seqc0a, 
	    const unsigned char *aa1a, const struct annot_str *annot1_p, char *seqc1a, 
	    int *score_delta,
	    struct dyn_string_str *annot_var_dyn,
#ifndef DROP_INTERN
	    void *f_arg,
#else
	    struct f_struct *f_arg,
#endif
	    void *pstat_void
	    );

int	/* returns lenc - length of aligment */
calc_code(const unsigned char *aa0, int n0,
	  const unsigned char *aa1, int n1,
	  struct a_struct *aln,
	  struct a_res_str *a_res,
	  struct pstruct *ppst,
	  struct dyn_string_str *align_code_dyn,
	  const unsigned char *ann_arr,
	  const unsigned char *aa0a,
	  const struct annot_str *annot0_p,
	  const unsigned char *aa1a, 
	  const struct annot_str *annot1_p,
	  struct dyn_string_str *ann_str_dyn,
	  int *score_delta,
#ifndef DROP_INTERN
	  void *f_arg,
#else
	  struct f_struct *f_arg,
#endif
	  void *pstat_void,
	  int code_fmt
	  );

int 	/* returns lenc - length of alignment */
calc_id(const unsigned char *aa0, int n0,
	const unsigned char *aa1, int n1,
	struct a_struct *aln, 
	struct a_res_str *a_res,
	struct pstruct *ppst,
	const struct annot_str *annot0_p,
	const struct annot_str *annot1_p,
	int *score_delta,
	struct dyn_string_str *annot_var_dyn,
#ifndef DROP_INTERN
	void *f_arg
#else
	struct f_struct *f_arg
#endif
	);

int 	/* returns lenc - length of alignment */
calc_idd(const unsigned char *aa0, int n0,
	const unsigned char *aa1, int n1,
	struct a_struct *aln, 
	struct a_res_str *a_res,
	struct pstruct *ppst,
	const struct annot_str *annot0_p,
	const struct annot_str *annot1_p,
	int *score_delta,
	struct dyn_string_str *annot_var_dyn,
#ifndef DROP_INTERN
	void *f_arg
#else
	struct f_struct *f_arg
#endif
	);
