
/* $Id: param.h 1233 2013-10-08 18:26:31Z wrp $ */
/* $Revision: 1233 $  */

#include <sys/types.h>

#ifndef P_STRUCT
#define P_STRUCT

#define MAXSQ 60

/* Concurrent read version */

struct fastr {
  int ktup;
  int cgap;
  int pgap;
  int pamfact;
  int scfact;
  /* these values will soon be abandoned */
  int bestoff;
  int bestscale;
  int bkfact;
  int bktup;
  int bestmax;
  /* statistics based scaling values */
  int use_E_thresholds;
  double E_join, E_band_opt;
  int altflag;
  int optflag;
  int iniflag;
  int optcut;
  int optcut_set;
  int optwid;
  int optwid_set;
};

struct prostr {
    int gopen;
    int gextend;
    int width;
};

/* must be identical in thr_bufs.h */
struct score_count_s {
  long s_cnt[3];
  long tot_scores;
};

struct pstruct		/* parameters */
{
  int n0;	/* length of query sequence, used for statistics */
  int gdelval;	/* value gap open (-10) */
  int ggapval;	/* value for additional residues in gap (-2) */
  int gshift;	/* frameshift for fastx, fasty */
  int gsubs;	/* nt substitution in fasty */
  int p_d_mat;	/* dna match penalty */
  int p_d_mis;	/* dna mismatch penalty */
  int p_d_set;	/* using match/mismatch */
  int n1_low;
  int n1_high;	/* sequence length limits */
  int score_ix;	/* index to sorted score */
  int show_ident;	/* flag - show identical lalign alignment */
  int nseq;	/* number of different sequences (for lalign) */
  int zsflag;	/* use scalebest() */
  int zsflag2;	/* statistics for best shuffle */
  int zsflag_f;	/* use scalebest() */
  int zs_win;	/* window shuffle size */
  int shuffle_dna3;	/* shuffle dna as codons */
  int histint;		/* histogram interval */
  unsigned char sq[MAXSQ+1];
  int hsq[MAXSQ+1];
  int nsq;		/* length of normal sq */
  /* int pamh1[MAXSQ+1]; */	/* identical match score (diagonal scores) */
  /* int *pamh2[MAXSQ+1]; */	/* ktup match score */
  int ext_sq_set;	/* flag for using extended alphabet */
  unsigned char sqx[MAXSQ+1];
  int hsqx[MAXSQ+1];
  int c_nt[MAXSQ+1];
  int nsqx;	/* length of extended sq */
  int nsq_e;	/* effective nsq */
  int dnaseq;	/* -1 = not set (protein); 0 = protein; 1 = DNA; 2 = other, 3 RNA */
  int nt_align;	/* DNA/RNA alignment = 1 */
  int debug_lib;
  int tr_type;	/* codon table */
  int sw_flag;
  char pamfile[MAX_FN];	/* pam file name */
  char pamfile_save[MAX_FN];  /* original pam file */
  char pam_name[MAX_FN];
  char pgpfile[MAX_FN];
  int pgpfile_type;
  float pamscale;	/* ln(2)/3 or ln(2)/2 */
  float ulambda;	/* ungapped lambda */
  float entropy;	/* bits/position */
  float tfract_id;	/* target fraction id */
  int pam_pssm;
  int pam_set;
  int pam_variable;
  int have_pam2;
  int **pam2[2];	/* set of 2D scoring matrices; [0] lower-case 'x', [1] upper/lower case */
  int **pam2p[2];
  int pamoff;	/* offset for pam values */
  int pam_l, pam_h, pam_xx, pam_xm;	/* lowest, highest pam value */
  int pam_x_set;
  int pam_x_id_sim;	/* =0 -> 'N,X' identical but not similar;
			   =1 -> 'N,X' identical+similar;
			   <0 -> 'N,X' not identical, not similar */
  int pam_ms;		/* use a Mass Spec pam matrix */
  void *fp_struct;	/* function specific parameters based on algorith/scoring matrix */
  int LK_set;
  double pLambda, pK, pH;	/* Karlin-Altscul parameters */
  int maxlen;
  int max_repeat;	/* used for repeat count in ssearch34/lalign */
  int repeat_thresh;
  char *other_info;
  double e_cut;		/* cutoff for scores */
  double e_cut_r; 	/* cutoff for multiple local alignments */
  double zs_off;	/* z-score offset from sampling */
  int do_rep;		/* enable multiple alignments */
  int can_pre_align;	/* flag for have_ares & 0x1 pre-alignments */
  long zdb_size; 	/* force database size */
  int zdb_size_set;	/* flag for user -Z */
  int pgm_id;
  int pseudocts;
  int shuff_node;
  union {
    struct fastr fa;
    struct prostr pr;
  } param_u;
};

#include "rstruct.h"

/* the seq_record has all the invariant data about a sequence -
   sequence length, libstr, sequence itself, etc.
   it does not have the results information
   we can have 1, 2, or 6 (obsolete tfasta) results records for a sequence,
   but there will still be only one sequence record.
*/

struct annot_str {
  /* information for conventional annotations */
  unsigned char *aa1_ann;	/* annotation string */
  /* information for "rich" annotations */
  int n_annot;	/* length of ann_arr_str array */
  int n_domains; 	/* length of domain_arr_p array */
  struct annot_entry *annot_arr_p;	/* array[n_annot] of annot_entry's for all annotations */
  struct annot_entry **s_annot_arr_p;	/* sorted version of annots */
};

/* ann_str keeps information on "rich" annotations, position, type, value */
struct annot_entry {
  long pos;
  long end;
  char label;	/* currently -V *#%!@ symbols, plus 'V' for variant */
  unsigned char value;	/* must be amino acid residue, binary encoded */
  char *comment;
  int target;	 /* 0 for query/ 1 for library */
};

/* domain_str keeps information on "rich" annotations, position, type, value */
struct domfeat_data {
  struct annot_entry *annot_entry_p;
  struct domfeat_data *next;
  long pos;	/* annotation position */
  long a_pos;	/* aligned annotation position */
  long end_pos;	/* domain annotation end */
  int score;	/* score of current region */
  int n_ident;	/* count for percent id */
  int n_alen; 	/* align len for percent id */
  int n_gaplen;	/* number of gap residues in alignment */
};

/* seq_record has the data required to do a calculation */
struct seq_record {
  unsigned char *aa1b;		/* sequence buffer */
  struct annot_str *annot_p;
  int n1;
  long l_offset;	/* q_offset/l_offset set outside getlib() based on chunks; 0-based */
  long l_off;		/* q_off/l_off comes from @C:123, and is 1-based */
  int index;		/* index in search */
#ifdef DEBUG
  long adler32_crc;
#endif
};

/* mseq_record has meta data not required to calculate score or alignment */
struct mseq_record {
  int *n1tot_p;
#ifdef USE_FSEEKO
  off_t lseek;
#else
  long lseek;
#endif
  struct lmf_str *m_file_p;
  int cont;
  char libstr[MAX_UID];
  char *bline;
  int bline_max;
  int annot_req_flag;
  int index		/* index in search */;
#ifdef DEBUG
  long adler32_crc;
#endif
};

struct seqr_chain {
  struct seq_record *seqr_base;
  struct mseq_record *mseqr_base;
  struct seqr_chain *next;
  /*   struct lib_seq_info *ldb_info; */
  int max_chain_seqs;
  int cur_seq_cnt;
  unsigned char *aa1b_base;
  int aa1b_size;
  int aa1b_next;
  int contiguous;
};

struct getlib_str {
  int lcont;		/* lcont save */
  int ocont;		/* ocont save */
  int eof;		/* done with this file */
#ifdef USE_FSEEKO
  off_t lseek;
#else
  long lseek;
#endif
  long loffset;		/* loffset save */
  char libstr[MAX_UID];	/* repository for libstr */
  int n_libstr;		/* length of libstr */
  unsigned char *aa1save;	/* overlapping sequence save */
  struct lib_struct *lib_list_p;
  int *n1tot_ptr, *n1tot_cur;
  int n1tot_cnt;
  int n1tot_v;
  long tot_memK;	/* cummulative amount of memory allocated for aa1b;
			   used to limit memory use */
  long max_memK;	/* allow separate memory limits for main,link
			   searches */
  long lost_memK;	/* check for waste */
  struct seqr_chain *start_seqr_chain;
  struct seqr_chain *cur_seqr_chain;
  int use_memory;
};

#endif	/* P_STRUCT */

#ifndef A_STRUCT
#include "aln_structs.h"
#endif
