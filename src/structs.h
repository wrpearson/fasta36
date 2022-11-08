/* $Id: structs.h 1259 2014-05-28 19:28:06Z wrp $ */

#include "rstruct.h"

#include "aln_structs.h"

#include "param.h"

struct hist_str {
  int histflg;
  int *hist_a;
  int histint, min_hist, max_hist, maxh;
  long entries;
  int z_calls;
  char stat_info[MAX_STR];
};

struct db_str {
  long entries;
  unsigned long length;
  int carry;
};

struct mng_thr {	/* structure to keep track of thread buffers */
  int max_work_buf;	/* number of threads buffers */
  int max_buf2_res;	/* number of results/thread buffer */
  int max_chain_seqs;	/* number of sequences/seqr chain */
  int seq_buf_size;	/* size of sequence buffer in each max_seq_cnt chain */
};

/* values required for successive calles to getlib() */
struct lib_seq_info {
  int ldnaseq;
  int term_code;
  int maxn;
  int dupn;
  int l_overlap;
  int maxt3;
};

/* used by help system to describe optional arguments */
struct opt_def_str {
  char opt_char;
  int has_arg;
  char *opt_str;
  char *opt_descr_s;
  char *opt_descr_l;
  int opt_rank;
  int fmt_type;
  int i_param1;
  int i_param2;
  double d_param1;
  double d_param2;
  char *s_param;
};

struct markx_str {
  /* values to be copied into m_msg to modify format */
  int markx;
  int nohist;
  int ashow;
  int show_code;
  int long_info;
  int aln_llen;
  int aln_llcntx;
  int aln_llcntx_set;
  int std_output;

  char *out_file;
  FILE *out_fd;
  struct markx_str *next;
};

struct mngmsg 		/* Message from host to manager */
{
  char pgm_name[MAX_FN];	/* program name from argv[0] */
  long max_memK;	/* maximum amount of sequence buffer memory */
  int cur_seqr_cnt;	/* current number of seqr chains */
  int n0;		/* query length ^qm_msg */
  int nm0;		/* number of segments ^qm_msg */
  int nmoff;		/* length of fastf segment */
  unsigned char *aa0a;	/* annotation array */
  struct annot_str *annot_p;	/* annot_str for query */
  unsigned char ann_arr[MAX_STR]; /* annotation characters */
  int ann_arr_n;	/* number of annotation characters */
  char *ann_arr_def[MAX_STR];	/* definitions of ann_arr characters */
  int ann_flg;		/* have annotation array, characters */
  int have_ann;		/* have annotation on this query */
  char tname[MAX_FN];	/* Query sequence name */
  int tnamesize;	/* Query name size */
  char lname[MAX_LSTR];	/* Library  file  name */
  char link_lname[MAX_LSTR]; /* link-library name */
  char annot0_sname[MAX_LSTR]; /* query annotation script name */
  char annot1_sname[MAX_LSTR]; /* library annotation script name */
  int max_tot;		/* function defined total sequence area */
  int qdnaseq;		/* query is protein (0)/dna (1) */
  int q_overlap;	/* overlap when segmenting long query sequence */
  long q_offset;	/* q_offset, l_offset are set outside getlib()
			   and report the number of residues that were
			   skipped/read-previously to get to this
			   version of aa0/aa1 (0-based) */
  long q_off;		/* q_off, l_off are set by getlib(); starting
			   at 1 but modified by @C: the position of
			   the first residue of aa0/aa1 in the
			   original sequence */
  int qframe;		/* number of possible query frames */
  int nframe;		/* frame for TFASTA */
  int nitt1;		/* nframe-1 */
  int thr_fact;		/* fudge factor for threads */
  int s_int;		/* sampling interval for statistics */
  int ql_start;		/* starting query sequence */
  int ql_stop;		/* ending query sequence */
  int pbuf_siz;		/* buffer size for sequences send in p2_complib */
  char qtitle[MAX_STR];	/* query title */
  char ltitle[MAX_STR];	/* library title */
  char flstr[MAX_FN];	/* FASTLIBS string */
  int std_output;	/* produce normal output */
  char outfile[MAX_FN];
  FILE *outfd;	/* indicates outfile is already open */
  char label [MAX_SSTR];	/* Output label, "opt", "s-w", "initn init1", "initn opt"  */
  char alabel[MAX_SSTR];	/* Output label, "Smith-Waterman", "banded Smith-Waterman", etc  */
  char f_id0[4];	/* function id for markx==10 */
  char f_id1[4];	/* function id for markx==10 */
  char sqnam[4];	/* "aa" or "nt" */ 
  char sqtype[10];	/* "DNA" or "protein" */
  int long_info;	/* long description flag*/
  int blast_ident;	/* calculate identities excluding gaps */
  long sq0off, sq1off;	/* virtual offset into aa0, aa1 */
  int markx;		/* alignment display type */
  int tot_markx;	/* markx as summ of all alternative markx */
  struct markx_str *markx_list;	/* list of alternate outputs */
  int nbr_seq;		/* number of library sequences */
  int n1_high;		/* upper limit on sequence length */
  int n1_low;		/* lower limit on sequence length */
  double e_cut;		/* e_value for display */
  double e_low;		/* e_value for display */
  int e_cut_set;	/* e_value deliberately set */
  int zsflag;		/* zsflag for all searches */
  int zsflag2;		/* zsflag2 for all searches */
  int pamd1;		/* 1st dimension of pam matrix */
  int pamd2;		/* 2nd dimension of pam matrix */
  int revcomp;		/* flag to do reverse complement */
  int quiet;		/* quiet option */
  int nrelv;		/* number of interesting scores */
  int srelv;		/* number of scores to show in showbest */
  int arelv;		/* number of scores to show at alignment */
  int z_bits;		/* z_bits==1: show bit score, ==0 show z-score */
  int tot_ident;	/* tot_ident=1 -> no mismatches for 100% identity */
  int gi_save;		/* do not remove gi|12345 in displays */
  char alab[3][24];	/* labels for alignment scores */
  int nohist;		/* no histogram option */
  int do_showbest;	/* do not showbest() */
  int nshow;		/* number shown in showbest() */
  int nskip;		/* number skipped by e_low */
  int mshow;		/* number of scores to show */
  int mshow_set;	/* mshow set with -b */
  int mshow_min;	/* at least mshow scores must be shown, but limited by e_cut */
  int ashow;		/* number of alignments to show */
  int ashow_set;	/* ashow set with -d */
  int nmlen;		/* length of name label */
  int show_code;	/* show alignment code in -m 9;  ==1 => identity only, ==2 alignment code*/
  int m8_show_annot;	/* show annotations only in -m 8CB output */
  int tot_show_code;	/* show alignment for all outputs */
  int pre_load_done;	/* set after pre_load_best() call */
  int align_done;	/* do_walign() called */
  unsigned char *aa1save_buf_b; /* buffer for re_getlib() sequences for alignments */
  char *bline_buf_b; /* buffer for re_getlib() sequences for alignments */
  int thold;		/* threshold */
  int max_repeat;	/* max number of non-intersecting local alignments */
  int last_calc_flg;	/* needs a last calculation stage */
  int qshuffle;		/* shuffle the query and do additional comparisons ^qm_msg*/
  int shuff_max;	/* number of shuffles to perform */
  int shuff_max_save;	/* number of shuffles to perform */
  int shuff_node;	/* number of shuffles/worker node */
  int shuff_wid;
  int stages;		/* number of stages */
  double Lambda, K, H;	/* (ungapped) Karlin-Altschul parameters */
  int escore_flg;	/* use escore calculated by do_work() */
  struct lib_seq_info ldb_info;	/* maxn, maxt, l_overlap, ldnaseq */
  struct db_str db;	/* the database size as read */
  struct db_str ldb;	/* the database size via save_best() */
  struct score_count_s s_info;	/* counts of different score types */
  struct score_count_s ss_info;	/* counts of different score types from shuffles */
  struct hist_str hist;
  void *pstat_void;	/* pstat structure for standard statistics */
  void *pstat_void2;	/* pstat structure for zsflag > 2 */
  struct a_struct aln;	/* has llen, llnctx, llnctx_flg, showall */
  struct a_res_str a_res; /* has individual alignment coordinates */
  char dfile [MAX_FN];	/* file for dumping scores to */
};

struct lib_struct {
  char *file_name; /* this library file */
  int lib_type;	   /* the library type can be specified here, for files in an indirect list */
  /* struct lib_mol_info *lib_params; */ /* parameters (ldnaseq, term_code, constant for all libraries */
  struct lmf_str *acc_file_p;
  struct lmf_str *m_file_p;	/* magic *m_file_p for reading */
  struct lib_struct *next;	/* next in the list */  
};
