/*  mm_file.h - defines m_file_str for mmap()ed files  */

/* $Id: mm_file.h 938 2012-06-04 16:15:06Z wrp $ */

/* copyright (c) 1999, 2014 by William R. Pearson and The Rector &
   Vistors of the University of Virginia */


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

#include <sys/types.h>

#ifndef FSEEK_T_DEF
#ifndef USE_FSEEKO
#define FSEEK_T_DEF
#define FSEEK fseek
#define FTELL ftell
typedef long fseek_t;
#else
#define FSEEK fseeko
#define FTELL ftello
typedef off_t fseek_t;
#endif
#endif

#ifdef HAS_INTTYPES
#include <inttypes.h>
#else
#ifdef WIN32
typedef __int64 int64_t;
typedef unsigned __int64 uint64_t;
#else
typedef long int64_t;
typedef unsigned long uint64_t;
#endif
#endif
typedef int64_t MM_OFF;

#ifdef MYSQL_DB
#include <mysql.h>
#endif
#ifdef PGSQL_DB
#include <libpq-fe.h>
#endif

#ifndef MAX_FN
#include "defs.h"
#endif

extern unsigned long adler32();

struct lmf_str {
  FILE *libf;		/* sequence file being read */
  char *lb_name;	/* file name */
  char opt_text[MAX_FN];	  /* text after filename */
  int lb_type;		/* library type */
  int *sascii;		/* ascii -> sq mapping */
  int *vascii; 		/* annotation to ann mapping */

  char *annot_sname;	/* annotation script name */

  /* used by flat files */
  char *lline;		/* last line read */
  int acc_off;		/* start of libstr (+1 for agetlib/fasta) */
  unsigned char *cpsave;	/* position in line for lgetlib() */
  fseek_t lpos;			/* position in file */

  /* blast2.0 stuff */
  FILE *hfile;		/* BLAST2.0 description file */
  int bl_format_ver;	/* blast formatdb version */
  int bl_lib_pos;	/* for ncbl2 */
  int pref_db;		/* preferred database */
  int have_oid_list;	/* we have an oid file, must read oid's */	
  unsigned int *oid_list;	/* oid list for subsets */
  int oid_seqs;		/* start offset for mask array */
  unsigned int max_oid;	/* start offset for mask array */

  /* Genbank Flat files */
  int lfflag;		/* flag for CRLF in EMBL CDROM files */

  /* stuff for GCG format files (5,6) */
  int gcg_binary;	/* flag for binary gcg format */
  long gcg_len;		/* length of GCG sequence */

  /* used when memory mapping */
  int mm_flg;		/* mmap worked */
  int mmap_fd;		/* mmap_fd */
  char *mmap_base;	/* base */
  char *mmap_addr;	/* current pos */
  long st_size;		/* file size */

  MM_OFF *d_pos_arr;	/* pointer to desc. offsets */
  MM_OFF *s_pos_arr;	/* pointer to seq. offsets */
  MM_OFF *b_pos_arr;	/* pointer to binary seq. offsets */
  MM_OFF *a_pos_arr;	/* pointer to aux offsets */

  /* currently available only for memory mapped files */
  int max_cnt;		/* # database entries */
  int64_t tot_len;	/* total residue length */
  long max_len;		/* maximum sequence lengh */
  long maxn;		/* maximum possible length */
  long mdup;		/* duplication for overlapping sequences */
  int lib_aa;		/* 0 = DNA, 1 = prot */
  char *tmp_buf;	/* temporary buffer */
  int tmp_buf_max;	/* max size */
  int (*sel_acc_p)(char *, int gi, void *);	/* used to select subset of library */
  void *sel_local;				/* local data structure for sel_acc_p() */

  /* used for SQL database queries */
  char *sql_db, *sql_query, *sql_getdesc, *sql_getseq, *sql_close_tables;
  int sql_reopen;
  char **sql_uid_arr;	/* indexed by lpos */
  /* used to get sequence data */
  char *sql_seqp;

#ifdef MYSQL_DB
  /* used to open the database */
  MYSQL *mysql_conn;
  MYSQL_RES *mysql_res;
  MYSQL_ROW mysql_row;
#endif

#ifdef PGSQL_DB
  /* used to open the database */
  PGconn *pgsql_conn;
  PGresult *pgsql_res;
#endif

  int (*getlib)(unsigned char *seq, int maxs, char *ann,
		int n_libstr,
		fseek_t *libpos,
		int *lcont,
		struct lmf_str *lm_fd,
		long *l_off);

  void (*ranlib)(char *str, int cnt,
		 fseek_t libpos, char *libstr,
		 struct lmf_str *lm_fd);

  int (*get_mmap_chain)(struct seqr_chain *cur_seqr_chain,
			struct lmf_str *m_fd, struct db_str *db);
};
