/* Concurrent read version */

/* $Id: defs.h 1261 2014-06-11 19:38:36Z wrp $ */
/* $Revision: 1261 $  */

#ifdef SUNOS
#include <sys/stdtypes.h>
#endif

#ifndef IS_BIG_ENDIAN
#if defined(__BIG_ENDIAN__) || defined(_BIG_ENDIAN)
#define IS_BIG_ENDIAN
#else
#undef IS_BIG_ENDIAN
#endif
#endif

#if !defined(MAX_WORKERS) && !defined(PCOMPLIB)
#define MAX_WORKERS 1
#endif
#if defined(PCOMPLIB) && !defined(MAXWRKR)
#define MAXWRKR	64
#endif

#define SAFE_STRNCPY(dest,src,dest_len) strncpy(dest,src,dest_len); dest[dest_len-1]='\0'
#define SAFE_STRNCAT(str,cat,str_len) strncat(str,cat,str_len-strlen(str)-1)

/* constants associated with displaying annotation links */
#ifndef DESCR_OFFSET
#define DESCR_OFFSET 20
#endif

#define NO_FILE_EXIT 4


/* 3-Oct-2003 - we can now have 2 nucleotide query types, DNA
   and RNA.  pst.dnaseq can also be SEQT_RNA.
   ldnaseq can only be DNA */

#define SEQT_DNA 1
#define SEQT_RNA 3	/* DNA and RNA seqtypes must be odd */

#define SEQT_PROT 0
#define SEQT_UNK -1
#define SEQT_OTHER 2

#ifndef DEF_NMLEN
#define DEF_NMLEN 6
#endif

#define DEF_MIN_BITS 40	/* minimum number of bits required, appropriate for swissprot */

/* unfortunately, there is an important relationship between MAXTRN and
   MAXLIB embedded here.  MAXTRN must be >= (MAXLIB)/3
   or it will be possible for a translated DNA sequence to be longer
   than the translation space available */

#define MAX_STR	512 /* standard label/message buffer */
#define MAX_SSTR 32 /* short string */
#define MAX_LSTR 4096 /* long label/message buffer */
#define MAX_FN  120 /* maximum size of a file name */
#define MAX_CH	40 /* maximum number of library choices */
#ifndef SMALLMEM
#define MAX_LF  2000 /* maximum numer of library files */
#else
#define MAX_LF  80 /* maximum numer of library files */
#endif

#ifndef MAX_MEMK
#if defined(BIG_LIB64) && (defined(COMP_THR) || defined(PCOMPLIB))
#define MAX_MEMK 16*1024*1024	/* 16 GB (<<10) for library in memory */
#else
#define MAX_MEMK 2*1024*1024	/* 2 GB (<<10) for library in memory */
#endif
#endif

/* padding at the end of sequences for ALTIVEC, other vector
   processors */
#define SEQ_PAD 16

#define MAX_UID 20 /* length of libstr, used for character keys with SQL */

#define BUF_MULT 2	/* increase to increase the number of buffers */
#define DEF_WORKER_BUF 6000000
#define AVE_AA_LEN 400
#define AVE_NT_LEN 1200

#define MAX_RSTATS 500		/* number of random shuffle stats */
#define MIN_LOCAL_LEN 33	/* minimum length for addn'l local alignments
				   (should be in pstruct)*/
#ifndef SMALLMEM
#define MAXTST	40000		/* longest query */
#define MAXLIB	150000		/* longest library sequence*/
#define MAXLIB_P 45000
#define MIN_RES 2000		/* minimum amount allocated for alignment */
#ifndef TFAST
#define MAXTRN  45000		/* buffer for fastx translation */
#else
#define MAXTRN 165000		/* buffer for tfastx translation, must be > 3 * MAXTST */
#endif
#define SEQDUP	150		/* future - overlap */
#ifndef PCOMPLIB
#ifndef MAX_BEST
#define MAX_BEST	60000	/* max number of best scores */
#endif
#define MAX_STATS 60000
#else
#ifndef MAX_BEST
#define MAX_BEST	60000	/* max number of best scores */
#endif
#define MAX_STATS 60000
#endif
#define BIGNUM  1000000000
#ifndef MAXINT
#define MAXINT 2147483647
#endif
#define MAXLN	120	/* size of a library name */
#else
#define MAXTST	1500
#define MAXLIB	10000
#define MAXLIB_P MAXLIB
#define MIN_RES 1000
#ifndef TFAST
#define MAXTRN  4500
#else
#define MAXTRN 11500
#endif
#define SEQDUP	300
#define MAX_BEST 2000
#define MAX_STATS 20000
#define BIGNUM  32767
#define MAXINT  32767
#define MAXLN	40	/* size of a library name */
#endif
#if !defined(TFAST)
#define MAXDIAG	(MAXTST+MAXLIB)
#else
#define MAXDIAG	(MAXTST+MAXTRN)
#endif

#define MAXPAM	600	/* maximum allowable size of the pam matrix */
#define PROF_MAX 500
#define ALF_MAX 30

#define max(a,b) (((a) > (b)) ? (a) : (b))
#define min(a,b) (((a) < (b)) ? (a) : (b))

#define MX_ATYPE 7	/* markx==0,1,2 7=> no alignment */
#define MX_ASEP  8	/* markx==3  - separate lines */
#define MX_AMAP  16	/* markx==4,5 - graphic map */
#define MX_HTML  32	/* markx==6  - HTML */
#define MX_M9SUMM 64	/* markx==9(c) */
#define MX_M10FORM 128	/* markx==10 - verbose output */
#define MX_M11OUT 256	/* markx==11 - lalign lav */
#define MX_M8OUT 512	/* markx==8 blast tabular (-outfmt=6) output */
#define MX_M8COMMENT 1024	/* markx==8 blast tabular (-outfmt=7) with comments  output */
#define MX_MBLAST 2048	/* markx=B blast alignment -outfmt=0 output */
#define MX_MBLAST2 4096	/* markx=BB blast best scores and alignment (-outfmt=0) output */
#define MX_ANNOT_COORD 16384 /* -m 0, use -m 0B for both */
#define MX_ANNOT_MID  32768 /* markx 0M, 1M, 2M annotations in middle */
#define MX_RES_ALIGN_SCORE (1<<20)  /* show residue alignment score, not alignment */
#define MX_M8_BTAB_LEN  (1<<21) /* show query/subject seq. lens in -m 8 output */

/* codes for -m 9, -m 8C? */
#define SHOW_CODE_ID	1	/* identity only */
#define SHOW_CODE_IDD   2	/* identity with domains */
#define SHOW_CODE_ALIGN 4	/* encoded alignment */
#define SHOW_CODE_CIGAR 8	/* CIGAR vs old encoded alignment */
#define SHOW_CODE_BTOP	16	/* BLAST BTOP encoding */
#define SHOW_CODE_MASK  12	/* use higher bits for annotation format */
#define SHOW_CODE_EXT   16	/* encode identity, mismatch state */
#define SHOW_ANNOT_FULL 32	/* show full-length annot in calc_code */
#define SHOW_CODE_DOMINFO 64    /* include raw domain info in btab/BTOP */

