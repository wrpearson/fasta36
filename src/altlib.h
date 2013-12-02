
/* $Id: altlib.h 905 2012-01-30 17:33:06Z wrp $ */
/* $Revision: 905 $  */

/* #ifdef UNIX */
/* ncbi blast 1.3 format */
/*
#define NCBIBL13 11
extern int ncbl_getliba();
extern void ncbl_ranlib();
void ncbl_closelib();
*/
#define NCBIBL20 12
/* #endif */

#ifdef MYSQL_DB
#define MYSQL_LIB 16
#define LASTLIB MYSQL_LIB+1
#endif

#ifdef PGSQL_DB
#define PGSQL_LIB 17
#define LASTLIB PGSQL_LIB+1
#endif

#if !defined (LASTLIB) && defined(NCBIBL20)
#define LASTLIB NCBIBL20+1
#endif
#if !defined (LASTLIB)
#define LASTLIB 10
#endif

#define FASTA_F 0
#define DEFAULT 0
#define FULLGB 1
#define UNIXPIR 2
#define EMBLSWISS 3
#define INTELLIG 4
#define VMSPIR 5
#define GCGBIN 6
#define FASTQ 7
#define LASTTXT 7
#define ACC_LIST 10

#include "mm_file.h"

/* pearson fasta format */
int agetlib(unsigned char *, int, char *, int, fseek_t *, int *,
	    struct lmf_str *, long *);
void aranlib(char *, int, fseek_t, char *, struct lmf_str *);
/* full uncompressed GB FULLGB*/
extern int lgetlib(unsigned char *, int, char *, int, fseek_t *, int *,
	    struct lmf_str *, long *);
extern void lranlib(char *, int, fseek_t, char *, struct lmf_str *);
/* PIR UNIX protein UNIXPIR */
extern int pgetlib(unsigned char *, int, char *, int, fseek_t *, int *,
	    struct lmf_str *, long *);
extern void pranlib(char *, int, fseek_t, char *, struct lmf_str *);
/* EMBL/SWISS-PROT EMBLSWISS */
extern int egetlib(unsigned char *, int, char *, int, fseek_t *, int *,
	    struct lmf_str *, long *);
extern void eranlib(char *, int, fseek_t, char *, struct lmf_str *);

/* Intelligenetics INTELLIG */
extern int igetlib(unsigned char *, int, char *, int, fseek_t *, int *,
	    struct lmf_str *, long *);
extern void iranlib(char *, int, fseek_t, char *, struct lmf_str *);
/* PIR VMS format */
extern int vgetlib(unsigned char *, int, char *, int, fseek_t *, int *,
	    struct lmf_str *, long *);
extern void vranlib(char *, int, fseek_t, char *, struct lmf_str *);
/* GCG 2bit format */
extern int gcg_getlib(unsigned char *, int, char *, int, fseek_t *, int *,
	    struct lmf_str *, long *);
extern void gcg_ranlib(char *, int, fseek_t, char *, struct lmf_str *);

/* FASTQ format (ignoring quality scores) */
int qgetlib(unsigned char *, int, char *, int, fseek_t *, int *,
	    struct lmf_str *, long *);
void qranlib(char *, int, fseek_t, char *, struct lmf_str *);

#ifdef NCBIBL20
/* ncbi blast 2.0 format */
extern int ncbl2_getliba(unsigned char *, int, char *, int, fseek_t *, int *,
	    struct lmf_str *, long *);
extern void ncbl2_ranlib(char *, int, fseek_t, char *, struct lmf_str *);
void ncbl2_closelib();
#endif

#ifdef MYSQL_DB
extern int mysql_getlib(unsigned char *, int, char *, int, fseek_t *, int *,
	    struct lmf_str *, long *);
extern void mysql_ranlib(char *, int, fseek_t, char *, struct lmf_str *);
int mysql_closelib();
#endif

int (*getliba[LASTLIB])(unsigned char *, int, char *, int, fseek_t *, int *,
	    struct lmf_str *, long *)={
  agetlib,lgetlib,pgetlib,egetlib,
  igetlib,vgetlib,gcg_getlib,qgetlib,
  agetlib,agetlib
#ifdef UNIX
  ,agetlib
#ifdef NCBIBL13
  ,ncbl_getliba
#else
  ,ncbl2_getliba
#endif
#ifdef NCBIBL20
  ,ncbl2_getliba
#endif
#ifdef MYSQL_DB
  ,agetlib
  ,agetlib
  ,agetlib
  ,mysql_getlib
#endif
#endif
};

void (*ranliba[LASTLIB])(char *, int, fseek_t, char *, struct lmf_str *)={
  aranlib,lranlib,pranlib,eranlib,
  iranlib,vranlib,gcg_ranlib,qranlib,
  aranlib,aranlib
#ifdef UNIX
  ,aranlib
#ifdef NCBIBL13
  ,ncbl_ranlib
#else
  ,ncbl2_ranlib
#endif
#ifdef NCBIBL20
  ,ncbl2_ranlib
#endif
#ifdef MYSQL_DB
  ,aranlib
  ,aranlib
  ,aranlib
  ,mysql_ranlib
#endif
#endif
};
