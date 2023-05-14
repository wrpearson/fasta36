/* $Id: nmgetlib.c 1251 2014-01-24 21:34:05Z wrp $ */
/* $Revision: 1251 $  */

/*  copyright (c) 1987, 1988, 1989, 1992, 1995, 2000, 2014 by 
    William  R. Pearson and The Rector & Vistors of the University of
    Virginia */

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

/*	May, June 1987	- modified for rapid read of database

	revised (split) version of nmgetaa.c -> renamed nmgetlib.c

	This version seeks to be a thread safe, no global, library
	reading program.  While adjusting the routines in this file
	should be relatively easy, ncbl2_mlib.c and mysql_lib.c may be
	more difficult.

	nmgetlib.c and mmgetaa.c are used together.  nmgetlib.c provides
	the same functions as nxgetaa.c if memory mapping is not used,
	mmgetaa.c provides the database reading functions if memory
	mapping is used. The decision to use memory mapping is made on
	a file-by-file basis.

	June 2, 1987 - added TFASTA
	March 30, 1988 - combined ffgetaa, fgetgb;
	April 8, 1988 - added PIRLIB format for unix
	Feb 4, 1989 - added universal subroutines for libraries
	December, 1995 - added range option file.name:1-1000
	September, 1999 - added option for mmap()ed files using ".xin"
*/

/*
	February 4, 1988 - this starts a major revision of the getaa
	routines.  The goal is to be able to seach the following format
	libraries:

	0 - normal FASTA format
	1 - full Genbank flatfile format
	2 - NBRF/PIR CODATA format
	3 - EMBL/Swiss-prot format
	4 - Intelligentics format
	5 - NBRF/PIR VMS format
	6 - GCG 2bit format
	7 - FASTQ format
	8 - accession script

	10 - list of gi/acc's
	11 - NCBI setdb/blastp (1.3.2) AA/NT
	12 - NCBI setdb/blastp (2.0) AA/NT
	13 - NCBI makeblastdb/blastp (3.0) AA/NT

	16 - mySQL queries
	
	see file altlib.h to confirm numbers

*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>

#ifdef __GNUC__
extern FILE *popen(const char *script, const char *mode);
extern int pclose(FILE *);
#endif

#include "defs.h"
#include "structs.h"

#ifndef SFCHAR
#define SFCHAR ':'
#endif

#define EOSEQ 0

#include "uascii.h"
/* #include "upam.h" */

#define LFCHAR '\015'  /* for MWC 5.5 */

#include "altlib.h"

#include <fcntl.h>
#ifndef O_RAW
#ifdef O_BINARY
#define O_RAW O_BINARY
#else
#define O_RAW 0
#endif		/* O_BINARY */
#endif		/* O_RAW */

#ifdef WIN32
#define RBSTR "rb"	/* read file in binary mode */
#else
#define RBSTR "r"
#endif

char *alloc_file_name(char *f_name);
struct lib_struct *get_lnames(char *tname, struct lib_struct *cur_lib_p);
struct lmf_str *load_mmap(FILE *, char *, int, int, struct lmf_str *);
struct lmf_str *ncbl2_openlib(struct lib_struct *, int ldnaseq);
struct lmf_str *ncbl2_reopen(struct lmf_str *);

static struct lmf_str *last_m_fptr=NULL;

int sel_acc_libstr(char *libstr, int gi, void *ptr);
void *sel_acc_libstr_init(FILE *libf, int *acc_off, char fmt_term);

int sel_acc_gi(char *libstr, int gi, void *ptr);
void *sel_acc_gi_init(FILE *libf, int *acc_off, char fmt_term);

int sel_hacc_libstr(char *libstr, int gi, void *ptr);
void *sel_hacc_libstr_init(FILE *libf, int *acc_off, char fmt_term);

int sel_hacc_gi(char *libstr, int gi, void *ptr);
void *sel_hacc_gi_init(FILE *libf, int *acc_off, char fmt_term);

#define MAX_ACC_TYPE 4
int (*sel_acc_arr[MAX_ACC_TYPE+1])(char *libstr, int gi, void *ptr) = {
  NULL, sel_acc_libstr, sel_acc_gi, sel_hacc_libstr, sel_hacc_gi
};

void *(*sel_acc_init[MAX_ACC_TYPE+1])(FILE *libf, int *acc_off, char fmt_term) = {
  NULL, sel_acc_libstr_init, sel_acc_gi_init, sel_hacc_libstr_init, sel_hacc_gi_init
};

unsigned int hash_func(char *key);
unsigned int fast_hash32 (unsigned int data);

#ifdef MYSQL_DB
struct lmf_str *mysql_openlib(char *, int, int *);
struct lmf_str *mysql_reopen(struct lmf_str *);
#endif

#ifdef PGSQL_DB
struct lmf_str *pgsql_openlib(char *, int, int *);
struct lmf_str *pgsql_reopen(struct lmf_str *);
#endif

extern int can_mmap(int lib_type);

int closelib(struct lmf_str *m_fptr, int force);
extern void newname(char *nname, char *oname, char *suff, int maxn);

/* a file name for openlib may include a library type suffix */

struct lmf_str *
open_lib(struct lib_struct *lib_p, int ldnaseq, int *sascii, int outtty)
{
  struct lmf_str *om_fptr;
  char rline[10], iname[MAX_FN];
  char *bp, *bp1, *bp2;
  char opt_text[MAX_FN];	/* save text after ':' */
  char f_line[MAX_STR];
  int wcnt, opnflg;
  int lib_type;
  int acc_ltype = 1;	/* def type is 1, not zero, so that the acc is read */
  struct lmf_str *acc_fptr;	/* file of subset accessions */
  char af_name[MAX_FN];
  FILE *libi=NULL;
  FILE *libf;
  int use_stdin;
  struct lmf_str *m_fptr=NULL;
  int acc_off=0;
  char fmt_term;
  char acc_script[MAX_LSTR];
  struct lib_struct *next_lib_p, *this_lib_p, *tmp_lib_p;

  om_fptr = lib_p->m_file_p;

  if (om_fptr != NULL && om_fptr->mm_flg) {
    om_fptr->lpos = 0;
    return om_fptr;
  }

  wcnt = 0;	/* number of times to ask for file name */

  /* check for library type */
  lib_type=0;
  if ((bp=strchr(lib_p->file_name,' '))!=NULL 
      || (bp=strchr(lib_p->file_name,'^'))!=NULL) {
    if (isdigit((int)(bp+1)[0])) {  /* check for number for lib_type */
	*bp='\0';
	sscanf(bp+1,"%d",&lib_type);
	if (lib_type<0 || lib_type >= LASTLIB) {
	  fprintf(stderr,"*** Warning [%s:%d] - invalid library type: %d (>%d)- resetting\n%s\n",
		  __FILE__, __LINE__, lib_type,LASTLIB,lib_p->file_name);
	  lib_type=0;
	}
    }  /* don't change lib_type if its not a number */
  }
  else if (lib_p->file_name[0] =='!') {  /* check for script */
    lib_type = lib_p->lib_type = ACC_SCRIPT;
  }

  /* check for stdin indicator '-' or '@'  (or ACC_SCRIPT) */
  if (lib_p->file_name[0] == '-' || lib_p->file_name[0] == '@'
      || lib_type == ACC_SCRIPT) {
    use_stdin = 1;
  }
  else use_stdin=0;

  if (use_stdin && !(lib_type ==0 || lib_type==ACC_SCRIPT)) {
    fprintf(stderr,"*** Warning [%s:%d] -  @/- STDIN libraries must be in FASTA format\n",__FILE__, __LINE__);
    return NULL;
  }

  opt_text[0]='\0';
  if (lib_type != ACC_SCRIPT) {
  /* check to see if there is a file option ":1-100" */
#ifndef WIN32
    if ((bp=strchr(lib_p->file_name,':'))!=NULL && *(bp+1)!='\0') {
#else
    if ((bp=strchr(lib_p->file_name+3,':'))!=NULL && *(bp+1)!='\0') {
#endif
      SAFE_STRNCPY(opt_text,bp+1,sizeof(opt_text));
      *bp = '\0';
    }
  }

  /* check to see if file can be open()ed? */
 l1:
  opnflg = 0;
  if (lib_type<=LASTTXT) {  /* FASTA, GCG, GENBANK, etc including ACC_SCRIPT */
    if (!use_stdin) {       /* read a file, reading an ACC_SCRIPT is similar to STDIN  */
      opnflg=((libf=fopen(lib_p->file_name,RBSTR))!=NULL);
    }
    else if (lib_type==ACC_SCRIPT) {  /* open a pipe to get the results of the script */
      bp = lib_p->file_name;
      if (lib_p->file_name[0] == '!') {	bp += 1;}
      SAFE_STRNCPY(acc_script, bp, sizeof(acc_script));

      /* convert '+' in annot_script to ' ' */
      bp = strchr(acc_script,'+');
      for ( ; bp; bp=strchr(bp+1,'+')) {
	*bp=' ';
      }
#ifndef WIN32
      libf=popen(acc_script,"r");
#else
      /* windows does not have popen(), but does have _popen() -- not tested */
      libf=_popen(acc_script,"r");
#endif      
      opnflg=1;

    }
    else {  /* just read STDIN */
      libf=stdin;
      lib_p->file_name = alloc_file_name("STDIN");
      opnflg=1;
    }
  } 
  else if (lib_type==ACC_LIST) {
    /* if we have already processed the acc_list file, 
       open the file, modify the acc_list stuff, and return it
    */
    if (lib_p->acc_file_p != NULL) {
      if ((acc_fptr = open_lib(lib_p, ldnaseq, sascii, outtty))==NULL) {
	fprintf(stderr, "*** Warning [%s:%d] - Cannot open %s library for ACC_LIST\n",__FILE__, __LINE__,lib_p->file_name);
	return NULL;
      }
      else {
	/* note that sel_acc_arr[0] must be NULL */
	acc_fptr->sel_acc_p = lib_p->acc_file_p->sel_acc_p;
	acc_fptr->acc_off = lib_p->acc_file_p->acc_off;
	return acc_fptr;
      }
    }

    /* open the file, read the first line, do an openlib on the first line */
    if (!use_stdin) {
      opnflg=((libf=fopen(lib_p->file_name,RBSTR))!=NULL);
    }
    else {
      libf=stdin;
      lib_p->file_name = alloc_file_name("STDIN");
      opnflg=1;
    }

    if (!opnflg) {
      fprintf(stderr, "*** Warning [%s:%d] - Cannot open %s library\n",__FILE__, __LINE__,lib_p->file_name);
      return NULL;
    }
    else {
      /* read in the file line */
      if (fgets(f_line, sizeof(f_line), libf)==NULL) {
	fprintf(stderr, "*** Warning [%s:%d] - Cannot read ACC_LIST file line\n",__FILE__, __LINE__);
	return NULL;
      }
      /* else parse the file line */
      if (f_line[0] != '<') {
	fprintf(stderr, "*** Warning [%s:%d] - missing < - %s\n",__FILE__, __LINE__,f_line); return NULL;
      }
      if ((bp=strchr(f_line+1,'\r'))!=NULL) {*bp = '\0';}
      if ((bp=strchr(f_line+1,'\n'))!=NULL) {*bp = '\0';}

      /* check for accession format */
      if ((bp=strchr(f_line+1,':'))!=NULL) {
	*bp = '\0';
	/* access string should be %d %d%c - acc_ltype, acc_off, fmt_term */
	sscanf(bp+1,"%d %d%c",&acc_ltype, &acc_off, &fmt_term);
	/* blank terminator is default */
	if (acc_off == 0) acc_off = 1;	/* always skip the '>' */
	if (fmt_term < ' ' || fmt_term > '~') fmt_term = ' ';
	if (acc_ltype > MAX_ACC_TYPE) {acc_ltype = MAX_ACC_TYPE;}
      }

      this_lib_p = get_lnames(f_line+1, NULL);

      /* this_lib_p now has the list of files specified by the
	 <@?acc_list_file. If there is only one file, then this
	 information, and the associated m_file_p, should be put into
	 lib_p.  If there is more than one file, the the first should
	 be put in lib_p, and the subsequent lib_struct's should be
	 linked from lib_p->next, and the end of list needs to point
	 to lib_p->next.
      */

      /* done whether list or not */
      lib_p->file_name = this_lib_p->file_name;
      /* deal with other list issues after we have a acc_fptr */

      /* check that we can open the library file */
      if ((acc_fptr = open_lib(this_lib_p, ldnaseq, sascii, outtty))==NULL) {
	fprintf(stderr, "*** Warning [%s:%d] - Cannot open %s library for ACC_LIST\n",__FILE__, __LINE__,f_line+1);
	free(this_lib_p);
	return NULL;
      }
      else {
	/* set up the auxiliary information for the current open file */
	this_lib_p->m_file_p = acc_fptr;
	/* note that sel_acc_arr[0] must be NULL */
	acc_fptr->sel_acc_p = sel_acc_arr[acc_ltype];
	acc_fptr->acc_off = acc_off;
	/* read in the data */
	acc_fptr->sel_local = sel_acc_init[acc_ltype](libf, &acc_fptr->acc_off, fmt_term);

	/* now handle the rest of the this_lib_p list */

	tmp_lib_p = this_lib_p->next;	/* skip over the first entry in this_lib_p */
	/* fill in the information up to the next to last entry in the chain */
	while (tmp_lib_p && tmp_lib_p->next) {
	  tmp_lib_p->acc_file_p = acc_fptr;
	  tmp_lib_p = tmp_lib_p->next;
	}
	if (tmp_lib_p) {
	  tmp_lib_p->acc_file_p = acc_fptr;	/* fill in last entry */
	  tmp_lib_p->next = lib_p->next;	/* continue the chain */
	}
	lib_p->next = this_lib_p->next;		/* insert the chain */
	return acc_fptr;
      }
    } /* done reading ACC_LIST file */
  } /* end of if (lib_type==ACC_LIST */
#ifdef NCBIBL13
  else if (lib_type==NCBIBL13) opnflg=(ncbl_openlib(lib_p->file_name,ldnaseq)!= -1);
#endif
#ifdef NCBIBL20
  else if (lib_type==NCBIBL20) {
    opnflg=((m_fptr=ncbl2_openlib(lib_p,ldnaseq))!=NULL);
  }
#endif

#ifdef MYSQL_DB
  /* a mySQL filename contains mySQL commands, not sequences */
  else if (lib_type==MYSQL_LIB) {
    opnflg=((m_fptr=mysql_openlib(lib_p->file_name,ldnaseq,sascii))!=NULL);
    m_fptr->get_mmap_chain = NULL;
  }
#endif
#ifdef PGSQL_DB
  /* a mySQL filename contains mySQL commands, not sequences */
  else if (lib_type==PGSQL_LIB) {
    opnflg=((m_fptr=pgsql_openlib(lib_p->file_name,ldnaseq,sascii))!=NULL);
    m_fptr->get_mmap_chain = NULL;
  }
#endif

  if (!opnflg) {	/* here if open failed */
    if (outtty) {
      fprintf(stderr,"\n cannot open %s library\n",lib_p->file_name);
      fprintf(stderr," enter new file name or <RET> to quit ");
      fflush(stderr);
      if (fgets(f_line,sizeof(f_line),stdin)==NULL) return NULL;
      if ((bp=strchr(f_line,'\n'))!=0) *bp='\0';
      if (strlen(f_line)==0) return NULL;
      if (++wcnt > 10) return NULL;
      lib_p->file_name = alloc_file_name(f_line);
      goto l1;
    }
    else return NULL;
  }	/* !openflg */

  if (lib_type <= LASTTXT) {
    /* modify to re-use the om_fptr if it exists */
    if (om_fptr != NULL) {
      m_fptr = om_fptr;
    }
    else {
      if ((m_fptr = calloc(1,sizeof(struct lmf_str)))==NULL) {
	fprintf(stderr,"\n*** Warning [%s:%d] - cannot allocate lmf_str (%ld) for %s\n",
		__FILE__, __LINE__,sizeof(struct lmf_str),lib_p->file_name);
	return NULL;
      }
      if ((m_fptr->lline = calloc(MAX_STR,sizeof(char)))==NULL) {
	fprintf(stderr,"\n *** Warning [%s:%d] - cannot allocate lline (%d) for %s\n",
		__FILE__, __LINE__,MAX_STR,lib_p->file_name);
	return NULL;
      }
    }

    m_fptr->lb_name = lib_p->file_name;
    SAFE_STRNCPY(m_fptr->opt_text,opt_text,MAX_FN);
    m_fptr->opt_text[MAX_FN-1]='\0';
    m_fptr->sascii = sascii;
    m_fptr->get_mmap_chain = NULL;

    m_fptr->libf = libf;
    m_fptr->lb_type = lib_type;
    m_fptr->acc_off = 1;	/* default for FASTA format */
    m_fptr->getlib = getliba[lib_type];
    m_fptr->ranlib = ranliba[lib_type];
    m_fptr->sel_acc_p = NULL;
    m_fptr->mm_flg = 0;
    m_fptr->tot_len = 0;
    m_fptr->max_len = 0;
    m_fptr->lib_aa = (ldnaseq==SEQT_PROT);
  }
  last_m_fptr = m_fptr;

#ifdef USE_MMAP
  /* check for possible mmap()ed files */
  if (!use_stdin && (lib_type <= LASTTXT) && can_mmap(lib_type)) {
    /* this is a file we can mmap() */
    /* look for .xin file */
    newname(iname,lib_p->file_name,"xin",sizeof(iname));
    if ((libi=fopen(iname,"r"))!=NULL) { /* have a *.xin file, use mmap */
      if (load_mmap(libi,lib_p->file_name,lib_type,ldnaseq,m_fptr)!=NULL) {
	fclose(libi);	/* close index file */
	return m_fptr;
      }
    fclose(libi);	/* memory mapping failed, but still must close file */
    }
  }
#endif

  if (lib_type <= LASTTXT) {
    m_fptr->lpos = 0;
    if (fgets(m_fptr->lline,MAX_STR,libf)==NULL) return NULL;
  }
  return m_fptr;
}

int
closelib(struct lmf_str *m_fptr,int force) {

  if (m_fptr == NULL) return 0;

#ifdef USE_MMAP
  if (!force && m_fptr->mm_flg) {
/* don't close memory mapped files
*/
    m_fptr->lpos = 0;
    return 0;
  }
#endif

  if (m_fptr->libf!=NULL && m_fptr->libf != stdin) {
    fclose(m_fptr->libf);
    m_fptr->libf = NULL;
    m_fptr->mm_flg = 0;
  }

  /* keep m_fptr->lline around for re-use -- alternatively, always
     allocate */
  /*
  if (m_fptr->lline != NULL) {
    free(m_fptr->lline);
    m_fptr->lline = NULL;
  }
  */

#ifdef NCBIBL13
  if (m_fptr->lb_type == NCBIBL13) ncbl_closelib(m_fptr);
#endif
#ifdef NCBIBL20
  if (m_fptr->lb_type == NCBIBL20) ncbl2_closelib(m_fptr);
#endif
#ifdef MYSQL_DB
  if (m_fptr->lb_type == MYSQL_LIB) mysql_closelib(m_fptr);
#endif

  if (m_fptr == last_m_fptr) {
    last_m_fptr = NULL;
  }

  return 1;
}

struct lmf_str *
  re_openlib(struct lmf_str *om_fptr, struct lib_struct *lib_p, int outtty)
{
  int opnflg;

  /* if its already open, return it */
  if (om_fptr == last_m_fptr) {
    return om_fptr;
  }
  else {
    if (om_fptr->mm_flg) {
      last_m_fptr = om_fptr;
      om_fptr->lpos = 0;
      return om_fptr;
    }
#ifdef MYSQL_DB
    /* if this is a mysql database - use it and return */
    else if (om_fptr->lb_type == MYSQL_LIB) {
      return om_fptr;
    }
#endif
    else {
      closelib(last_m_fptr,1);
    }
  }

  last_m_fptr = om_fptr;


  /* data is available, but file is closed or not memory mapped, open it */
  /* no longer check to memory map - because we could not do it before */

  opnflg = 1;
  if (om_fptr->lb_type<=LASTTXT && om_fptr->libf==NULL)
    opnflg=((om_fptr->libf=fopen(om_fptr->lb_name,RBSTR))!=NULL);
#ifdef NCBIBL20
  else if (om_fptr->lb_type==NCBIBL20) {
    opnflg=((om_fptr=ncbl2_reopen(om_fptr))!=NULL);
  }
#endif
#ifdef MYSQL_DB
  /* a mySQL filename contains mySQL commands, not sequences */
  else if (om_fptr->lb_type==MYSQL_LIB) 
    opnflg=(mysql_reopen(om_fptr)!=NULL);
#endif

  if (!opnflg) {
    fprintf(stderr,"*** Warning [%s:%d] - could not re_open %s\n",__FILE__, __LINE__,om_fptr->lb_name);
    return NULL;
  }

  /* use the old buffer for the opened text file */
  return om_fptr;
}

void sf_sort(int *, int);

int
agetlib(unsigned char *seq, int maxs,
	char *libstr, int n_libstr,
	fseek_t *libpos,
	int *lcont,
	struct lmf_str *lm_fd,
	long *l_off)
{
  int i;
  register unsigned char *cp, *seqp, *seqb;
  register int *ap;
  unsigned char *seqm, *seqm1;
  /* int ic, l_start, l_stop, l_limit, rn; */
  char *bp, *bp1, *bpa, *tp;
  int sel_status;

  seqp = seqb = seq;
  seqm = &seq[maxs-9];
  seqm1 = seqm-1;

  ap = lm_fd->sascii;

  *l_off = 1;
  if (*lcont==0) {

    while (lm_fd->lline[0]=='#') {
      if (lm_fd->libf != stdin) lm_fd->lpos = FTELL(lm_fd->libf);
      if (fgets(lm_fd->lline,MAX_STR,lm_fd->libf)==NULL) return (-1);
    }

  start_seq:
    while (lm_fd->lline[0]!='>' && lm_fd->lline[0]!=';' ) {
      if (lm_fd->libf != stdin) lm_fd->lpos = FTELL(lm_fd->libf);
      if (fgets(lm_fd->lline,MAX_STR,lm_fd->libf)==NULL) return (-1);
    }

    /* get l_off coordinate from @C:123 */
    if ((bp=strchr(lm_fd->lline,'@'))!=NULL && !strncmp(bp+1,"C:",2)) {
      sscanf(bp+3,"%ld",l_off);
    }

    SAFE_STRNCPY(libstr,lm_fd->lline+lm_fd->acc_off,n_libstr);

    if ((lm_fd->sel_acc_p != NULL) &&
	(sel_status = (lm_fd->sel_acc_p)(libstr, 0, lm_fd->sel_local)) <= 0) {
      if (sel_status < 0) return (-1);
      while (strchr((char *)lm_fd->lline,'\n')==NULL) {
	if (fgets(lm_fd->lline,MAX_STR,lm_fd->libf)==NULL) return (-1);
      }
      if (fgets(lm_fd->lline,MAX_STR,lm_fd->libf)==NULL) return (-1);
      goto start_seq;
    }

    if ((bp=strchr(libstr,'\r'))!=NULL) *bp='\0';
    if ((bp=strchr(libstr,'\n'))!=NULL) *bp='\0';

    if (n_libstr > MAX_UID) {
      tp = libstr;
      while (*tp++) if (*tp == '\001' || *tp== '\t') *tp = ' ';
    }

    *libpos = lm_fd->lpos;

    /* make certain we have the end of the line */
    while (strchr((char *)lm_fd->lline,'\n')==NULL) {
      if (strlen(lm_fd->lline)<MAX_STR/2) 
	fgets(&lm_fd->lline[strlen(lm_fd->lline)],MAX_STR/2,lm_fd->libf);
      else 
	fgets(&lm_fd->lline[MAX_STR/2],MAX_STR/2,lm_fd->libf);
    }
    lm_fd->lline[MAX_STR-1]='\0';
  }

  lm_fd->lline[0]='\0';
  while (seqb<seqm1 && fgets((char *)seqb,(size_t)(seqm-seqb),lm_fd->libf)!=NULL) {
    if (*seqb=='>') goto new;
    if (*seqb=='#' || *seqb==';') {
      if (strchr((char *)seqb,'\n')==NULL) goto cont;
      continue;
    }

    /* removed - used for @P:1-n 
       if (l_limit) {
       for (cp=seqp; seqp<seqm1 && rn < l_stop && (ic=ap[*cp++])<EL; )
       if (ic < NA && ++rn > l_start) *seqp++ = (unsigned char)ic;
       if (rn > l_stop) goto finish;
       }
       else {
    */
    seqp = seqb;
    for (cp=seqp; seqp<seqm1; ) {
      if ((*seqp++=ap[*cp++])<NA &&
	  (*seqp++=ap[*cp++])<NA &&
	  (*seqp++=ap[*cp++])<NA &&
	  (*seqp++=ap[*cp++])<NA &&
	  (*seqp++=ap[*cp++])<NA &&
	  (*seqp++=ap[*cp++])<NA &&
	  (*seqp++=ap[*cp++])<NA &&
	  (*seqp++=ap[*cp++])<NA &&
	  (*seqp++=ap[*cp++])<NA) continue;
      if (*(--seqp)>NA) break;
    }
    seqb = seqp;
    if (*seqp==ES) goto done;
    if (lm_fd->libf != stdin) lm_fd->lpos = FTELL(lm_fd->libf);
  }
  goto done;
 new:
  SAFE_STRNCPY(lm_fd->lline,(char *)seqp,MAX_STR);
  /* be certain to get complete line, if possible */
  if (strchr(lm_fd->lline,'\n')==NULL)
    fgets(&lm_fd->lline[strlen(lm_fd->lline)],MAX_STR-strlen(lm_fd->lline),lm_fd->libf);
  lm_fd->lline[MAX_STR-1]='\0';
  if (strchr(lm_fd->lline,'\n')==NULL && strchr((char *)seqp,'\n')!=NULL)
    lm_fd->lline[strlen(lm_fd->lline)-1]='\n';
  goto done;

  /* removed - used for @P:1-n
finish: 
   while (lm_fd->lline[0]!='>' && 
	  fgets(lm_fd->lline,MAX_STR,lm_fd->libf)!=NULL) {
     if (lm_fd->libf != stdin) lm_fd->lpos = FTELL(lm_fd->libf);
   }
   goto done;
*/
 cont:
  fgets(lm_fd->lline,MAX_STR,lm_fd->libf);
  seqm1 = seqp;
 done:
  if (seqp>=seqm1) (*lcont)++;
  else {
    *lcont=0;
  }

  *seqp = EOSEQ;
  /*  if ((int)(seqp-seq)==0) return 1; */
  return (int)(seqp-seq);
}

void
aranlib(char *str, int cnt, fseek_t seek, char *libstr, struct lmf_str *lm_fd)
{
  char *bp;

  if (lm_fd->libf != stdin) {
    FSEEK(lm_fd->libf, seek, 0);
    fgets(lm_fd->lline,MAX_STR,lm_fd->libf);

    if (lm_fd->lline[0]=='>' || lm_fd->lline[0]==';') {
      SAFE_STRNCPY(str,lm_fd->lline+lm_fd->acc_off,cnt);
      if ((bp = strchr(str,'\r'))!=NULL) *bp='\0';
      if ((bp = strchr(str,'\n'))!=NULL) *bp='\0';
      /*
	if ((bp = strchr(str,SFCHAR))!=NULL) *bp='\0';
	else if ((bp = strchr(str,'\001'))!=NULL) *bp='\0';
	else if ((bp = strchr(str,'\n'))!=NULL) *bp='\0';
	else str[cnt-1]='\0';
      */
      bp = str;
      while (*bp++) if (*bp=='\001' || *bp=='\t') *bp=' ';
    }
    else {
      str[0]='\0';
    }
  }
  else str[0]='\0';
}

int
qgetlib(unsigned char *seq, int maxs,
	char *libstr, int n_libstr,
	fseek_t *libpos,
	int *lcont,
	struct lmf_str *lm_fd,
	long *l_off)
{
  int i;
  register unsigned char *cp, *seqp, *seqb;
  register int *ap;
  unsigned char *seqm, *seqm1;
  /* int ic, l_start, l_stop, l_limit, rn; */
  char *bp, *bp1, *bpa, *tp;
  int sel_status;

  seqp = seqb = seq;
  seqm = &seq[maxs-9];
  seqm1 = seqm-1;

  ap = lm_fd->sascii;

  *l_off = 1;
  if (*lcont==0) {

    while (lm_fd->lline[0]!='@') {
      if (lm_fd->libf != stdin) lm_fd->lpos = FTELL(lm_fd->libf);
      if (fgets(lm_fd->lline,MAX_STR,lm_fd->libf)==NULL) return (-1);
    }

    SAFE_STRNCPY(libstr,lm_fd->lline+lm_fd->acc_off,n_libstr);

    if ((bp=strchr(libstr,'\r'))!=NULL) *bp='\0';
    if ((bp=strchr(libstr,'\n'))!=NULL) *bp='\0';

    *libpos = lm_fd->lpos;

    /* make certain we have the end of the line */
    while (strchr((char *)lm_fd->lline,'\n')==NULL) {
      if (strlen(lm_fd->lline)<MAX_STR/2) 
	fgets(&lm_fd->lline[strlen(lm_fd->lline)],MAX_STR/2,lm_fd->libf);
      else 
	fgets(&lm_fd->lline[MAX_STR/2],MAX_STR/2,lm_fd->libf);
    }
    lm_fd->lline[MAX_STR-1]='\0';
  }

  lm_fd->lline[0]='\0';
  while (seqb<seqm1 && fgets((char *)seqb,(size_t)(seqm-seqb),lm_fd->libf)!=NULL) {
    if (*seqb=='+') goto new;

    seqp = seqb;
    for (cp=seqp; seqp<seqm1; ) {
      if ((*seqp++=ap[*cp++])<NA &&
	  (*seqp++=ap[*cp++])<NA &&
	  (*seqp++=ap[*cp++])<NA &&
	  (*seqp++=ap[*cp++])<NA &&
	  (*seqp++=ap[*cp++])<NA &&
	  (*seqp++=ap[*cp++])<NA &&
	  (*seqp++=ap[*cp++])<NA &&
	  (*seqp++=ap[*cp++])<NA &&
	  (*seqp++=ap[*cp++])<NA) continue;
      if (*(--seqp)>NA) break;
    }
    seqb = seqp;
    if (*seqp==ES) goto done;
    if (lm_fd->libf != stdin) lm_fd->lpos = FTELL(lm_fd->libf);
  }
  goto done;
 new:
  SAFE_STRNCPY(lm_fd->lline,(char *)seqp,MAX_STR);
  /* be certain to get complete line, if possible */
  if (strchr(lm_fd->lline,'\n')==NULL)
    fgets(&lm_fd->lline[strlen(lm_fd->lline)],MAX_STR-strlen(lm_fd->lline),lm_fd->libf);
  lm_fd->lline[MAX_STR-1]='\0';
  if (strchr(lm_fd->lline,'\n')==NULL && strchr((char *)seqp,'\n')!=NULL)
    lm_fd->lline[strlen(lm_fd->lline)-1]='\n';

 done:
  if (seqp>=seqm1) (*lcont)++;
  else {
    *lcont=0;
  }

  *seqp = EOSEQ;
  /*  if ((int)(seqp-seq)==0) return 1; */
  return (int)(seqp-seq);
}

void
qranlib(char *str, int cnt, fseek_t seek, char *libstr, struct lmf_str *lm_fd)
{
  char *bp;

  if (lm_fd->libf != stdin) {
    FSEEK(lm_fd->libf, seek, 0);
    fgets(lm_fd->lline,MAX_STR,lm_fd->libf);

    if (lm_fd->lline[0]=='@') {
      SAFE_STRNCPY(str,lm_fd->lline+lm_fd->acc_off,cnt);
      if ((bp = strchr(str,'\r'))!=NULL) *bp='\0';
      if ((bp = strchr(str,'\n'))!=NULL) *bp='\0';
      bp = str;
    }
    else {
      str[0]='\0';
    }
  }
  else str[0]='\0';
}

void lget_ann(struct lmf_str *, char *, int);

int
lgetlib(unsigned char *seq, int maxs,
	char *libstr,
	int n_libstr,
	fseek_t *libpos,
	int *lcont,
	struct lmf_str *lm_fd,
	long *l_off)
{
  register unsigned char *cp, *seqp;
  register int *ap;
  unsigned char *seqm, *seqm1;
  char *bp, *bp_gid;

  *l_off = 1;

  seqp = seq;
  seqm = &seq[maxs-11];
  seqm1 = seqm-1;

  ap = lm_fd->sascii;

  if (*lcont==0) {
    while (lm_fd->lline[0]!='L' || lm_fd->lline[1]!='O' || 
	   strncmp(lm_fd->lline,"LOCUS",5)) { /* find LOCUS */
      lm_fd->lpos = FTELL(lm_fd->libf);
      if (fgets(lm_fd->lline,MAX_STR,lm_fd->libf)==NULL) return (-1);
      if (lm_fd->lfflag) getc(lm_fd->libf);
    }
    *libpos= lm_fd->lpos;

    if (n_libstr <= 21) {
      SAFE_STRNCPY(libstr,&lm_fd->lline[12],12);
    }
    else {
      lget_ann(lm_fd,libstr,n_libstr);
      fgets(lm_fd->lline,MAX_STR,lm_fd->libf);
    }

    while (lm_fd->lline[0]!='O' || lm_fd->lline[1]!='R' ||
	   strncmp(lm_fd->lline,"ORIGIN",6)) { /* find ORIGIN */
      if (fgets(lm_fd->lline,MAX_STR,lm_fd->libf)==NULL) return (-1);
      if (lm_fd->lfflag) getc(lm_fd->libf);
    }
  }
  else {
    for (cp= lm_fd->cpsave; seqp<seqm1; ) {
      if ((*seqp++=ap[*cp++])<NA) continue;
      if (*(--seqp)>NA) break;
    }
  }

  lm_fd->lline[0]='\0';
  while (seqp<seqm1 && fgets(lm_fd->lline,MAX_STR,lm_fd->libf)!=NULL) {
    if (lm_fd->lfflag) getc(lm_fd->libf);
    if (lm_fd->lline[0]=='/') goto new;
    for (cp= (unsigned char *)&lm_fd->lline[10]; seqp<seqm1; ) {
      if ((*seqp++=ap[*cp++])<NA &&
	  (*seqp++=ap[*cp++])<NA &&
	  (*seqp++=ap[*cp++])<NA &&
	  (*seqp++=ap[*cp++])<NA &&
	  (*seqp++=ap[*cp++])<NA &&
	  (*seqp++=ap[*cp++])<NA &&
	  (*seqp++=ap[*cp++])<NA &&
	  (*seqp++=ap[*cp++])<NA &&
	  (*seqp++=ap[*cp++])<NA &&
	  (*seqp++=ap[*cp++])<NA &&
	  (*seqp++=ap[*cp++])<NA) continue;
      if (*(--seqp)>NA) break;
    }
  }
  goto done;
new:
  lm_fd->lpos = FTELL(lm_fd->libf);
  fgets(lm_fd->lline,MAX_STR,lm_fd->libf);
  if (lm_fd->lfflag) getc(lm_fd->libf);

done:
  if (seqp>=seqm1) {
    lm_fd->cpsave = cp;
    (*lcont)++;
  }
  else *lcont=0;

  *seqp = EOSEQ;
  /*  if ((int)(seqp-seq)==0) return 1; */
  return (int)(seqp-seq);
}

void
lget_ann(struct lmf_str *lm_fd, char *libstr, int n_libstr) {
  char *bp, *bp_gid, locus[120], desc[120], acc[120], ver[120];

  /* copy in locus from lm_fd->lline */
  SAFE_STRNCPY(locus,&lm_fd->lline[12],sizeof(locus));
  if ((bp=strchr(locus,' '))!=NULL) *(bp+1) = '\0';

  /* get description */
  fgets(desc,sizeof(desc),lm_fd->libf);
  while (desc[0]!='D' || desc[1]!='E' || strncmp(desc,"DEFINITION",10))
    fgets(desc,sizeof(desc),lm_fd->libf);
  if ((bp = strchr(&desc[12],'\n'))!=NULL) *bp='\0';

  /* get accession */
  fgets(acc,sizeof(acc),lm_fd->libf);
  while (acc[0]!='A' || acc[1]!='C' || strncmp(acc,"ACCESSION",9)) {
    fgets(acc,sizeof(acc),lm_fd->libf);
    if (acc[0]=='O' && acc[1]=='R' && strncmp(acc,"ORIGIN",6)==0)
      break;
  }
  if ((bp = strchr(&acc[12],'\n'))!=NULL) *bp='\0';
  if ((bp = strchr(&acc[12],' '))!=NULL) *bp='\0';

  /* get version */
  fgets(ver,sizeof(ver),lm_fd->libf);
  while (ver[0]!='V' || ver[1]!='E' || strncmp(ver,"VERSION",7)) {
    fgets(ver,sizeof(ver),lm_fd->libf);
    if (ver[0]=='O' && ver[1]=='R' && strncmp(ver,"ORIGIN",6)==0)
      break;
  }
  if ((bp = strchr(&ver[12],'\n'))!=NULL) *bp='\0';

      /* extract gi:123456 from version line */
  bp_gid = strchr(&ver[12],':');
  if (bp_gid != NULL) {
    if ((bp=strchr(bp_gid+1,' '))!=NULL) *bp='\0';
    bp_gid++;
  }
  if ((bp = strchr(&ver[12],' '))!=NULL) *bp='\0';

      /* build up FASTA header line */
  if (bp_gid != NULL) {
    SAFE_STRNCPY(libstr,"gi|",n_libstr);
    SAFE_STRNCAT(libstr,bp_gid,n_libstr);
    SAFE_STRNCAT(libstr,"|gb|",n_libstr);
  }
  else {libstr[0]='\0';}

  /* if we have a version number, use it, otherwise accession, 
	 otherwise locus/description */

  if (ver[0]=='V') {
    SAFE_STRNCAT(libstr,&ver[12],n_libstr);
    SAFE_STRNCAT(libstr,"|",n_libstr);
  }
  else if (acc[0]=='A') {
    SAFE_STRNCAT(libstr,&acc[12],n_libstr);
    SAFE_STRNCAT(libstr," ",n_libstr);
  }

  SAFE_STRNCAT(libstr,locus,n_libstr);
  SAFE_STRNCAT(libstr,&desc[11],n_libstr);
  libstr[n_libstr-1]='\0';
}

/* this code seeks to provide both the various accession numbers
   necessary to identify the sequence, and also some description.

   Unfortunately, the various contributors to Genbank use three
   slightly different formats for including the accession number.

(1)LOCUS       HSJ214M20  107422 bp    DNA             HTG       16-JUN-2000
   DEFINITION  Homo sapiens chromosome 6 clone RP1-214M20 map p12.1-12.3, ***
               SEQUENCING IN PROGRESS ***, in unordered pieces.
   ACCESSION   AL121969

(2)LOCUS       AL359201   117444 bp    DNA             HTG       15-JUN-2000
   DEFINITION  Homo sapiens chromosome 1 clone RP4-671C13 map p13.2-21.1, ***
               SEQUENCING IN PROGRESS ***, in unordered pieces.
   ACCESSION   AL359201

(3)LOCUS       BB067000      280 bp    mRNA            EST       19-JUN-2000
   DEFINITION  BB067000 RIKEN full-length enriched, 15 days embryo male testis Mus
               musculus cDNA clone 8030456L01 3', mRNA sequence.
   ACCESSION   BB067000

This makes it more difficult to both provide the accession number in a
standard location and to conserve definition space
*/

void
lranlib(char *str,
	int cnt,
	fseek_t seek,
	char *libstr,
	struct lmf_str *lm_fd)
{
  char *bp, acc[MAX_STR], desc[MAX_STR];

  FSEEK(lm_fd->libf, seek, 0);
  fgets(lm_fd->lline,MAX_STR,lm_fd->libf);
  if (lm_fd->lfflag) getc(lm_fd->libf);

  lget_ann(lm_fd, str, cnt);
  str[cnt-1]='\0';

  FSEEK(lm_fd->libf,seek,0);
  fgets(lm_fd->lline,MAX_STR,lm_fd->libf);
  if (lm_fd->lfflag) getc(lm_fd->libf);
}

int
pgetlib(unsigned char *seq, int maxs,
	char *libstr,
	int n_libstr,
	fseek_t *libpos,
	int *lcont,
	struct lmf_str *lm_fd,
	long *l_off)
{
  int ic;
  register unsigned char *cp, *seqp;
  register int *ap;
  unsigned char *seqm, *seqm1;

  *l_off = 1;

  seqp = seq;
  seqm = &seq[maxs-11];
  seqm1 = seqm-1;

  ap = lm_fd->sascii;

  if (*lcont==0) {
    while (lm_fd->lline[0]!='E' || lm_fd->lline[1]!='N' || strncmp(lm_fd->lline,"ENTRY",5))
      { /* find ENTRY */
	lm_fd->lpos = FTELL(lm_fd->libf);
	if (fgets(lm_fd->lline,MAX_STR,lm_fd->libf)==NULL) return (-1);
      }
    SAFE_STRNCPY(libstr,&lm_fd->lline[16],8);
    *libpos = lm_fd->lpos;
    while (lm_fd->lline[2]!='Q' || lm_fd->lline[0]!='S' || strncmp(lm_fd->lline,"SEQUENCE",8))
      { /* find SEQUENCE */
	if (fgets(lm_fd->lline,MAX_STR,lm_fd->libf)==NULL) return (-1);
      }
    fgets(lm_fd->lline,MAX_STR,lm_fd->libf); /* get the extra line */
  }
  else {
    for (cp= lm_fd->cpsave; seqp<seqm1; ) {
      if ((*seqp++=ap[*cp++])<NA) continue;
      if (*(--seqp)>NA) break;
    }
    if (*seqp==ES) goto done;
  }

  lm_fd->lline[0]='\0';
  while (seqp<seqm1 && fgets(lm_fd->lline,MAX_STR,lm_fd->libf)!=NULL) {
    if (lm_fd->lline[0]=='/') goto new;
    for (cp= (unsigned char *)&lm_fd->lline[8]; seqp<seqm1; ) {
      if ((*seqp++=ap[*cp++])<NA) continue;
      if (*(--seqp)>NA) break;
    };
    if (*seqp==ES) goto done;
  }
  goto done;
new:
  lm_fd->lpos = FTELL(lm_fd->libf);
  fgets(lm_fd->lline,MAX_STR,lm_fd->libf);

done:
  if (seqp>=seqm1) {
    lm_fd->cpsave = cp;
    (*lcont)++;
  }
  else *lcont=0;

  *seqp = EOSEQ;
  /*  if ((int)(seqp-seq)==0) return 1; */
  return (int)(seqp-seq);
}

void
pranlib(char *str,
	int cnt,
	fseek_t seek,
	char *libstr,
	struct lmf_str *lm_fd)
{
  char *bp;

  FSEEK(lm_fd->libf, seek, 0);
  fgets(lm_fd->lline,MAX_STR,lm_fd->libf);

  SAFE_STRNCPY(str,&lm_fd->lline[16],8);
  fgets(lm_fd->lline,MAX_STR,lm_fd->libf);
  while (lm_fd->lline[0]!='T' || lm_fd->lline[1]!='I' || strncmp(lm_fd->lline,"TITLE",5))
    fgets(lm_fd->lline,MAX_STR,lm_fd->libf);

  SAFE_STRNCAT(str,&lm_fd->lline[16],cnt);
  if ((bp = strchr(str,'\n'))!=NULL) *bp='\0';

  FSEEK(lm_fd->libf,seek,0);
  fgets(lm_fd->lline,MAX_STR,lm_fd->libf);
}

int
egetlib(unsigned char *seq, int maxs,
	char *libstr,
	int n_libstr,
	fseek_t *libpos,
	int *lcont,
	struct lmf_str *lm_fd,
	long *l_off)
{
  int ll;
  int ic;
  register unsigned char *cp, *seqp;
  register int *ap;
  unsigned char *seqm, *seqm1;
  int sel_status;
  char id[11];  /* Holds Identifier */

  *l_off=1;

  seqp = seq;
  seqm = &seq[maxs-11];
  seqm1 = seqm-1;

  ap = lm_fd->sascii;

  if (*lcont==0) {
  start_seq:
    while (lm_fd->lline[0]!='I' || lm_fd->lline[1]!='D') { /* find ID */
      lm_fd->lpos = FTELL(lm_fd->libf);
      if (fgets(lm_fd->lline,MAX_STR,lm_fd->libf)==NULL) return (-1);
      if (lm_fd->lfflag) getc(lm_fd->libf);
    }
    sscanf(&lm_fd->lline[5],"%s",id);
    sprintf(libstr,"%-12.12s",id);
    libstr[12]='\0';

    if ((lm_fd->sel_acc_p != NULL) &&
	(sel_status = (lm_fd->sel_acc_p)(libstr, 0, lm_fd->sel_local)) <= 0) {
      if (sel_status < 0) return (-1);
      if (fgets(lm_fd->lline,MAX_STR,lm_fd->libf)==NULL) return (-1);
      goto start_seq;
    }

    *libpos = lm_fd->lpos;
    while (lm_fd->lline[0]!='S' || lm_fd->lline[1]!='Q') { /* find ORIGIN */
      if (fgets(lm_fd->lline,MAX_STR,lm_fd->libf)==NULL) return (-1);
      if (lm_fd->lfflag) getc(lm_fd->libf);
    }
    sscanf(&lm_fd->lline[14],"%ld",&lm_fd->gcg_len);
  }
  else {
    for (cp= lm_fd->cpsave; seqp<seqm1; ) {
      if ((*seqp++=ap[*cp++])<NA) continue;
      if (*(--seqp)>NA) break;
    }
    if (*seqp==ES) goto done;
  }

  lm_fd->lline[0]='\0';
  while (seqp<seqm1 && fgets(lm_fd->lline,MAX_STR,lm_fd->libf)!=NULL) {
    if (lm_fd->lfflag) getc(lm_fd->libf);
    if (lm_fd->lline[0]=='/') goto new;
    lm_fd->lline[70]='\0';
    for (cp= (unsigned char *)&lm_fd->lline[5]; seqp<seqm1; ) {
      if ((*seqp++=ap[*cp++])<NA &&
	  (*seqp++=ap[*cp++])<NA &&
	  (*seqp++=ap[*cp++])<NA &&
	  (*seqp++=ap[*cp++])<NA &&
	  (*seqp++=ap[*cp++])<NA &&
	  (*seqp++=ap[*cp++])<NA &&
	  (*seqp++=ap[*cp++])<NA &&
	  (*seqp++=ap[*cp++])<NA &&
	  (*seqp++=ap[*cp++])<NA &&
	  (*seqp++=ap[*cp++])<NA &&
	  (*seqp++=ap[*cp++])<NA) continue;
      if (*(--seqp)>NA) break;
    }
    if (*seqp==ES) goto done;
  }
  goto done;
new:	lm_fd->lpos = FTELL(lm_fd->libf);
fgets(lm_fd->lline,MAX_STR,lm_fd->libf);
if (lm_fd->lfflag) getc(lm_fd->libf);
goto done;

done:	if (seqp>=seqm1) {
  lm_fd->cpsave = cp;
  (*lcont)++;
  lm_fd->gcg_len -= (long)(seqp-seq);
}
else *lcont=0;

*seqp = EOSEQ;
/* if ((int)(seqp-seq)==0) return 1; */
/*	if (*lcont==0 && (long)(seqp-seq)!=lm_fd->gcg_len)
	printf("%s read %d of %d\n",libstr,(int)(seqp-seq),lm_fd->gcg_len);
	*/
return (int)(seqp-seq);
}

void
eranlib(char *str,
	int cnt,
	fseek_t seek,
	char *libstr,
	struct lmf_str *lm_fd)
{
  char *bp;
  char id[14];  /* Holds Identifier */

  FSEEK(lm_fd->libf, seek, 0);
  fgets(lm_fd->lline,MAX_STR,lm_fd->libf);
  if (lm_fd->lfflag) getc(lm_fd->libf);

  sscanf(&lm_fd->lline[5],"%s",id);
  sprintf(str,"%-12.12s ",id);
  fgets(lm_fd->lline,MAX_STR,lm_fd->libf);
  if (lm_fd->lfflag) getc(lm_fd->libf);

  while (lm_fd->lline[0]!='D' || lm_fd->lline[1]!='E') fgets(lm_fd->lline,MAX_STR,lm_fd->libf);

  SAFE_STRNCAT(str,&lm_fd->lline[5],cnt);
  if ((bp = strchr(str,'\r'))!=NULL) *bp='\0';
  if ((bp = strchr(str,'\n'))!=NULL) *bp='\0';

  FSEEK(lm_fd->libf,seek,0);
  fgets(lm_fd->lline,MAX_STR,lm_fd->libf);
  if (lm_fd->lfflag) getc(lm_fd->libf);
}

int
igetlib(unsigned char *seq, int maxs,
	char *libstr,
	int n_libstr,
	fseek_t *libpos,
	int *lcont,
	struct lmf_str *lm_fd,
	long *l_off)
{
	register unsigned char *cp, *seqp;
	register int *ap;
	unsigned char *seqm, *seqm1;
	char *bp;

	*l_off = 1;

	seqp = seq;
	seqm = &seq[maxs-9];
	seqm1 = seqm-1;

	ap = lm_fd->sascii;

	if (*lcont==0) {
		while (lm_fd->lline[0]==';') {
			lm_fd->lpos = FTELL(lm_fd->libf);
			if (fgets(lm_fd->lline,MAX_STR,lm_fd->libf)==NULL) return (-1);
			}
		*libpos = lm_fd->lpos;
		while (lm_fd->lline[0]==';') {
		  fgets(lm_fd->lline,MAX_STR,lm_fd->libf);
		}
		SAFE_STRNCPY(libstr,lm_fd->lline+1,12);
		if((bp=strchr(libstr,'\n'))!=NULL) *bp='\0';
		}

	lm_fd->lline[0]='\0';
	while (seqp<seqm1 && fgets((char *)seqp,(size_t)(seqm-seqp),lm_fd->libf)!=NULL) {
		if (*seqp=='>') goto new;
		if (*seqp==';') {
			if (strchr((char *)seqp,'\n')==NULL) goto cont;
			continue;
			}
		for (cp=seqp; seqp<seqm1; ) {
			if ((*seqp++=ap[*cp++])<NA &&
			    (*seqp++=ap[*cp++])<NA &&
		            (*seqp++=ap[*cp++])<NA &&
		            (*seqp++=ap[*cp++])<NA &&
			    (*seqp++=ap[*cp++])<NA &&
			    (*seqp++=ap[*cp++])<NA &&
			    (*seqp++=ap[*cp++])<NA &&
			    (*seqp++=ap[*cp++])<NA &&
			    (*seqp++=ap[*cp++])<NA) continue;
			    if (*(--seqp)>NA) break;
			    }
		if (*seqp==ES) goto done;
		lm_fd->lpos = FTELL(lm_fd->libf);
		}
	goto done;
new:	SAFE_STRNCPY(lm_fd->lline,(char *)seqp,MAX_STR);
	if (strchr((char *)seqp,'\n')==NULL)
	    fgets(lm_fd->lline,MAX_STR-strlen(lm_fd->lline),lm_fd->libf);
	goto done;

cont:
	fgets(lm_fd->lline,MAX_STR,lm_fd->libf);
	seqm1 = seqp;

done:	if (seqp>=seqm1) {
		(*lcont)++;
		}
	else {
	*lcont=0;
		}


	*seqp = EOSEQ;
	/*	if ((int)(seqp-seq)==0) return 1; */
	return (int)(seqp-seq);
	}

void
iranlib(char *str,
	int cnt,
	fseek_t seek,
	char *libstr,
	struct lmf_str *lm_fd)
{
	char *bp;
	char tline[MAX_FN];

	FSEEK(lm_fd->libf, seek, 0);
	fgets(lm_fd->lline,MAX_STR,lm_fd->libf);

	if (lm_fd->lline[0]=='>' || lm_fd->lline[0]==';') {
		SAFE_STRNCPY(tline,lm_fd->lline+1,sizeof(tline));
		if ((bp = strchr(tline,'\n'))!=NULL) *bp='\0';
		}
	else {
		tline[0]='\0';
		}

	while (lm_fd->lline[0]==';') {
	  fgets(lm_fd->lline,MAX_STR,lm_fd->libf);
	}
	if ((bp=strchr(lm_fd->lline,'\n'))!=NULL) *bp=0;
	if ((bp=strchr(lm_fd->lline,' '))!=NULL) *bp=0;
	SAFE_STRNCPY(str,lm_fd->lline,cnt);
	SAFE_STRNCAT(str,"  ",cnt);
	SAFE_STRNCAT(str,tline,cnt);
	
	FSEEK(lm_fd->libf,seek,0);
	fgets(lm_fd->lline,MAX_STR,lm_fd->libf);
	}

int
vgetlib(unsigned char *seq, int maxs,
	char *libstr,
	int n_libstr,
	fseek_t *libpos,
	int *lcont,
	struct lmf_str *lm_fd,
	long *l_off)
{
  int i, ich;
  register unsigned char *cp, *seqp;
  register int *ap;
  unsigned char *seqm, *seqm1;
  char *bp, *tp;

  *l_off = 1;

  seqp = seq;
  seqm = &seq[maxs-9];
  seqm1 = seqm-1;

  ap = lm_fd->sascii;

  if (*lcont==0) {
    while (lm_fd->lline[0]!='>' && lm_fd->lline[0]!=';') {
      lm_fd->lpos = FTELL(lm_fd->libf);
      if (fgets(lm_fd->lline,MAX_STR,lm_fd->libf)==NULL) return (-1);
      if (lm_fd->lfflag) getc(lm_fd->libf);
    }

    if ((bp=strchr(lm_fd->lline,'\n'))!=NULL) *bp='\0';
    SAFE_STRNCPY(libstr,&lm_fd->lline[4],12);
    if ((bp=strchr(libstr,' '))!=NULL) *bp='\0';
    if ((bp=strchr(libstr,'\n'))!=NULL) *bp='\0';
    
    fgets(lm_fd->lline,MAX_STR,lm_fd->libf);
    if (lm_fd->lfflag) getc(lm_fd->libf);

    if (n_libstr > 21) {
      strcat(libstr," ");
      SAFE_STRNCAT(libstr,lm_fd->lline,n_libstr);
      if ((bp=strchr(libstr,'\n'))!=NULL) *bp='\0';
      libstr[n_libstr-1]='\0';
    }
    *libpos = lm_fd->lpos;
  }

  lm_fd->lline[0]='\0';
  while (seqp<seqm1 && fgets((char *)seqp,(size_t)(seqm-seqp),lm_fd->libf)!=NULL) {
    if (lm_fd->lfflag && (ich=getc(lm_fd->libf))!=LFCHAR) ungetc(ich,lm_fd->libf);
    if (*seqp=='>') goto new;
    if (*seqp==';') {
      if (strchr((char *)seqp,'\n')==NULL) goto cont;
      continue;
    }
    for (cp=seqp; seqp<seqm1; ) {
      if ((*seqp++=ap[*cp++])<NA &&
	  (*seqp++=ap[*cp++])<NA &&
	  (*seqp++=ap[*cp++])<NA &&
	  (*seqp++=ap[*cp++])<NA &&
	  (*seqp++=ap[*cp++])<NA &&
	  (*seqp++=ap[*cp++])<NA &&
	  (*seqp++=ap[*cp++])<NA &&
	  (*seqp++=ap[*cp++])<NA &&
	  (*seqp++=ap[*cp++])<NA) continue;
      if (*(--seqp)>NA) break;
    }
    if (*seqp==ES) goto done;
    lm_fd->lpos = FTELL(lm_fd->libf);
  }
  goto done;
new:
  SAFE_STRNCPY(lm_fd->lline,(char *)seqp,MAX_STR);
  if (strchr((char *)seqp,'\n')==NULL) {
    fgets(&lm_fd->lline[strlen(lm_fd->lline)],MAX_STR-strlen(lm_fd->lline),lm_fd->libf);
    if (lm_fd->lfflag && (ich=getc(lm_fd->libf))!=LFCHAR) ungetc(ich,lm_fd->libf);
  }
  goto done;

cont:
  fgets(lm_fd->lline,MAX_STR,lm_fd->libf);
  if (lm_fd->lfflag && (ich=getc(lm_fd->libf))!=LFCHAR) ungetc(ich,lm_fd->libf);
  seqm1 = seqp;

done:
  if (seqp>=seqm1) {
    (*lcont)++;
  }
  else {
    *lcont=0;
  }

  *seqp = EOSEQ;
  /*   if ((int)(seqp-seq)==0) return 1;*/
  return (int)(seqp-seq);
}

void
vranlib(char *str,
	int cnt,
	fseek_t seek,
	char *libstr,
	struct lmf_str *lm_fd)
{
  char *bp, *llp;

  FSEEK(lm_fd->libf, seek, 0);
  fgets(lm_fd->lline,MAX_STR,lm_fd->libf);
  if (lm_fd->lfflag) getc(lm_fd->libf);

  if (lm_fd->lline[0]=='>'&&(lm_fd->lline[3]==';'|| lm_fd->lline[3]=='>')	/* GCG ascii */
      ) {
    SAFE_STRNCPY(str,&lm_fd->lline[4],cnt-1);

    if (lm_fd->lline[3]=='>' && (bp = strchr(str,' '))!=NULL) *bp='\0';

    if ((bp = strchr(str,':'))!=NULL) *bp='\0';
    if ((bp=strchr(str,'\r'))!=NULL) *bp='\0';
    else if ((bp = strchr(str,'\n'))!=NULL) *bp='\0';
    else str[cnt-1]='\0';

    fgets(lm_fd->lline,MAX_STR,lm_fd->libf);
    if (lm_fd->lfflag) getc(lm_fd->libf);

    /* skip over redundant stuff */
    for (llp=lm_fd->lline,bp=str; *llp==*bp; llp++,bp++);
    if ((int)(llp-lm_fd->lline)<5) llp = lm_fd->lline;

    if ((bp=strchr(llp,'\r'))!=NULL) *bp=' ';
    if ((bp=strchr(llp,'\n'))!=NULL) *bp='\0';
    SAFE_STRNCAT(str," ",cnt);
    SAFE_STRNCAT(str,llp,cnt);
  }
  else {
    str[0]='\0';
  }

  FSEEK(lm_fd->libf,seek,0);
  fgets(lm_fd->lline,MAX_STR,lm_fd->libf);
  if (lm_fd->lfflag) getc(lm_fd->libf);
}

static int gcg_bton[4]={2,4,1,3};

int
gcg_getlib(unsigned char *seq, int maxs,
	   char *libstr,
	   int n_libstr,
	   fseek_t *libpos,
	   int *lcont,
	   struct lmf_str *lm_fd,
	   long *l_off)
{
  char dummy[20];
  char gcg_date[10];
  register unsigned char *cp, *seqp, stmp;
  register int *ap;
  char gcg_type[10];
  unsigned char *seqm, *seqm1;
  long r_block, b_block;
  char *bp;

  *l_off = 1;

  seqp = seq;
  seqm = &seq[maxs-9];
  seqm1 = seqm-1;

  ap = lm_fd->sascii;

  if (*lcont==0) {
    while (lm_fd->lline[0]!='>' && lm_fd->lline[0]!=';') {
      lm_fd->lpos = FTELL(lm_fd->libf);
      if (fgets(lm_fd->lline,MAX_STR,lm_fd->libf)==NULL) return (-1);
    }
    sscanf(&lm_fd->lline[4],"%s %s %s %s %ld",
	   libstr,gcg_date,gcg_type,dummy,&(lm_fd->gcg_len));

    lm_fd->gcg_binary = (gcg_type[0]=='2');

    fgets(lm_fd->lline,MAX_STR,lm_fd->libf);
    while (strchr((char *)lm_fd->lline,'\n')==NULL) {
      if (strlen(lm_fd->lline)<MAX_STR/2) 
	fgets(&lm_fd->lline[strlen(lm_fd->lline)],MAX_STR/2,lm_fd->libf);
      else 
      fgets(&lm_fd->lline[strlen(lm_fd->lline)-MAX_STR/2],MAX_STR/2,lm_fd->libf);
    }
    lm_fd->lline[MAX_STR-1]='\0';
    if (n_libstr <= 21) {
      libstr[12]='\0';
    }
    else {
      SAFE_STRNCAT(libstr," ",n_libstr);
      SAFE_STRNCAT(libstr,lm_fd->lline,n_libstr);
      if ((bp = strchr(libstr,'\n'))!=NULL) *bp='\0';
      libstr[n_libstr-1]='\0';
    }
    *libpos = lm_fd->lpos;
  }

  lm_fd->lline[0]='\0';

  r_block = b_block = min((size_t)(seqm-seqp),lm_fd->gcg_len);
  if (lm_fd->gcg_binary) { r_block = (r_block+3)/4; }

  fread((char *)seqp,(size_t)r_block,(size_t)1,lm_fd->libf);
  if (!lm_fd->gcg_binary) 
    for (cp=seqp; seqp<seq+r_block; ) *seqp++ = ap[*cp++];
  else if (lm_fd->gcg_binary) {
    seqp = seq + r_block;
    cp = seq + 4*r_block;
    while (seqp > seq) {
      stmp = *--seqp;
      *--cp = gcg_bton[stmp&3];
      *--cp = gcg_bton[(stmp >>= 2)&3];
      *--cp = gcg_bton[(stmp >>= 2)&3];
      *--cp = gcg_bton[(stmp >>= 2)&3];
    }
  }
  if (4 * r_block >= lm_fd->gcg_len) {
    fgets(lm_fd->lline,MAX_STR,lm_fd->libf);
    *lcont = 0;
  }
  else {
    if (lm_fd->gcg_binary) b_block = 4*r_block;
    lm_fd->gcg_len -= b_block;
    (*lcont)++;
  }

  seq[b_block] = EOSEQ;
  /*   if (b_block==0) return 1; else */
  return b_block;
}

void
gcg_ranlib(char *str,
	   int cnt,
	   fseek_t seek,
	   char *libstr,
	   struct lmf_str *lm_fd)
{
  char *bp, *bp1, *llp;

  FSEEK(lm_fd->libf, seek, 0);
  fgets(lm_fd->lline,MAX_STR,lm_fd->libf);

  if (lm_fd->lline[0]=='>'&&(lm_fd->lline[3]==';'||lm_fd->lline[3]=='>')) {
    SAFE_STRNCPY(str,&lm_fd->lline[4],cnt-1);
    if ((bp = strchr(str,' '))!=NULL) *bp='\0';
    else if ((bp=strchr(str,'\r'))!=NULL) *bp='\0';
    else if ((bp = strchr(str,'\n'))!=NULL) *bp='\0';
    else str[cnt-1]='\0';

    fgets(lm_fd->lline,MAX_STR,lm_fd->libf);

    /* check beginning of line it is a duplicate */
    for (llp=lm_fd->lline,bp=str; *llp == *bp; llp++,bp++);
    if ((int)(llp-lm_fd->lline)<5) llp = lm_fd->lline;

    /* here we would like to skip over some species stuff */
	/*
    if ((bp1 = strchr(llp,';'))!=NULL && (int)(bp1-llp)<50) {
      if ((bp2 = strchr(bp1+1,';'))!=NULL && (int)(bp2-bp1)<50) {
	*(bp2+1)='\0'; bp1 = bp2+2;
      }
      else {bp1=llp;}
    }
    else if ((bp1=strchr(llp,'.'))!=NULL && *(bp1+1)==' ') {
      *(bp1+1) = '\0'; bp1 += 2;}
    else bp1 = llp;
    */
    
    bp1 = llp;
    if ((bp=strchr(bp1,'\r'))!=NULL) *bp='\0';
    if ((bp=strchr(bp1,'\n'))!=NULL) *bp='\0';
    SAFE_STRNCAT(str," ",cnt);
    SAFE_STRNCAT(str,bp1,cnt);
    if (bp1!=llp) SAFE_STRNCAT(str,llp,cnt);
  }
  else {
    str[0]='\0';
  }

  FSEEK(lm_fd->libf,seek,0);
  fgets(lm_fd->lline,MAX_STR,lm_fd->libf);
}

/* **************************************************************** */
/* the following section contains the functions used to initialize,
   read, hash, and lookup acc_lists, either from a sorted list, or a
   hash table.

   The functions are:

   void *sel_acc_libstr_init(FILE *libf, int *acc_off, char fmt_term)
   -- allocate space, read the accessions from an already open file of
      accessions

   int sel_acc_libstr(char *libstr, int gi, void *ptr)
   -- compare libstr to current accession, returning 1 and
      incrementing cur_entry if match is found.
      Requires sorted list; does not access gi

   int sel_acc_gi(char *libstr, int gi, void *ptr)
   -- compare gi to gi_list, returning 1 and incrementing if match.
      if (gi <= 0), get gi from libstr

   int sel_hacc_libstr(char *libstr, int gi, void *ptr)
   -- check to see whether libstr is in hash table

   int sel_hacc_gi(char *libstr, int gi, void *ptr)
   -- check to see if gi in hash32 table
*/
/* **************************************************************** */

struct sel_acc_str {
  int curr_entry;
  int max_entry;
  char fmt_term;
  char *acc_buff;
  char **acc_list;
  int *gi_list;
  int *acc_hash;
  int *acc_hash_link;
  int hash_mask;
};

/* allocate space, read the accessions from an already open file of
   accessions */

void *sel_acc_libstr_init(FILE *libf, int *acc_off, char fmt_term) {
  struct sel_acc_str *sel_acc_ptr;
  char acc_line[MAX_STR];
  char *bp, *bp1;
  char *acc_buff;
  char *acc_buff_max;	/* end of buffer */
  char *new_buff;	/* reallocated buffer size */
  char *acc_buff_p;
  char **acc_list;
  int acc_cnt, i;
  int new_buff_siz;
  int abuff_siz;		/* allocated buffer size */
  int buff_siz;		/* fread buff_siz */

  if ((sel_acc_ptr = (struct sel_acc_str *)calloc(1,sizeof(struct sel_acc_str)))==NULL) {
    fprintf(stderr, "*** Warning [%s:%d] - Cannot allocate struct sel_acc_str\n",__FILE__, __LINE__);
    return NULL;
  }

  /*
  if (fmt && *fmt != '\0') {
    sel_acc_ptr->fmt = (char *)calloc(strlen(fmt)+1,sizeof(char));
    SAFE_STRNCPY(sel_acc_ptr->fmt, fmt, strlen(fmt)+1);
  }
  */

  sel_acc_ptr->fmt_term = fmt_term;

  /* allocate some space for the ACC's */

  abuff_siz = new_buff_siz = 640000;

  if ((acc_buff = (char *)calloc(abuff_siz*10, sizeof(char)))==NULL) {
    fprintf(stderr, "*** Warning [%s:%d] - Cannot allocate acc buff %d\n",__FILE__, __LINE__,abuff_siz*10);
    free(sel_acc_ptr);
    return NULL;
  }

  /* iteratively read and reallocate space for buffer until its all read */

  acc_buff_p = acc_buff;
  while ((buff_siz = fread(acc_buff_p, sizeof(char), new_buff_siz, libf))==new_buff_siz) {
    if ((new_buff = realloc(acc_buff, (size_t)(abuff_siz+new_buff_siz)))==NULL) {
      fprintf(stderr, "*** Warning [%s:%d] - cannot reallocate for acc_buf[%d]\n",__FILE__, __LINE__,abuff_siz);
      break;
    }
    else {
      acc_buff = new_buff;
      acc_buff_p = acc_buff + abuff_siz;
      abuff_siz += new_buff_siz;
    }
  }
  fclose(libf);

  acc_buff_max = acc_buff_p + buff_siz;

  /* convert all the ACC lines (with \n) to null-terminated and
     count the number of aacc's */

  acc_cnt = 0;
  acc_buff_p = acc_buff;
  while (acc_buff_p < acc_buff_max && (bp = strchr(acc_buff_p,'\n'))!=NULL) {
    *bp = '\0';
    /*  also remove '\r'); */
    if ((bp1=strchr(acc_buff_p,'\r'))!=NULL) {*bp1 = '\0';}
    acc_cnt++;
    acc_buff_p = bp+1;
  }

  /* allocate the acc_list */
  if ((acc_list=(char **)calloc(acc_cnt+1, sizeof(char **)))==NULL) {
    fprintf(stderr,"*** Warning [%s:%d] - cannot allocate acc_list[%d]\n",__FILE__, __LINE__,acc_cnt+1);
    free(sel_acc_ptr);
    return NULL;
  }

  /* now load acc_list[] */
  for (i=0, acc_buff_p=acc_buff; i<acc_cnt; i++) {
    acc_list[i] = acc_buff_p;
    acc_buff_p += strlen(acc_buff_p)+1;
  }

  /* finally put everything in the structure to be returned */
  sel_acc_ptr->acc_buff = acc_buff;
  sel_acc_ptr->acc_list = acc_list;
  sel_acc_ptr->curr_entry = 0;
  sel_acc_ptr->max_entry = acc_cnt;
  return (void *)sel_acc_ptr;
}

int sel_acc_libstr(char *libstr, int gi, void *ptr) {
  struct sel_acc_str *sel_acc_ptr;
  char *curr_acc;
  char acc[MAX_SSTR], *acc_p, *bp;

  sel_acc_ptr = (struct sel_acc_str *)ptr;
  
  if (sel_acc_ptr->curr_entry >= sel_acc_ptr->max_entry) return -1;

  if ((bp = strchr(libstr,sel_acc_ptr->fmt_term))!=NULL) {
    *bp = '\0';
  }

  curr_acc = sel_acc_ptr->acc_list[sel_acc_ptr->curr_entry];
  if (libstr[2] == curr_acc[2] && libstr[1] == curr_acc[1] && 
      strncmp(libstr,curr_acc,MAX_UID)==0) {
    sel_acc_ptr->curr_entry++;
    return 1;
  }
  else {
    return 0;
  }
}

void *sel_acc_gi_init(FILE *libf, int *acc_off, char fmt_term) {
  struct sel_acc_str *sel_acc_ptr;
  char acc_line[MAX_STR];
  char *bp;
  int *gi_list;
  int *new_buff;	/* reallocated buffer size */
  int *acc_buff_p;
  int acc_cnt, i;
  int new_buff_siz;
  int abuff_siz;		/* allocated buffer size */
  int buff_siz;		/* fread buff_siz */

  if ((sel_acc_ptr = (struct sel_acc_str *)calloc(1,sizeof(struct sel_acc_str)))==NULL) {
    fprintf(stderr, "*** Warning [%s:%d] - Cannot allocate struct sel_acc_str\n",__FILE__, __LINE__);
    return NULL;
  }

  /*
  if (fmt != NULL && *fmt != '\0') {
    sel_acc_ptr->fmt = (char *)calloc(strlen(fmt)+1,sizeof(char));
    SAFE_STRNCPY(sel_acc_ptr->fmt, fmt, strlen(fmt)+1);
  }
  else {
    sel_acc_ptr->fmt = NULL;
  }
  */

  sel_acc_ptr->fmt_term = fmt_term;

  /* now allocate some space for the ACC's */

  abuff_siz = new_buff_siz = 64000;

  if ((gi_list = (int *)calloc(abuff_siz, sizeof(int)))==NULL) {
    fprintf(stderr, "*** Warning [%s:%d] - Cannot allocate acc buff %d\n",__FILE__, __LINE__,abuff_siz);
    free(sel_acc_ptr);
    return NULL;
  }

  /* now iteratively read and reallocate space for buffer until its all read */

  acc_cnt = 0;
  while (fgets(acc_line, sizeof(acc_line), libf)!=NULL) {
    gi_list[acc_cnt++] = atoi(acc_line);
    if (acc_cnt >= abuff_siz) {
      if ((new_buff = realloc(gi_list,(abuff_siz + new_buff_siz)*sizeof(int)))==NULL) {
	fprintf(stderr,"*** Warning [%s:%d] - cannot realloc gi_list[%d]\n",__FILE__, __LINE__,abuff_siz+new_buff_siz);
	break;
      }
      else {
	abuff_siz += new_buff_siz;
	gi_list = new_buff;
      }
    }
  }
  fclose(libf);

  /* finally put everything in the structure to be returned */
  sel_acc_ptr->gi_list = gi_list;
  sel_acc_ptr->curr_entry = 0;
  sel_acc_ptr->max_entry = acc_cnt;
  return (void *)sel_acc_ptr;
}

int sel_acc_gi(char *libstr, int gi, void *ptr) {
  struct sel_acc_str *sel_acc_ptr;
  char *bp;

  sel_acc_ptr = (struct sel_acc_str *)ptr;
  
  if (sel_acc_ptr->curr_entry >= sel_acc_ptr->max_entry) return -1;

  if (gi <= 0) {
    if (libstr) {
      /*
      if (sel_acc_ptr->fmt) {
	sscanf(libstr,sel_acc_ptr->fmt,&gi);
      }
      */
      gi = atoi(libstr);
    }
  }

  if (gi == sel_acc_ptr->gi_list[sel_acc_ptr->curr_entry]) {
    sel_acc_ptr->curr_entry++;
    return 1;
  }
  else {
    return 0;
  }
}

/* this version of the selection algorithm does not require sorted
   lists.  It hashes the initial list, and uses the hash table to
   lookup the library sequences.

*/

#define HASH_TABLE_MULT 8

void *sel_hacc_libstr_init(FILE *libf, int *acc_off, char fmt_term) {
  struct sel_acc_str *sel_acc_ptr;
  char acc_line[MAX_STR];
  char *bp, *bp1;
  char *acc_buff;
  char *acc_buff_max;	/* end of buffer */
  char *new_buff;	/* reallocated buffer size */
  char *acc_buff_p;
  char **acc_list;
  int acc_cnt, i;
  int new_buff_siz;
  int abuff_siz;		/* allocated buffer size */
  int buff_siz;		/* fread buff_siz */
  int hash_mask, hash_max;
  int hash_val, *acc_hash, *acc_hash_link;
  int link_save;

  if ((sel_acc_ptr = (struct sel_acc_str *)calloc(1,sizeof(struct sel_acc_str)))==NULL) {
    fprintf(stderr, "*** Warning [%s:%d] - Cannot allocate struct sel_acc_str\n",__FILE__, __LINE__);
    return NULL;
  }

  /*
  if (fmt && *fmt != '\0') {
    sel_acc_ptr->fmt = (char *)calloc(strlen(fmt)+1,sizeof(char));
    SAFE_STRNCPY(sel_acc_ptr->fmt, fmt, strlen(fmt)+1);
  }
  */

  sel_acc_ptr->fmt_term = fmt_term;

  /* now allocate some space for the ACC's */

  abuff_siz = new_buff_siz = 640000;

  if ((acc_buff = (char *)calloc(abuff_siz, sizeof(char)))==NULL) {
    fprintf(stderr, "*** Warning [%s:%d] - Cannot allocate acc buff %d\n",__FILE__, __LINE__,abuff_siz);
    free(sel_acc_ptr);
    return NULL;
  }

  /* now iteratively read and reallocate space for buffer until its all read */

  acc_buff_p = acc_buff;
  while ((buff_siz = fread(acc_buff_p, sizeof(char), new_buff_siz, libf))==new_buff_siz) {
    if ((new_buff = realloc(acc_buff, (size_t)(abuff_siz+new_buff_siz)))==NULL) {
      fprintf(stderr, "*** Warning [%s:%d] - cannot reallocate for acc_buf[%d]\n",__FILE__, __LINE__,abuff_siz);
      break;
    }
    else {
      acc_buff = new_buff;
      acc_buff_p = acc_buff + abuff_siz;
      abuff_siz += new_buff_siz;
    }
  }
  fclose(libf);

  acc_buff_max = acc_buff_p + buff_siz;

  /* now convert all the ACC lines (with \n) to null-terminated and
     count the number of aacc's */

  acc_cnt = 0;
  acc_buff_p = acc_buff;
  while (acc_buff_p < acc_buff_max && (bp = strchr(acc_buff_p,'\n'))!=NULL) {
    *bp = '\0';
    /*  also remove '\r'); */
    if ((bp1=strchr(acc_buff_p,'\r'))!=NULL) {*bp1 = '\0';}
    acc_cnt++;
    acc_buff_p = bp+1;
  }

  /* allocate the acc_list */
  if ((acc_list=(char **)calloc(acc_cnt+1, sizeof(char **)))==NULL) {
    fprintf(stderr,"*** Warning [%s:%d] - cannot allocate acc_list[%d]\n",__FILE__, __LINE__,acc_cnt+1);
    free(sel_acc_ptr);
    return NULL;
  }

  /* now load acc_list[] */
  for (i=0, acc_buff_p=acc_buff; i<acc_cnt; i++) {
    acc_list[i] = acc_buff_p;
    acc_buff_p += strlen(acc_buff_p)+1;
  }

  /* allocate the hash for the acc_list - we want a table that is
     about 4X acc_cnt and a power of 2 */
  for (hash_max = 4096; hash_max <= HASH_TABLE_MULT * acc_cnt; hash_max *= 2);
  hash_mask = hash_max - 1;

  if ((acc_hash = (int *)calloc(hash_max, sizeof(int)))==NULL) {
    fprintf(stderr, "*** Warning [%s:%d] - cannot allocate acc_hash[%ld]\n",__FILE__, __LINE__,hash_max*sizeof(int));
    return NULL;
  }

  /* allocate the acc_list link table */
  if ((acc_hash_link = (int *)calloc(acc_cnt+1,sizeof(char *)))==NULL) {
    fprintf(stderr, "*** Warning [%s:%d] - cannot allocate acc_hash_link[%ld]\n",__FILE__, __LINE__,acc_cnt*sizeof(char *));
    return NULL;
  }

  for (i=0; i<acc_cnt; i++) {
    hash_val = hash_func(acc_list[i]) & hash_mask;
    if ((link_save = acc_hash[hash_val]) != 0) {
      acc_hash_link[i+1]=link_save;
    }
    acc_hash[hash_val] = i+1;
  }

  /* finally put everything in the structure to be returned */
  sel_acc_ptr->acc_buff = acc_buff;
  sel_acc_ptr->acc_list = acc_list;
  sel_acc_ptr->curr_entry = 0;
  sel_acc_ptr->max_entry = acc_cnt;
  /* hash stuff */
  sel_acc_ptr->acc_hash = acc_hash;
  sel_acc_ptr->acc_hash_link = acc_hash_link;
  sel_acc_ptr->hash_mask = hash_mask;

  return (void *)sel_acc_ptr;
}

int sel_hacc_libstr(char *libstr, int gi, void *ptr) {
  int i, i1;
  struct sel_acc_str *sel_acc_ptr;
  char *bp,  **acc_list;
  int hash_val;

  sel_acc_ptr = (struct sel_acc_str *)ptr;
  
  if (sel_acc_ptr->curr_entry >= sel_acc_ptr->max_entry) return -1;

  if ((bp = strchr(libstr,sel_acc_ptr->fmt_term))!=NULL) {
    *bp = '\0';
  }

  hash_val = hash_func(libstr) & sel_acc_ptr->hash_mask;

  if (sel_acc_ptr->acc_hash[hash_val] == 0) return 0;

  acc_list = sel_acc_ptr->acc_list;

  for (i=sel_acc_ptr->acc_hash[hash_val]; i > 0;
       i = sel_acc_ptr->acc_hash_link[i]) {
    i1 = i-1;
    if (libstr[2] == acc_list[i1][2] && libstr[1] == acc_list[i1][1] && 
      strncmp(libstr,acc_list[i1],MAX_UID)==0) {
      return 1;
    }
  }
  return 0;
}

void *sel_hacc_gi_init(FILE *libf, int *acc_off, char fmt_term) {
  struct sel_acc_str *sel_acc_ptr;
  char acc_line[MAX_STR];
  char *bp;
  int *gi_list;
  int *new_buff;	/* reallocated buffer size */
  int *acc_buff_p;
  int acc_cnt, i;
  int new_buff_siz;
  int abuff_siz;		/* allocated buffer size */
  int buff_siz;		/* fread buff_siz */
  int hash_val, hash_max, hash_mask, link_save;
  int *gi_hash, *gi_hash_link;

  if ((sel_acc_ptr = (struct sel_acc_str *)calloc(1,sizeof(struct sel_acc_str)))==NULL) {
    fprintf(stderr, "*** Warning [%s:%d] - Cannot allocate struct sel_acc_str\n",__FILE__, __LINE__);
    return NULL;
  }

  /*
  if (fmt != NULL && *fmt != '\0') {
    sel_acc_ptr->fmt = (char *)calloc(strlen(fmt)+1,sizeof(char));
    SAFE_STRNCPY(sel_acc_ptr->fmt, fmt, strlen(fmt)+1);
  }
  else {
    sel_acc_ptr->fmt = NULL;
  }
  */

  sel_acc_ptr->fmt_term = fmt_term;

  /* now allocate some space for the ACC's */

  abuff_siz = new_buff_siz = 64000;

  if ((gi_list = (int *)calloc(abuff_siz, sizeof(int)))==NULL) {
    fprintf(stderr, "*** Warning [%s:%d] - Cannot allocate acc buff %d\n",__FILE__, __LINE__,abuff_siz);
    free(sel_acc_ptr);
    return NULL;
  }

  /* now iteratively read and reallocate space for buffer until its all read */

  acc_cnt = 0;
  while (fgets(acc_line, sizeof(acc_line), libf)!=NULL) {
    gi_list[acc_cnt++] = atoi(acc_line);
    if (acc_cnt >= abuff_siz) {
      if ((new_buff = realloc(gi_list,(abuff_siz + new_buff_siz)*sizeof(int)))==NULL) {
	fprintf(stderr,"*** Warning [%s:%d] - cannot realloc gi_list[%d]\n",__FILE__, __LINE__,abuff_siz+new_buff_siz);
	break;
      }
      else {
	abuff_siz += new_buff_siz;
	gi_list = new_buff;
      }
    }
  }
  fclose(libf);

  /* allocate the hash for the gi_list - we want a table that is
     about 4X acc_cnt and a power of 2 */
  for (hash_max = 4096; hash_max <= HASH_TABLE_MULT * acc_cnt; hash_max *= 2);
  hash_mask = hash_max - 1;

  if ((gi_hash = (int *)calloc(hash_max, sizeof(int)))==NULL) {
    fprintf(stderr, "*** Warning [%s:%d] - cannot allocate gi_hash[%ld]\n",__FILE__, __LINE__,hash_max*sizeof(int));
    return NULL;
  }

  /* allocate the gi_list link table */
  if ((gi_hash_link = (int *)calloc(acc_cnt+1,sizeof(char *)))==NULL) {
    fprintf(stderr, "*** Warning [%s:%d] - cannot allocate gi_hash_link[%ld]\n",__FILE__, __LINE__,acc_cnt*sizeof(char *));
    return NULL;
  }

  for (i=0; i<acc_cnt; i++) {
    hash_val = fast_hash32(gi_list[i]) & hash_mask;
    if (gi_hash[hash_val] != 0) {
      link_save = gi_hash[hash_val];
      gi_hash_link[i+1]=link_save;
    }
    gi_hash[hash_val] = i+1;
  }

  /* finally put everything in the structure to be returned */
  sel_acc_ptr->gi_list = gi_list;
  sel_acc_ptr->curr_entry = 0;
  sel_acc_ptr->max_entry = acc_cnt;

  sel_acc_ptr->acc_hash = gi_hash;
  sel_acc_ptr->acc_hash_link = gi_hash_link;
  sel_acc_ptr->hash_mask = hash_mask;

  return (void *)sel_acc_ptr;
}

int sel_hacc_gi(char *libstr, int gi, void *ptr) {
  struct sel_acc_str *sel_acc_ptr;
  int hash_val;
  int *gi_list, i;
  char *bp;

  sel_acc_ptr = (struct sel_acc_str *)ptr;
  
  if (sel_acc_ptr->curr_entry >= sel_acc_ptr->max_entry) return -1;

  if (gi <= 0) {
    if (libstr) {
      /*
      if (sel_acc_ptr->fmt) {
	sscanf(libstr,sel_acc_ptr->fmt,&gi);
      }
      */
      gi = atoi(libstr);
    }
  }

  hash_val = fast_hash32(gi) & sel_acc_ptr->hash_mask;
  if (sel_acc_ptr->acc_hash[hash_val]==0) return 0;

  gi_list = sel_acc_ptr->gi_list;

  for (i=sel_acc_ptr->acc_hash[hash_val]; i > 0;
       i = sel_acc_ptr->acc_hash_link[i]) {
    if (gi_list[i-1] == gi) return 1;
  }
  return 0;
}


/* adapted from	http://burtleburtle.net/bob/hash/doobs.html */
unsigned int
hash_func(char *key)
{
  unsigned int hash;

  hash = 0;

  while (*key) {
    hash += *key++;
    hash += (hash << 10);
    hash ^= (hash >> 6);
  }

  hash += (hash << 3);
  hash ^= (hash >> 11);
  hash += (hash << 15);

  return hash;
}

unsigned int
fast_hash32 (unsigned int data) {
  int tmp, hash;


  hash  = data >> 16;
  tmp    = ((data & 0xFFFF) << 11) ^ hash;
  hash   = (hash << 16) ^ tmp;
  hash  += hash >> 11;

  /* Force "avalanching" of final 127 bits */
  hash ^= hash << 3;
  hash += hash >> 5;
  hash ^= hash << 4;
  hash += hash >> 17;
  hash ^= hash << 25;
  hash += hash >> 6;

  return hash;
}

/* takes a file name string and allocates space for it, returning
   pointer to space */
char *alloc_file_name(char *f_name) {
  int fn_len;
  char *alloc_f_name;

  fn_len = strlen(f_name);
  if ((alloc_f_name = calloc(fn_len+1,sizeof(char)))==NULL) {
    fprintf(stderr, "*** Warning [%s:%d] - Cannot allocate %d space for %s\n",
	    __FILE__, __LINE__,fn_len+1, f_name);
    exit(1);
  }
  else {
    SAFE_STRNCPY(alloc_f_name,f_name,fn_len);
    return alloc_f_name;
  }
}
