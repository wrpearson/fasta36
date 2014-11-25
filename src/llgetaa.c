/* $Id: llgetaa.c 625 2011-03-23 17:21:38Z wrp $ */

/* copyright (c) 1996, 1997, 1998, 1999, 2007 by William R. Pearson and
   The Rector & Visitors of the University of Virginia */

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

/*
   Feb, 1998 - version for prss 

   March, 2001 - modifications to support comp_thr.c: use libpos to indicate
   whether the score is shuffled==1 or unshuffled==0.  This simplifies
   complib.c and makes comp_thr.c possible

   modified version of nxgetaa.c that generates random sequences
   for a library
*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "defs.h"
#include "mm_file.h"

#include "uascii.h"
#include "structs.h"

#define XTERNAL
#include "upam.h"
#undef XTERNAL

#define YES 1
#define NO 0
#define MAXLINE 512

#ifndef min
#define min(x,y) ((x) > (y) ? (y) : (x))
#endif

static int use_stdin=0;
static char llibstr0[256];
static char llibstr1[256];
static char o_line[256];

#define NO_FORMAT 0
#define FASTA_FORMAT 1
#define GCG_FORMAT 2
static int seq_format=NO_FORMAT;
static char seq_title[200];

extern int irand(int);
extern void shuffle(unsigned char *from, unsigned char *to, int n);
extern void wshuffle(unsigned char *from, unsigned char *to, int n, int wsiz, int *ieven);

int
getseq(char *filen, int *qascii,
       unsigned char *seq, int maxs, char *libstr,
       int n_libstr, long *sq0off)
{
  FILE *fptr;
  char line[512],*bp;
  int i, j, n;
  int ic;
  int sstart, sstop, sset=0;
  int have_desc = 0;
  int desc_complete = 0;
  int llen, l_offset;

  seq_title[0]='\0';

  sstart = sstop = -1;
#ifndef DOS
  if ((bp=strchr(filen,':'))!=NULL) {
#else
  if ((bp=strchr(filen+3,':'))!=NULL) {
#endif
    *bp='\0';
    if (*(bp+1)=='-') sscanf(bp+2,"%d",&sstop);
    else sscanf(bp+1,"%d-%d",&sstart,&sstop);
    sset=1;
  }

  if (strcmp(filen,"-") && strcmp(filen,"@")) {
    if ((fptr=fopen(filen,"r"))==NULL) {
      fprintf(stderr," could not open %s\n",filen);
      return 0;
    }
  }
  else {
    fptr = stdin;
    use_stdin++;
  }

  if (use_stdin > 1) {
    have_desc = 1;
    if ((bp=strchr(o_line,'\001'))!=NULL) *bp='\0';
    strncpy(llibstr1,o_line,sizeof(llibstr1));
    strncpy(libstr,o_line,n_libstr);
    libstr[n_libstr-1]='\0';
    l_offset = 0;
  }

  if (sset==1) {
    filen[strlen(filen)]=':';
    if (*sq0off==1 || sstart>1) *sq0off = sstart;
  }

  desc_complete = 0;
  n=0;
  while(fgets(line,sizeof(line),fptr)!=NULL) {
    if (line[0]=='>') {
      if (have_desc) {
	strncpy(o_line,line,sizeof(o_line));
	goto last;
      }
      l_offset = 0;
      seq_format = FASTA_FORMAT;

      if ((bp=(char *)strchr(line,'\n'))!=NULL) {
	*bp='\0';				/* have newline */
	desc_complete = 1;
      }

      if ((bp=strchr(line+1,'\001'))!=NULL) *bp='\0';
      strncpy(seq_title,line+1,sizeof(seq_title));
      strncpy(llibstr0,line+1,sizeof(llibstr0));
      if (n_libstr <= 20) {
	if ((bp=(char *)strchr(line,' '))!=NULL) *bp='\0';
      }
      strncpy(libstr,line+1,n_libstr);
      libstr[n_libstr-1]='\0';

      if (!desc_complete) {
	while (fgets(line, sizeof(line), fptr) != NULL) {
	  if (strchr(line,'\n') != NULL) {
	    line[0]='>';
	    break;
	  }
	}
	desc_complete = 1;
      }
    }
    else if (seq_format==NO_FORMAT) {
      seq_format = GCG_FORMAT;
      qascii['*'] = qascii['X'];
      l_offset = 10;
      llen = strlen(line);
      while (strncmp(&line[llen-3],"..\n",(size_t)3) != 0) {
	if (fgets(line,sizeof(line),fptr)==NULL) return 0;
	llen = strlen(line);
      }
      if (n_libstr <= 20) {
	if ((bp=(char *)strchr(line,' '))!=NULL) *bp='\0';
	else if ((bp=(char *)strchr(line,'\n'))!=NULL) *bp='\0';
      }
      strncpy(libstr,line,n_libstr);
      libstr[n_libstr-1]='\0';
      if (fgets(line,sizeof(line),fptr)==NULL) return 0;
    }

    if (seq_format==GCG_FORMAT && strlen(line)<l_offset) continue;

    if (line[0]!='>'&& line[0]!=';') {
      for (i=l_offset; (n<maxs)&&
	     ((ic=qascii[line[i]&AAMASK])<EL); i++)
	if (ic<NA) seq[n++]= ic;
      if (ic == ES) break;
    }
    else {
      if (have_desc) {
	strncpy(o_line,line,sizeof(o_line));
	goto last;
      }
      else {
	have_desc = 1;
      }
    }
  }

 last:
  if (n==maxs) {
    fprintf(stderr," sequence may be truncated %d %d\n",n,maxs);
    fflush(stderr);
  }
  if ((bp=strchr(libstr,'\n'))!=NULL) *bp = '\0';
  if ((bp=strchr(libstr,'\r'))!=NULL) *bp = '\0';
  seq[n]= EOSEQ;

  if (fptr!=stdin) fclose(fptr);

  if (sset) {
    if (sstart <= 0) sstart = 1;
    if (sstop <= 0) sstop = n;
    sstart--;
    sstop--;
    for (i=0, j=sstart; j<=sstop; i++,j++)
      seq[i] = seq[j];
    n = sstop - sstart +1;
    seq[n]=EOSEQ;
  }

  return n;
}

int
gettitle(filen,title,len)
  char *filen, *title; int len;
{
  FILE *fptr;
  char line[512];
  char *bp;
  int ll,sset;
#ifdef WIN32
  char *strpbrk();
#endif
  sset = 0;

  if (use_stdin) {
    if (use_stdin == 1) {
      /*      use_stdin++; */
      strncpy(title,llibstr0,len);
    }
    else {
      strncpy(title,llibstr1,len);
    }
    if ((bp=strchr(title,'\001'))!=NULL) *bp='\0';
    return strlen(title);
  }

  if ((bp=strchr(filen,':'))!=NULL) { *bp='\0'; sset=1;}
	  
  if ((fptr=fopen(filen,"r"))==NULL) {
    fprintf(stderr," file %s was not found\n",filen);
    fflush(stderr);
    return 0;
  }

  if (sset==1) filen[strlen(filen)]=':';

  while(fgets(line,sizeof(line),fptr)!=0) {
    if (line[0]=='>'|| line[0]==';') goto found;
  }
  fclose(fptr);
  title[0]='\0';
  return 0;

 found:
  if ((bp=strchr(line,'\001'))!=NULL) *bp = 0;
#ifdef WIN32
  bp = strpbrk(line,"\n\r");
#else
  bp = strchr(line,'\n');
#endif
  if (bp!=NULL) *bp = 0;
  strncpy(title,line,len);
  title[len-1]='\0';
  fclose(fptr);
  return strlen(title);
}	

FILE *libf=NULL;

long lpos;
char lline[MAXLINE];
int lfflag=0;	/* flag for CRLF in EMBL CDROM files */
#define LFCHAR '\015'  /* for MWC 5.5 */

int agetlib(); void aranlib();	/* pearson fasta format */

/*	the following is from fgetgb.c */

/* a file name for open_lib may now include a library type suffix */
/* only opens fasta format files */

static char libn_save[MAX_FN];
static int ldna_save=0;
static int do_shuffle;
static int shuff_cnt=10;
static int w_flag = 0;
#ifdef DEBUG
static FILE *dfile=NULL;
#endif
static unsigned char *aa_save;
static int n1_save;
static int i_even;

/* lmf_str * is used here for compatibility with the "normal" open_lib,
   but is largely unnecessary */

void 
set_shuffle(struct mngmsg m_msg) {
  char dfname[MAX_FN];

  if (m_msg.shuff_wid > 0) w_flag = m_msg.shuff_wid;
  if (m_msg.shuff_max > shuff_cnt) shuff_cnt = m_msg.shuff_max;

#ifdef DEBUG
  if (m_msg.dfile[0]!='\0') {
    strncpy(dfname,m_msg.dfile,sizeof(dfname));
    strncat(dfname,"_rlib",sizeof(dfname));
    dfile = fopen(dfname,"w");
  }
#endif
}

struct lmf_str *
open_lib(char *lname, int ldnaseq, int *sascii, int quiet, struct lmf_str *m_fd)
{
  char rline[10],libn[MAX_FN], *bp;
  int wcnt, ll, opnflg;
  int libtype;
  struct lmf_str *m_fptr;

  wcnt = 0;
  libtype = 0;

  strncpy(libn_save,lname,sizeof(libn_save));

  /* now allocate a buffer for the opened text file */
  if ((m_fptr = calloc(1,sizeof(struct lmf_str)))==NULL) {
    fprintf(stderr," cannot allocate lmf_str (%ld) for %s\n",
	    sizeof(struct lmf_str),lname);
    return NULL;
  }

  strncpy(m_fptr->lb_name,lname,MAX_FN);
  m_fptr->lb_name[MAX_FN-1]='\0';

  m_fptr->sascii = sascii;
  m_fptr->getlib = agetlib;
  m_fptr->ranlib = aranlib;
  m_fptr->mm_flg = 0;

  do_shuffle = 0;
  irand(0);		/* initialize the random number generator */

  return m_fptr;
}

void
closelib()
{
  if (libf!=NULL) {
    fclose(libf);
    libf = NULL;
  }
#ifdef DEBUG
  if (dfile) fclose(dfile);
#endif
}

static int ieven=0;
static char *desc_save;

int
agetlib(unsigned char *seq, 
	int maxs,
	char *libstr,
	int n_libstr,
	fseek_t *libpos,
	int *lcont, 
	struct lmf_str *lf_fd,
	long *l_off)
{
  long sq1_off;
  char lib_desc[120];
  int i;

  *l_off = 1;

  if (!do_shuffle) {
    do_shuffle = 1;
    
    if ((n1_save = getseq(libn_save,lf_fd->sascii,
			  seq,maxs,lib_desc,sizeof(lib_desc),&sq1_off)) < 1)
      return n1_save;

    strncpy(libstr,lib_desc,n_libstr);
    libstr[n_libstr-1]='\0';

    if ((aa_save = (unsigned char *)calloc(n1_save+1,sizeof(unsigned char)))==
	NULL) fprintf(stderr," cannot allocate %d for saved sequence\n",
		       n1_save);
    memcpy((void *)aa_save,(void *)seq,n1_save);

    if ((desc_save =
	 (char *)calloc(strlen(lib_desc)+1,sizeof(char)))== NULL) {
      fprintf(stderr," cannot allocate saved desciption [%d]\n",
	      strlen(lib_desc)+1);
    }
    else {
      strncpy (desc_save,lib_desc,strlen(lib_desc));
      desc_save[strlen(lib_desc)]=='\0';
    }

    *libpos = 0;
    return n1_save;
  }
  else {	/* return a shuffled sequence - here we need a window size; */
    strncpy(libstr,desc_save,n_libstr);
    libstr[n_libstr-1]='\0';

    if (shuff_cnt-- <= 0 ) return -1;
    if (w_flag > 0) wshuffle(aa_save,seq,n1_save,w_flag,&ieven);
    else shuffle(aa_save,seq,n1_save);
    seq[n1_save] = EOSEQ;
#ifdef DEBUG
    if (dfile!=NULL) {
      fprintf(dfile,">%d\n",shuff_cnt);
      for (i=0; i<n1_save; i++) {
	if (aa[seq[i]]>0) fputc(aa[seq[i]],dfile);
	else {fprintf(stderr,"error aa0[%d]: %d %d\n",
		      i,seq[i],aa[seq[i]]);}
	if (i%60 == 59) fputc('\n',dfile);
      }
      fputc('\n',dfile);
    }
#endif
    *libpos = 1;
    return n1_save;
  }
}

void
aranlib(char *str,
	int cnt,
	fseek_t seek,
	char *libstr,
	struct lmf_str *lm_fd)
{
  char *bp;
  int ll;

  if (use_stdin == 2) {
    if (llibstr1[0]=='>' || llibstr1[0]==';') {
      strncpy(str,llibstr1+1,cnt);
    }
    else {
      strncpy(str,llibstr1,cnt);
    }
  }
  else {
    strncpy(str,desc_save,cnt);
  }
  str[cnt-1]='\0';
  if ((bp = strchr(str,'\001'))!=NULL) *bp='\0';
  else if ((bp = strchr(str,'\n'))!=NULL) *bp='\0';
  else str[cnt-1]='\0';
}

/*
void
revcomp(unsigned char *seq, int n, int *c_nt)
{
  unsigned char tmp;
  int i, ni;


  for (i=0, ni = n-1; i< n/2; i++,ni--) {
    tmp = c_nt[seq[i]];
    seq[i] = c_nt[seq[ni]];
    seq[ni] = tmp;
  }
  if ((n%2)==1) {
    i = n/2;
    seq[i] = c_nt[seq[i]];
  }
}
*/

struct lmf_str *
re_openlib(struct lmf_str *om_fptr, int outtty)
{
  return om_fptr;
}

int re_getlib(unsigned char *aa1, int n1, int maxt3, int loff, int cont,
	      int term_code, long *loffset, long *l_off,
	      struct lmf_str *m_file_p)
{
  *loffset = 0;
  *l_off = 1;
  return n1;
}

