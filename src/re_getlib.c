/* re_getlib.c - re-acquire a sequence given lseek, lcont */

/* $Id: re_getlib.c 1227 2013-09-26 19:19:28Z wrp $  */

/* copyright (C) 2005, 2008, 2014 by William R. Pearson and The Rector &
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

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "defs.h"
#include "structs.h"
#include "mm_file.h"
#define XTERNAL
#include "uascii.h"

#define GETLIB (m_fptr->getlib)

/* modified Feb, 2008 to provide aa1a - annotation string */
extern int ann_scan(unsigned char *aa0, int n0, unsigned char **aa0a_p, int seqtype);

int
re_getlib(unsigned char *aa1,
	  struct annot_str **annot_p,
	  int maxn,	/* longest aa1 */
	  int maxt3,	/* alternate maxn */
	  int loff,	/* overlap */
	  int lcont,
	  int term_code,
	  long *loffset,	/* offset from real start of sequence */
	  long *l_off_p,	/* coordinate of sequence start */
	  struct lmf_str *m_fptr) {

  unsigned char *aa1ptr;
  int *sascii_save;
  int icont, maxt, ccont, n1;
  char libstr[MAX_UID];
  fseek_t lmark; 
  int sstart, sstop, is, id;
  
  aa1ptr = aa1;
  icont=0;

  /* no longer do selection */
  m_fptr->sel_acc_p = NULL;

  *loffset = 0l;
  maxt = maxn;
  n1 = -1;

  /* to process sequences in pieces properly, if lcont > 0, then we
     must read all but the last sequence using the scanning sascii,
     and then read the last piece using the ann_ascii */

  if (lcont > 1) {
    for (ccont=0; ccont<lcont-1; ccont++) {

      n1= GETLIB(aa1ptr,maxt,libstr,sizeof(libstr),&lmark,&icont,m_fptr,l_off_p);

      if (term_code && m_fptr->lib_aa && aa1ptr[n1-1]!=term_code) {
	aa1ptr[n1++]=term_code;
	aa1ptr[n1]=0;
      }

      if (aa1ptr!=aa1) n1 += loff;

      if (icont) {
	maxt = maxt3;
	memcpy(aa1,&aa1[n1-loff],loff);
	aa1ptr= &aa1[loff];
	*loffset += n1 - loff;
      }
      else {
	maxt = maxn;
	aa1ptr=aa1;
      }
    }
  }

  /* for the last one, replace m_fptr->sascii with ann_ascii[], and
     read the sequence */

  /* change sascii matrix only if there are annotations - otherwise
     l_ann_ascii is not initialized */
  /* cannot scan for annotations in memory mapped files */
  if (annot_p != NULL && !m_fptr->mm_flg) {
    if (*annot_p || (*annot_p = (struct annot_str *)calloc(1,sizeof(struct annot_str)))!=NULL) {
      sascii_save = m_fptr->sascii;
      m_fptr->sascii = l_ann_ascii;
      n1= GETLIB(aa1ptr,maxt,libstr,sizeof(libstr),&lmark,&icont,m_fptr,l_off_p);
      m_fptr->sascii = sascii_save;
      n1 = ann_scan(aa1ptr,n1,&((*annot_p)->aa1_ann),0);
    }
    else {
      fprintf(stderr,"re_getlib.c: cannot allocate annot_p\n");
      n1= GETLIB(aa1ptr,maxt,libstr,sizeof(libstr),&lmark,&icont,m_fptr,l_off_p);
    }
  }
  else {
    n1= GETLIB(aa1ptr,maxt,libstr,sizeof(libstr),&lmark,&icont,m_fptr,l_off_p);
  }    

  if (term_code && m_fptr->lib_aa && aa1ptr[n1-1]!=term_code) {
    aa1ptr[n1++]=term_code;
    aa1ptr[n1]=0;
  }

    /* check for subset */
    if (m_fptr->opt_text[0]!='\0') {
      if (m_fptr->opt_text[0]=='-') {
	sstart=0; sscanf(&m_fptr->opt_text[1],"%d",&sstop);
      }
      else {
	sstart = 0; sstop = -1;
	sscanf(&m_fptr->opt_text[0],"%d-%d",&sstart,&sstop);
	sstart--;
	if (sstop <= 0 ) sstop = BIGNUM;
      }

      n1 = min(n1, sstop);
      for (id=0,is=sstart; is<n1; ) {
	aa1ptr[id++]=aa1ptr[is++];
      }
      aa1ptr[id]='\0';
      n1 -= sstart;
      *l_off_p += sstart;
    }

  if (aa1ptr!=aa1) n1 += loff;

  return n1;
}
