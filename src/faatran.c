/* $Id: faatran.c $ */

/* copyright (c) 1996, 1997, 1998, 1999, 2014 by William R. Pearson and
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

/*	aatran.c	translates from nt to aa, 1 char codes */
/*	modified July 2, 1987 for all 6 frames */
/*	23 Jan 1991	fixed bug for short sequences */

/* 	this mapping is not alphabet independent */

#define XTERNAL
#include <stdio.h>
#include <stdlib.h>

#include "upam.h"
#include "uascii.h"

/*
1. The Standard Code (transl_table=1)

By default all transl_table in GenBank flatfiles are equal to id 1, and this
is not shown. When transl_table is not equal to id 1, it is shown as a
qualifier on the CDS feature.

*/
static
char *AA1="FFLLSSSSYY**CCUWLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG";
/*
  Starts = ---M---------------M---------------M----------------------------
  Base1  = TTTTTTTTTTTTTTTTCCCCCCCCCCCCCCCCAAAAAAAAAAAAAAAAGGGGGGGGGGGGGGGG
  Base2  = TTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGG
  Base3  = TCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAG

2. The Vertebrate Mitochondrial Code (transl_table=2)
*/
static
char *AA2 ="FFLLSSSSYY**CCWWLLLLPPPPHHQQRRRRIIMMTTTTNNKKSS**VVVVAAAADDEEGGGG";
/*
  Starts = --------------------------------MMMM---------------M------------
  Base1  = TTTTTTTTTTTTTTTTCCCCCCCCCCCCCCCCAAAAAAAAAAAAAAAAGGGGGGGGGGGGGGGG
  Base2  = TTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGG
  Base3  = TCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAG

3. The Yeast Mitochondrial Code (transl_table=3)
*/
static
char *AA3 ="FFLLSSSSYY**CCWWTTTTPPPPHHQQRRRRIIMMTTTTNNKKSSRRVVVVAAAADDEEGGGG";
/*
  Starts = -----------------------------------M----------------------------
  Base1  = TTTTTTTTTTTTTTTTCCCCCCCCCCCCCCCCAAAAAAAAAAAAAAAAGGGGGGGGGGGGGGGG
  Base2  = TTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGG
  Base3  = TCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAG

4. The Mold, Protozoan, and Coelenterate Mitochondrial Code and the
Mycoplasma/Spiroplasma Code (transl_table=4)
*/
static
char *AA4 ="FFLLSSSSYY**CCWWLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG";
/*
  Starts = --MM---------------M------------MMMM---------------M------------
  Base1  = TTTTTTTTTTTTTTTTCCCCCCCCCCCCCCCCAAAAAAAAAAAAAAAAGGGGGGGGGGGGGGGG
  Base2  = TTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGG
  Base3  = TCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAG

5. The Invertebrate Mitochondrial Code (transl_table=5)
*/
static
char *AA5 ="FFLLSSSSYY**CCWWLLLLPPPPHHQQRRRRIIMMTTTTNNKKSSSSVVVVAAAADDEEGGGG";
/*
  Starts = ---M----------------------------MMMM---------------M------------
  Base1  = TTTTTTTTTTTTTTTTCCCCCCCCCCCCCCCCAAAAAAAAAAAAAAAAGGGGGGGGGGGGGGGG
  Base2  = TTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGG
  Base3  = TCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAG

6. The Ciliate, Dasycladacean and Hexamita Nuclear Code (transl_table=6)
*/
static
char *AA6 ="FFLLSSSSYYQQCC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG";
/*
  Starts = -----------------------------------M----------------------------
  Base1  = TTTTTTTTTTTTTTTTCCCCCCCCCCCCCCCCAAAAAAAAAAAAAAAAGGGGGGGGGGGGGGGG
  Base2  = TTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGG
  Base3  = TCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAG

9. The Echinoderm Mitochondrial Code (transl_table=9)
*/
static
char *AA9 ="FFLLSSSSYY**CCWWLLLLPPPPHHQQRRRRIIIMTTTTNNNKSSSSVVVVAAAADDEEGGGG";
/*
  Starts = -----------------------------------M----------------------------
  Base1  = TTTTTTTTTTTTTTTTCCCCCCCCCCCCCCCCAAAAAAAAAAAAAAAAGGGGGGGGGGGGGGGG
  Base2  = TTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGG
  Base3  = TCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAG

10. The Euplotid Nuclear Code (transl_table=10)
*/
static
char *AA10="FFLLSSSSYY**CCCWLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG";
/*
  Starts = -----------------------------------M----------------------------
  Base1  = TTTTTTTTTTTTTTTTCCCCCCCCCCCCCCCCAAAAAAAAAAAAAAAAGGGGGGGGGGGGGGGG
  Base2  = TTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGG
  Base3  = TCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAG

11. The Bacterial "Code" (transl_table=11)
*/
static
char *AA11="FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG";
/*
  Starts = ---M---------------M------------MMMM---------------M------------
  Base1  = TTTTTTTTTTTTTTTTCCCCCCCCCCCCCCCCAAAAAAAAAAAAAAAAGGGGGGGGGGGGGGGG
  Base2  = TTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGG
  Base3  = TCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAG

12. The Alternative Yeast Nuclear Code (transl_table=12)
*/
static
char *AA12 ="FFLLSSSSYY**CC*WLLLSPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG";
/*
  Starts = -------------------M---------------M----------------------------
  Base1  = TTTTTTTTTTTTTTTTCCCCCCCCCCCCCCCCAAAAAAAAAAAAAAAAGGGGGGGGGGGGGGGG
  Base2  = TTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGG
  Base3  = TCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAG

13. The Ascidian Mitochondrial Code (transl_table=13)
*/
static
char *AA13="FFLLSSSSYY**CCWWLLLLPPPPHHQQRRRRIIMMTTTTNNKKSSGGVVVVAAAADDEEGGGG";
/*
  Starts = -----------------------------------M----------------------------
  Base1  = TTTTTTTTTTTTTTTTCCCCCCCCCCCCCCCCAAAAAAAAAAAAAAAAGGGGGGGGGGGGGGGG
  Base2  = TTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGG
  Base3  = TCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAG

14. The Flatworm Mitochondrial Code (transl_table=14)
*/
static
char *AA14 ="FFLLSSSSYYY*CCWWLLLLPPPPHHQQRRRRIIIMTTTTNNNKSSSSVVVVAAAADDEEGGGG";
/*
  Starts = -----------------------------------M----------------------------
  Base1  = TTTTTTTTTTTTTTTTCCCCCCCCCCCCCCCCAAAAAAAAAAAAAAAAGGGGGGGGGGGGGGGG
  Base2  = TTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGG
  Base3  = TCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAG

15. Blepharisma Nuclear Code (transl_table=15)
*/
static
char *AA15="FFLLSSSSYY*QCC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG";
/*
  Starts = -----------------------------------M----------------------------
  Base1  = TTTTTTTTTTTTTTTTCCCCCCCCCCCCCCCCAAAAAAAAAAAAAAAAGGGGGGGGGGGGGGGG
  Base2  = TTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGG
  Base3  = TCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAG
*/

static
char *AA16 ="FFLLSSSSYY*LCC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG";
/* 
  id 16 ,
  name "Chlorophycean Mitochondrial" ,
  sncbieaa "-----------------------------------M----------------------------"
  -- Base1  TTTTTTTTTTTTTTTTCCCCCCCCCCCCCCCCAAAAAAAAAAAAAAAAGGGGGGGGGGGGGGGG
  -- Base2  TTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGG
  -- Base3  TCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAG
*/

static
char *AA21 ="FFLLSSSSYY**CCWWLLLLPPPPHHQQRRRRIIMMTTTTNNNKSSSSVVVVAAAADDEEGGGG";
/*
  name "Trematode Mitochondrial" ,
  id 21 ,
  sncbieaa "-----------------------------------M---------------M------------"
  -- Base1  TTTTTTTTTTTTTTTTCCCCCCCCCCCCCCCCAAAAAAAAAAAAAAAAGGGGGGGGGGGGGGGG
  -- Base2  TTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGG
  -- Base3  TCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAG
*/

static
char *AA22 ="FFLLSS*SYY*LCC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG";
/*
  name "Scenedesmus obliquus Mitochondrial" ,
  id 22 ,
  sncbieaa "-----------------------------------M----------------------------"
  -- Base1  TTTTTTTTTTTTTTTTCCCCCCCCCCCCCCCCAAAAAAAAAAAAAAAAGGGGGGGGGGGGGGGG
  -- Base2  TTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGG
  -- Base3  TCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAG
*/

static
char *AA23 ="FF*LSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG";
/*
  name "Thraustochytrium Mitochondrial" ,
  id 23 ,
  sncbieaa "--------------------------------M--M---------------M------------"
  -- Base1  TTTTTTTTTTTTTTTTCCCCCCCCCCCCCCCCAAAAAAAAAAAAAAAAGGGGGGGGGGGGGGGG
  -- Base2  TTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGG
  -- Base3  TCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAG
*/


static char aacmap[64]={
  'K','N','K','N','T','T','T','T','R','S','R','S','I','I','M','I',
  'Q','H','Q','H','P','P','P','P','R','R','R','R','L','L','L','L',
  'E','D','E','D','A','A','A','A','G','G','G','G','V','V','V','V',
  '*','Y','*','Y','S','S','S','S','*','C','W','C','L','F','L','F'
};

static int aamap[64];	/* integer aa values */
static int aamapr[64]; /* reverse sequence map */

/* tnt is used only by aatran.c. It must be consistent with lascii and
the nt alphabet. It uses 3,3 because T and U are considered separately
*/
static int tnt[]={0,0,1,2,3,3,0,1,0,0,1,2,0,0,0,1,0,0,
		    0,1,2,3,3,0,1,0,0,1,2,0,0,0,1,0,0};

static int debug_set;

int
aatran(const unsigned char *ntseq, unsigned char *aaseq, int maxs, int frame)
{
  int iaa, im, nna, i;
  register int *nnp;
  const unsigned char *nts0;
  register int *aamp;
  register unsigned char *aap;

  iaa=nna=(maxs-(frame<3?frame:frame-3))/3;
  if (nna <= 3 ) {
    aaseq[0]=EOSEQ;
    return 0;
  }

  nnp = tnt;

  if (frame < 3) {
    aamp = aamap;
    nts0 = &ntseq[frame];
    aap = aaseq;
    while (nna--) {
      im = nnp[*nts0++]<<4;
      im += nnp[*nts0++]<<2;
      im += nnp[*nts0++];
      *aap++ = aamp[im];

      /* this check is included because of a bug in tfasty 
         which occurs only during the alignment process */

#ifdef DEBUG
      if (debug_set && aamp[im] > MAXUC) {
	fprintf(stderr,"faatran: %d %d %d %d %d?%d\n",
		*(nts0-3),*(nts0-2),*(nts0-1), im, aamp[im],aamap[im]);

	/* this allows recovery, but should not be done frequently */
	for (i=0; i<64; i++) {
	  aamap[i]=aascii[aacmap[i]];
	  aamapr[i]=aascii[aacmap[(~i)&63]];
	}
	*(aap-1) = aamp[im];
      }
#endif
    }
  }
  else {
    aamp = aamapr;
    nts0 = &ntseq[maxs-(frame-3)];
    aap = aaseq;
    while (nna--) {
      im = nnp[*--nts0]<<4;
      im += nnp[*--nts0]<<2;
      im += nnp[*--nts0];
      *aap++ = aamp[im];
      /* this check is included because of a bug in tfasty 
         which occurs only during the alignment process */

#ifdef DEBUG
      if (debug_set && aamp[im] > MAXUC) {
	fprintf(stderr,"faatran: %d %d %d %d %d?%d\n",
		*(nts0-3),*(nts0-2),*(nts0-1), im, aamp[im],aamap[im]);

	/* this allows recovery, but should not be done frequently */
	for (i=0; i<64; i++) {
	  aamap[i]=aascii[aacmap[i]];
	  aamapr[i]=aascii[aacmap[(~i)&63]];
	}
	*(aap-1) = aamp[im];
      }
#endif
    }
  }
  aaseq[iaa]=EOSEQ;
  return iaa;
}

/* slower version that masks out NNN,XXX */

/*                - A C G T U R Y M W S K D H V B N X */
static int snt[]={0,0,1,2,3,3,0,1,0,0,4,4,4,4,4,4,4,4};

int
saatran(const unsigned char *ntseq,
	unsigned char *aaseq, int maxs, int frame)
{
  int iaa, im, it, nna, xflag;
  register int *nnp;
  const unsigned char *nts0;
  register int *aamp;
  register unsigned char *aap;

  iaa=nna=(maxs-(frame<3?frame:frame-3))/3;
  if (nna <= 3 ) {
    aaseq[0]=EOSEQ;
    return 0;
  }

  nnp = snt;
  if (frame < 3) {
    aamp = aamap;
    nts0 = &ntseq[frame];
    aap = aaseq;
    while (nna--) {
      xflag = 0;
      if ((it=nnp[*nts0++])<4) {im = it<<4;}
      else {xflag = 1; im=0;}
      if ((it=nnp[*nts0++])<4) {im += it<<2;}
      else xflag = 1;
      if ((it=nnp[*nts0++])<4) {im += it;}
      else xflag = 1;
      if (xflag) *aap++ = aascii['X'];
      else *aap++ = aamp[im];
    }
  }
  else {
    aamp = aamapr;
    nts0 = &ntseq[maxs-(frame-3)];
    aap = aaseq;
    while (nna--) {
      xflag = 0;
      if ((it=nnp[*--nts0]) < 4) im = it<<4;
      else {xflag = 1; im=0;}
      if ((it=nnp[*--nts0]) < 4) im += it<<2;
      else xflag = 1;
      if ((it=nnp[*--nts0]) < 4) im += it;
      else xflag = 1;
      if (xflag) *aap++ = aascii['X'];
      else *aap++ = aamp[im];
    }
  }
  aaseq[iaa]=EOSEQ;
  return iaa;
}

void
aainit(int tr_type, int debug)
{
  int i,j;
  char *aasmap;
  int ascii_star;
  int imap[4]={3,1,0,2}, i0, i1, i2, ii;

  debug_set = debug;

  aasmap = AA1;

  ascii_star = aascii['*'];
  aascii['*'] = TERM;

  if (tr_type > 0) {
    /* need to put in a new translation table */
    switch (tr_type) {
    case 1: aasmap = AA1; break;
    case 2: aasmap = AA2; break;
    case 3: aasmap = AA3; break;
    case 4: aasmap = AA4; break;
    case 5: aasmap = AA5; break;
    case 6: aasmap = AA6; break;
    case 9: aasmap = AA9; break;
    case 10: aasmap = AA10; break;
    case 11: aasmap = AA11; break;
    case 12: aasmap = AA12; break;
    case 13: aasmap = AA13; break;
    case 14: aasmap = AA14; break;
    case 15: aasmap = AA15; break;
    case 16: aasmap = AA16; break;
    case 21: aasmap = AA21; break;
    case 22: aasmap = AA22; break;
    case 23: aasmap = AA23; break;

    default: aasmap = AA1; break;
    }

    if (debug) fprintf(stderr," codon table: %d\n     new old\n",tr_type);
    for (i0 = 0; i0 < 4; i0++)
      for (i1 = 0; i1 < 4; i1++)
	for (i2 = 0; i2 < 4; i2++) {
	  ii = (imap[i0]<<4) + (imap[i1]<<2) + imap[i2];
	  if (debug &&  aacmap[ii] != *aasmap) {
	    fprintf(stderr," %c%c%c: %c - %c\n",
		    nt[imap[i0]+1],nt[imap[i1]+1],nt[imap[i2]+1],
		    *aasmap,aacmap[ii]);
	  }
	  aacmap[ii]= *aasmap++;
	}

    if (debug) {
      for (i=0; i<64; i++) {
	fprintf(stderr,"'%c',",aacmap[i]);
	if ((i%16)==15) fputc('\n',stderr);
      }
      fputc('\n',stderr);
    }
  }
  for (i=0; i<64; i++) {
    aamap[i]=aascii[aacmap[i]];
    if (aamap[i] > TERM) {
      fprintf(stderr," *** ERROR - codon out of range: %d %d (%c)\n",i,aamap[i], NCBIstdaa_l[aamap[i]] );
    }
    aamapr[i]=aascii[aacmap[(~i)&63]];
    if (aamapr[i] > TERM) {
      fprintf(stderr," *** ERROR - codon_r out of range: %d %d (%c)\n",i,aamapr[i], NCBIstdaa_l[aamapr[i]]);
    }
  }
  aascii['*'] = ascii_star;
}

void
aagetmap(char *to, int n) 
{
  int i;
  for (i=0; i<n; i++) to[i] = aacmap[i];
}
