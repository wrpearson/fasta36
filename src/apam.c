/* apam.c	19-June-86 */

/* $Id: apam.c 1281 2014-08-21 17:32:06Z wrp $ */

/* copyright (c) 1987, 2014 by William R. Pearson and The Rector &
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

/*
	read in the alphabet and pam matrix data
	designed for universal matcher

	This version reads BLAST format (square) PAM files
*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>

#include "defs.h"
#include "param.h"

#define XTERNAL
#include "uascii.h"
#include "upam.h"
#undef XTERNAL

extern void alloc_pam (int d1, int d2, struct pstruct *ppst);
extern void init_altpam(struct pstruct *ppst);

/* pam_opts -- modify PAM matrix (pamoff, -MS) if -MS or +off is part
   of PAM matrix name, e.g. MD20-MS or BL50+2 
*/

void
pam_opts(char *smstr, struct pstruct *ppst) {
  char *bp;

  ppst->pam_ms = 0;
  ppst->pamoff = 0;

  if ((bp=strchr(smstr,'-'))!=NULL) {
    if (!strncmp(bp+1,"MS",2) || !strncmp(bp+1,"ms",2)) {
      ppst->pam_ms = 1;
    }
    else {
      ppst->pamoff=atoi(bp+1);
    }
    *bp = '\0';
   }
  else if ((bp=strchr(smstr,'+'))!=NULL) {
    ppst->pamoff= -atoi(bp+1);
    *bp = '\0';
  }
}

/* modified 13-Oct-2005 to accomodate asymmetrical matrices */
/* modified 15-Jul-2010 to ensure constant NCBIstdaa encoding */
/* ensure that all entries in NCBIstdaa have values */

int
initpam (char *mfname, struct pstruct *ppst)
{
   char    line[512], *lp;
   int     i, j, iaa, pval, p_i, p_j;
   int l_nsq;
   unsigned char l_sq[MAXSQ+1];
   int ess_tmp, max_val, min_val;
   int have_es = 0;
   FILE   *fmat;

   pam_opts(mfname, ppst);

   if ((fmat = fopen (mfname, "r")) == NULL)
   {
      printf ("***WARNING*** cannot open scoring matrix file %s\n", mfname);
      fprintf (stderr,"***WARNING*** cannot open scoring matrix file %s\n", mfname);
      return 0;
   }

/* removed because redundant, and causes crash under MacOSX -- because copying on top of itself */
/*
   SAFE_STRNCPY (ppst->pamfile, mfname, MAX_FN);
*/
   SAFE_STRNCPY(ppst->pam_name, ppst->pamfile, MAX_FN);

   if (ppst->pam_ms) {
     SAFE_STRNCAT(ppst->pam_name,"-MS",MAX_FN-strlen(ppst->pam_name));
   }

   /* 
      the size of the alphabet is determined in advance 
   */
   ppst->nt_align = (ppst->dnaseq == SEQT_DNA || ppst->dnaseq == SEQT_RNA);

   /* 
     look for alphabet line, skipping the comments, alphabet ends up in line[]
   */
   while (fgets (line, sizeof(line), fmat) != NULL && line[0]=='#');

   /* transfer the residue line into l_sq[] */
   l_nsq = 1;
   l_sq[0] = '\0';
   for (i=0; i<strlen(line); i++) {
     if (isalpha(line[i]) || line[i] == '*') {
       l_sq[l_nsq++] = line[i];
     }
   }

   /* if we have a DNA matrix, various defaults must be updated,
      particularly pascii, which is used to map the residue ordering
      in the matrix file to the residue ordering used by the
      program */

   if (l_nsq < 20) {
     if (ppst->dnaseq <= SEQT_PROT) {
       ppst->dnaseq = SEQT_DNA;
     }
     ppst->nt_align=1;
     pascii = nascii;	/* use correct DNA mapping, NCBIstdaa by default */
   }

   /* we no-longer re-initialize sascii[], we either use NCBIstdaa
      mapping for protein, or nascii for DNA */

   /* 11-July-2014 -- need to check that alphabet is consistent with pascii */
   /* 
   for (i=0; i < l_nsq; i++) {
   }
   */

   /* check for 2D pam  - if not found, allocate it */
   if (!ppst->have_pam2) {
     alloc_pam (MAXSQ+1, MAXSQ+1, ppst);
     ppst->have_pam2 = 1;
   }

   max_val = -1;
   min_val =  1;
   ppst->pam2[0][0][0] = -BIGNUM;
   /* make certain the [0] boundaries are -BIGNUM */
   for (j=1; j < l_nsq; j++) {
     p_j = pascii[l_sq[j]];
     ppst->pam2[0][0][p_j] = ppst->pam2[0][p_j][0] = -BIGNUM;
   }

   /*  read the scoring matrix values */
   for (iaa = 1; iaa < l_nsq; iaa++) {	/* read pam value line */
     p_i = pascii[l_sq[iaa]];
     if (p_i > MAXSQ) {
       fprintf(stderr,"*** error [%s:%d] - residue character %c out of range %d\n",
	       __FILE__, __LINE__, l_sq[iaa], p_i);
       p_i = pascii['X'];
     }
     if (fgets(line,sizeof(line),fmat)==NULL) {
       fprintf (stderr," error reading pam line: %s\n",line);
       exit (1);
     }
     /*     fprintf(stderr,"%d/%d %s",iaa,nsq,line); */
     strtok(line," \t\n");		/* skip the letter (residue) */

     for (j = 1; j < l_nsq; j++) {
       p_j = pascii[l_sq[j]];
       lp=strtok(NULL," \t\n");		/* get the number string */
       pval=ppst->pam2[0][p_i][p_j]=atoi(lp);	/* convert to integer */
       if (pval > max_val) max_val = pval;
       if (pval < min_val) min_val = pval;
     }
   }
   ppst->pam_h = max_val;
   ppst->pam_l = min_val;

   if (ppst->dnaseq==0) {
     pam_sq = apam_sq;
     pam_sq_n = apam_sq_n;
     init_altpam(ppst);
   }
   else {
     pam_sq = npam_sq;
     pam_sq_n = npam_sq_n;
   }

   /* is protein but do not have '*' in alphabet*/
   p_i = pascii['*'];
   p_j = pascii['X'];
   if (!ppst->nt_align && strchr((char *)l_sq,'*')==NULL) {
      /* add it */
     for (i=0; i< l_nsq; i++) {
       ppst->pam2[0][p_i][i] = ppst->pam2[0][p_j][i];
       ppst->pam2[0][i][p_i] = ppst->pam2[0][i][p_j];
     }
   }

   /* make sure that X:X is < 0 if -S */
   if (ppst->ext_sq_set && ppst->pam2[0][p_j][p_j] >= 0) {
     ppst->pam2[0][p_j][p_j] = -1;
   }

   fclose (fmat);
   return 1;
}

/* make a DNA scoring from +match/-mismatch values */

void mk_n_pam(int *arr,int siz, int mat, int mis)
{
  int i, j, k;
  /* current default match/mismatch values */
  int max_mat = +5;
  int min_mis = -4;
  float f_val, f_scale;
  
  f_scale = (float)(mat - mis)/(float)(max_mat - min_mis);
  
  k = 0;
  for (i = 0; i<nnt-1; i++)
    for (j = 0; j <= i; j++ ) {
      if (arr[k] == max_mat) arr[k] = mat;
      else if (arr[k] == min_mis) arr[k] = mis;
      else if (arr[k] != -1) { 
	f_val = (arr[k] - min_mis)*f_scale + 0.5;
	arr[k] = f_val + mis;
      }
      k++;
    }
}

int
standard_pam(char *smstr, struct pstruct *ppst, int del_set, int gap_set) {

  struct std_pam_str *std_pam_p;

  pam_opts(smstr, ppst);

  for (std_pam_p = std_pams; std_pam_p->abbrev[0]; std_pam_p++ ) {
    if (strcmp(smstr,std_pam_p->abbrev)==0) {
      pam = std_pam_p->pam;
      strncpy(ppst->pam_name,std_pam_p->name,MAX_FN);
      ppst->pam_name[MAX_FN-1]='\0';
      if (ppst->pam_ms) {
	strncat(ppst->pam_name,"-MS",MAX_FN-strlen(ppst->pam_name)-1);
      }
      ppst->pam_name[MAX_FN-1]='\0';
#ifdef OLD_FASTA_GAP
      if (!del_set) ppst->gdelval = std_pam_p->gdel+std_pam_p->ggap;
#else
      if (!del_set) ppst->gdelval = std_pam_p->gdel;
#endif
      if (!gap_set) ppst->ggapval = std_pam_p->ggap;
      ppst->pamscale = std_pam_p->scale;
      return 1;
    }
  }
  return 0;
}

/* scan through the sorted (by decreasing entropy) std_pams[] array
   for a scoring matrix with enough entropy
*/
int
min_pam_bits(int n0_eff, double bit_thresh, struct pstruct *ppst, int del_set, int gap_set) {
  struct std_pam_str *std_pam_p;
  int curr_pam_idx = 0;

  pam_opts(ppst->pamfile, ppst);

  /* get the index for the current (standard) pam file */
  for (curr_pam_idx = 0; std_pams[curr_pam_idx].abbrev[0]; curr_pam_idx++) {
    if (strcmp(ppst->pamfile,std_pams[curr_pam_idx].abbrev)==0) break;
  }

  /* only use matrices from the VT series */
  for ( ; curr_pam_idx > 0 ; curr_pam_idx-- ) {
    if ((strncmp(std_pams[curr_pam_idx].name,"VT",2)!=0) &&
	(strcmp(ppst->pamfile, std_pams[curr_pam_idx].abbrev) != 0)) continue;
    if (n0_eff * std_pams[curr_pam_idx].entropy >= bit_thresh) goto new_pam;
  }
  return 0;

 new_pam:
  std_pam_p = &std_pams[curr_pam_idx];

  pam = std_pam_p->pam;
  strncpy(ppst->pam_name,std_pam_p->name,MAX_FN);
  if (ppst->pam_ms) {
    strncat(ppst->pam_name,"-MS",MAX_FN-strlen(ppst->pamfile)-1);
  }
  ppst->pam_name[MAX_FN-1]='\0';
#ifdef OLD_FASTA_GAP
  if (!del_set) ppst->gdelval = std_pam_p->gdel+std_pam_p->ggap;
#else
  if (!del_set) ppst->gdelval = std_pam_p->gdel;
#endif
  if (!gap_set) ppst->ggapval = std_pam_p->ggap;
  ppst->pamscale = std_pam_p->scale;
  return 1;
}

/* build_xascii is only used for SEQT_UNK - it replaces the default
   input mapping (aascii[])  with a mapping that preserves any letter
   in either the aax[], ntx[], or othx[] alphabets, or in
   save_str[]. othx[] was added to support letters that are mapped,
   but are not (yet) in aax[], e.g. 'OoUu'.  Because build_xascii
   makes a qascii[] that is all ascii characters with values > '@',
   these values must be replaced using either aascii[] or nascii[] and
   the initial query sequence re-coded.
 */
void
build_xascii(int *qascii, char *save_str) {
  int i, max_save;
  int comma_val, term_val;
  int save_arr[MAX_SSTR];

  comma_val = qascii[','];
  term_val = qascii['*'];

  /* preserve special characters */
  for (i=0; i < MAX_SSTR && save_str[i]; i++ ) {
    save_arr[i] = qascii[save_str[i]];
  }
  max_save = i;

  for (i=1; i<128; i++) {
    qascii[i]=NA;
  }

  /* range of values in aax, ntx is from 1..naax,nntx - 
     do not zero-out qascii[0] - 9 Oct 2002 */

  for (i=1; i<NCBIstdaa_ext_n; i++) {
    qascii[NCBIstdaa_ext[i]]=NCBIstdaa_ext[i];
  }

  for (i=1; i<nntx; i++) {
    qascii[ntx[i]]=ntx[i];
  }

  /* put it letters that are not in other alphabets because they are
     mapped -- now included in NCBIstdaa */
  /*
  for (i=1; i<nothx; i++) {
    qascii[othx[i]]=othx[i];
  }
  */

  qascii['\n'] = qascii['\r'] = EL;

  qascii[','] = comma_val;
  qascii['*'] = term_val;
  qascii[0] = ES;

  for (i=0; i < max_save; i++) {
    qascii[save_str[i]]=save_arr[i];
  }
}


/* init_ascii0 -- initializes an ascii mapping from a sequence
ordering
*/
void 
init_ascii0(int *xascii, char *sq_map, int n_sq_map, struct pstruct *ppst) {
  int i;

  /* first map everything as non-sequence */
  for (i=0; i<128; i++) {
    xascii[i] = NA;
  }

  /* then map the actual sequence letters */
  for (i = 1; i < n_sq_map; i++) {
    xascii[sq_map[i]] = i;
    if (n_sq_map <= MAXUC) { /* only uppercase */
      xascii[tolower(sq_map[i])] = i;	/* map lowercase */
    }
  }

  ppst->nsq = n_sq_map;
  for (i=1; i < n_sq_map; i++) {
    ppst->sq[i] = sq_map[i];
  }
  ppst->sq[0] = 0;

  /* then map the other stuff, EL  etc */
  xascii[0] = ES;
  xascii[10] = EL;
  xascii[13] = EL;
}

/* init_ascii()

   checks for lower case letters in *sq array;
   if not present, map lowercase to upper

*/
void
init_ascii(int is_ext, int *xascii, int p_nsq, int is_dna) {

  int isq, have_lc;
  char *sq, term_char;
  int nsq;
  
  if (is_dna==SEQT_UNK) return;

  term_char = xascii['*'];

  if (is_dna==SEQT_DNA || is_dna == SEQT_RNA) {
    if (is_ext) { 
      sq = &ntx[0];
      nsq = nntx;
    }
    else {sq = &nt[0]; nsq = nnt;}
  }
  else {
    if (is_ext) { sq = NCBIstdaa_ext; nsq = NCBIstdaa_ext_n; }
    else {sq = NCBIstdaa; nsq = NCBIstdaa_n;}
  }

  /* initialize xascii from sq[], checking for lower-case letters */
  /* this code guarantees that all characters in sq are represented in
     xascii[], but it does not guarantee that everything else in
     xascii[], particularly xascii[O,U,J], have appropriate values
 */
  have_lc = 0;
  for (isq = 1; isq <= nsq; isq++) {
     xascii[sq[isq]] = isq;
     if (sq[isq] >= 'a' && sq[isq] <= 'z') have_lc = 1;
  }
  
  /* no lower case letters in alphabet, map lower case to upper */
  if (have_lc != 1) { 
    for (isq = 1; isq <= nsq; isq++) {
      if (sq[isq] >= 'A' && sq[isq] <= 'Z') xascii[sq[isq]-'A'+'a'] = isq;
    }
    if (is_dna==1) xascii['u'] = xascii['t'];
  }

  xascii['*']=term_char;
  xascii[0] = ES;
}

void
validate_novel_aa(int *xascii, int p_nsq, int dnaseq) {
  int isq, err_p_nsq_limit;
  /* these checks need to be done after xascii[] has been
     re-initialized */

  if (dnaseq != SEQT_DNA && dnaseq!=SEQT_RNA) {
    if (xascii['O'] > p_nsq || xascii['o'] > p_nsq) { xascii['O'] = xascii['K']; xascii['o'] = xascii['k'];}
    if (xascii['U'] > p_nsq || xascii['u'] > p_nsq) { xascii['U'] = xascii['C']; xascii['u'] = xascii['c'];}
    if (xascii['J'] > p_nsq || xascii['j'] > p_nsq) { xascii['J'] = xascii['L']; xascii['j'] = xascii['l'];}
  }

  /* one final check for characters out of range (> nsq)*/
  err_p_nsq_limit = 0;
  for (isq = 'A'; isq <= 'Z'; isq++) {
    if (xascii[isq] < NA && xascii[isq] > p_nsq) {
      fprintf(stderr, " *** ERROR *** xascii['%c']:%d > %d\n",isq, xascii[isq], p_nsq);
      err_p_nsq_limit = 1;
    }
  }
  if (err_p_nsq_limit) {exit(1);}
}


void 
print_pam(struct pstruct *ppst) {
  int i, nsq, ip;
  unsigned char *sq;

  fprintf(stderr," ext_sq_set: %d\n",ppst->ext_sq_set);

  nsq = ppst->nsq;
  ip = 0;
  sq = ppst->sq;

  fprintf(stderr," sq[%d]: %s\n",nsq, sq);

  if (ppst->ext_sq_set) {
    nsq = ppst->nsqx;
    ip = 1;
    sq = ppst->sqx;
    fprintf(stderr," sq[%d]: %s\n",nsq, sq);
  }

  for (i=1; i<=nsq; i++) {
    fprintf(stderr," %c:%c - %3d\n",sq[i], sq[i], ppst->pam2[ip][i][i]);
  }
}
