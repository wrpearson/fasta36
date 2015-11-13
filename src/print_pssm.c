/* print_pssm.c - 21-Jan-2005  */

/* $Id: print_pssm.c 1111 2013-01-09 18:46:57Z wrp $  */
/* $Revision: 1111 $  */

/* copyright (c) 2005, 2014 - William R. Pearson and The Rector &
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
   read a binary PSSM checkpoint file from blastpgp, and produce an ascii
   formatted file
*/

#include <stdio.h>
#include <stdlib.h>
#include <sys/types.h>
#include <math.h>
#include <string.h>
#ifdef UNIX
#include <getopt.h>
#else
extern int optind;		/* used by getopt() */
extern char *optarg;
#endif

#include "defs.h"
#include "param.h"

#include "uascii.h"
#include "upam.h"

void initenv(int, char **, struct pstruct *, char *);
void read_pssm(unsigned char *aa0, int n0, int nsq, double pamscale, 
	       FILE *fp, int pgpf_type, struct pstruct *ppst);
int
read_asn_pssm(unsigned char *aa0, int n0, int nsq,
	      double pamscale, FILE *fp, struct pstruct *ppst);
void alloc_pam();
int **alloc_pam2p();
void initpam2();
void init_ascii0(int *xascii, char *sq_map, int n_sq_map, struct pstruct *ppst);
void fill_pam();
double get_lambda(int **pam2p, int n0, int nsq, char *aa0);

int standard_pam(char *smstr, struct pstruct *ppst, int del_set, int gap_set);
int getseq(char *filen, int *qascii, unsigned char *seq, int maxs, char *libstr,
	   int n_libstr, long *sq0off);
int initpam (char *mfname, struct pstruct *ppst);
int karlin(int , int, double *, double *, double *);

int
main(int argc, char **argv) {

  unsigned char *aa0;
  char libstr[MAX_FN];
  char qname[MAX_FN];
  long sq0off;
  int i, n0;
  FILE *fp;
  struct pstruct pst, *ppst;

  /* stuff from initfa.c/h_init() */

  memcpy(qascii,aascii,sizeof(qascii));

  /* initialize a pam matrix */
  ppst = &pst;
  strncpy(ppst->pamfile,"BP62",MAX_FN);
  standard_pam(ppst->pamfile,ppst,0,0);

  /* this is always protein by default */
  ppst->nsq = naa;
  ppst->nsqx = naax;
  for (i=0; i<=ppst->nsqx; i++) {
    ppst->sq[i] = NCBIstdaa[i];
  }
  ppst->sq[ppst->nsqx+1] = ppst->sqx[ppst->nsqx+1] = '\0';

  if ((aa0 = calloc(MAXTST,sizeof(char)))==NULL) {
    fprintf(stderr,"Cannot allocate aa0\n");
    exit(1);
  }

  initenv(argc, argv, &pst, qname);

  if (argc < 1) {
    fprintf(stderr,"usage -- print_pssm -P \"pssm.asn 2\" query_file\n");
    exit(1);
  }

  alloc_pam(pst.nsq+1,pst.nsq+1, &pst);
  initpam2(&pst);
  init_ascii0(qascii, NCBIstdaa_ext, NCBIstdaa_ext_n, &pst);

  n0 = getseq (qname, qascii, aa0, MAXTST, libstr, sizeof(libstr), &sq0off);

  if (!pst.pam_pssm) {
    fprintf(stderr," ** ERROR ** No -P PSSM provided\n");
  }
  else {
    pst.pam2p[0] = alloc_pam2p(n0,pst.nsq);
    pst.pam2p[1] = alloc_pam2p(n0,pst.nsq);

    if ((pst.pgpfile_type == 0) && (fp=fopen(pst.pgpfile,"rb"))) {
      read_pssm(aa0, n0, pst.nsq, pst.pamscale, fp, 0, &pst);
    }
    else if ((pst.pgpfile_type == 1) && (fp=fopen(pst.pgpfile,"r"))) {
      read_pssm(aa0, n0, pst.nsq, pst.pamscale, fp, 1, &pst);
    }
    else if ((pst.pgpfile_type == 2) && (fp=fopen(pst.pgpfile,"rb"))) {
      if (read_asn_pssm(aa0, n0, pst.nsq, pst.pamscale, fp, &pst)<=0)
	fprintf(stderr," Could not parse PSSM file: %s\n",pst.pgpfile);
	pst.pam_pssm = 0;
	exit(1);
    }
  }
}

void
initenv(int argc, char **argv, struct pstruct *ppst, char *qname) {
  char copt;
  char *bp;

  pascii = aascii;

  while ((copt = getopt(argc, argv, "P:s:"))!=EOF) {
    switch (copt) {
      case 'P':
	if ((bp=(strchr(optarg,' ')))!=NULL) {
	  *bp = '\0';
	  ppst->pgpfile_type = atoi(bp+1);
	}
	strncpy(ppst->pgpfile,optarg,MAX_FN);
	ppst->pgpfile[MAX_FN-1]='\0';
	ppst->pam_pssm = 1;
	break;

      case 's':
	strncpy (ppst->pamfile, optarg, 120);
	ppst->pamfile[120-1]='\0';
	if (!standard_pam(ppst->pamfile,ppst,0, 0)) {
	  initpam (ppst->pamfile, ppst);
	}
	ppst->pam_set=1;
	break;
    }
  }
  optind--;

  if (argc - optind > 1) strncpy(qname, argv[optind+1], MAX_FN);
}


/*
   *aa0 - query sequence
   n0   - length
   pamscale - scaling for pam matrix - provided by apam.c, either
              0.346574 = ln(2)/2 (P120, BL62) or
	      0.231049 = ln(2)/3 (P250, BL50) 
*/

#define N_EFFECT 20

void
read_pssm(unsigned char *aa0, int n0, int nsq, double pamscale, FILE *fp, int pgf_type, struct pstruct *ppst) {
  int i, j, len;
  int qi, rj;
  int **pam2p;
  int first, too_high;
  char *query;
  double freq, **freq2d, lambda, new_lambda;
  double scale, scale_high, scale_low;

  pam2p = ppst->pam2p[0];

  if(1 != fread(&len, sizeof(int), 1, fp)) {
    fprintf(stderr, "error reading from checkpoint file: %d\n", len);
    exit(1);
  }

  if(len != n0) {
    fprintf(stderr, "profile length (%d) and query length (%d) don't match!\n",
	    len,n0);
    exit(1);
  }

  /* read over query sequence stored in BLAST profile */
  if(NULL == (query = (char *) calloc(len, sizeof(char)))) {
    fprintf(stderr, "Couldn't allocate memory for query!\n");
    exit(1);
  }

  if(len != fread(query, sizeof(char), len, fp)) {
    fprintf(stderr, "Couldn't read query sequence from profile: %s\n", query);
    exit(1);
  }

  printf("%d\n%s\n",len,query);

  /* currently we don't do anything with query; ideally, we should
     check to see that it actually matches aa0 ... */

  /* quick 2d array alloc: */
  if((freq2d = (double **) calloc(n0, sizeof(double *))) == NULL) {
    fprintf(stderr, "Couldn't allocate memory for frequencies!\n");
    exit(1);
  }

  if((freq2d[0] = (double *) calloc(n0 * N_EFFECT, sizeof(double))) == NULL) {
    fprintf(stderr, "Couldn't allocate memory for frequencies!\n");
    exit(1);
  }

  /* a little pointer arithmetic to fill out 2d array: */
  for (qi = 1 ; qi < n0 ; qi++) {
    freq2d[qi] = freq2d[0] + (N_EFFECT * qi);
  }

  printf(" rrtotal: %ld; pamscale=%0.5f\n", rrtotal, pamscale);

  printf("   ");
  for (rj = 1; rj <= N_EFFECT; rj++ ) {
    printf(" %c    ",NCBIstdaa[rj]);
  }
  printf("\n");

  for (qi = 0 ; qi < n0 ; qi++) {
    printf("%c",query[qi]);
    for (rj = 0 ; rj < N_EFFECT ; rj++) {
      if(1 != fread(&freq, sizeof(double), 1, fp)) {
	fprintf(stderr, "Error while reading frequencies!\n");
	exit(1);
      }
      /* printf(" %8.7g",freq*10.0); */

      if (freq > 1e-12) {
	freq = log(freq /((double) (rrcounts[rj+1])/(double) rrtotal));
	freq /= pamscale; /* this gets us close to originial pam scores */
	freq2d[qi][rj] = freq;
      }
      else {freq2d[qi][rj] = freq;}
      printf(" %5.2f",freq);
    }
    printf("\n");
  }


  /* now figure out the right scale */
  scale = 1.0;
  lambda = get_lambda(ppst->pam2[0], 20, 20, "\0ARNDCQEGHILKMFPSTWYV");

  /* should be near 1.0 because of our initial scaling by ppst->pamscale */
  fprintf(stderr, "real_lambda: %g\n", lambda);

  /* get initial high/low scale values: */
  first = 1;
  while (1) {
    fill_pam(pam2p, n0, 20, freq2d, scale);
    new_lambda = get_lambda(pam2p, n0, 20, query); 
    fprintf(stderr, "new_lambda: %g; scale: %g\n",new_lambda,scale);

    if (new_lambda > lambda) {
      if (first) {
	first = 0;
	scale = scale_high = 1.0 + 0.05;
	scale_low = 1.0;
	too_high = 1;
      } else {
	if (!too_high) break;
	scale = (scale_high += scale_high - 1.0);
      }
    } else if (new_lambda > 0) {
      if (first) {
	first = 0;
	scale_high = 1.0;
	scale = scale_low = 1.0 - 0.05;
	too_high = 0;
      } else {
	if (too_high) break;
	scale = (scale_low += scale_low - 1.0);
      }
    } else {  /* new_lambda <= 0 */
      fprintf(stderr, "new_lambda (%g) <= 0; matrix has positive average score", new_lambda);
      exit(1);
    }
  }

  /* now do binary search between low and high */
  for (i = 0 ; i < 10 ; i++) {
    scale = 0.5 * (scale_high + scale_low);
    fill_pam(pam2p, n0, 20, freq2d, scale);
    new_lambda = get_lambda(pam2p, n0, 20, query);
    fprintf(stderr, "it: %d - new_lambda: %g; scale: %g\n",i,new_lambda,scale);
    
    if (new_lambda > lambda) scale_low = scale;
    else scale_high = scale;
  }

  scale = 0.5 * (scale_high + scale_low);
  fill_pam(pam2p, n0, 20, freq2d, scale);

  fprintf(stderr, "final scale: %g\n", scale);

  fprintf(stderr,"         ");
  for (rj = 1; rj <= N_EFFECT+3; rj++ ) {
    fprintf(stderr,"   %c",NCBIstdaa[rj]);
  }
  fprintf(stderr,"\n");

  for (qi = 0 ; qi < n0 ; qi++) {
    fprintf(stderr, "%4d %c:  ", qi+1, query[qi]);
    for (rj = 1 ; rj <= N_EFFECT+3 ; rj++) {
      fprintf(stderr, "%4d", pam2p[qi][rj]);
    }
    fprintf(stderr, "\n");
  }

  free(freq2d[0]);
  free(freq2d);

  free(query);
}

/*
 * alloc_pam(): allocates memory for the 2D pam matrix as well
 * as for the integer array used to transmit the pam matrix
 */
void
alloc_pam (int d1, int d2, struct pstruct *ppst)
{
  int     i, *d2p;


  if ((ppst->pam2[0] = (int **) malloc (d1 * sizeof (int *))) == NULL) {
     fprintf(stderr,"Cannot allocate 2D pam matrix: %d",d1);
     exit(1);
  }

  if ((ppst->pam2[1] = (int **) malloc (d1 * sizeof (int *))) == NULL) {
     fprintf(stderr,"Cannot allocate 2D pam matrix: %d",d1);
     exit(1);
  }

  if ((d2p = (int *) malloc (d1 * d2 * sizeof (int))) == NULL) {
     fprintf(stderr,"Cannot allocate 2D pam matrix: %d",d1);
     exit(1);
   }

   for (i = 0; i < d1; i++, d2p += d2)
      ppst->pam2[0][i] = d2p;

   if ((d2p= (int *) malloc (d1 * d2 * sizeof (int))) == NULL) {
     fprintf(stderr,"Cannot allocate 2d pam matrix: %d",d2);
     exit(1);
   }

   for (i = 0;  i < d1; i++, d2p += d2)
      ppst->pam2[1][i] = d2p;
}

void
fill_pam(int **pam2p, int n0, int nsq, double **freq2d, double scale) {
  int i, j, n_j;
  double freq;

  /* fprintf(stderr, "scale: %g\n", scale); */
  
  /* now fill in the pam matrix: */
  for (j = 1 ; j <=nsq ; j++) {
    n_j = qascii[pssm_aa[j]];
    for (i = 0 ; i < n0 ; i++) {
      freq = scale * freq2d[i][j-1];
      if ( freq < 0.0) freq -= 0.5;
      else freq += 0.5;
      pam2p[i][n_j] = (int)freq;
    }
  }
}

/*
 *  initpam2(struct pstruct pst): Converts 1-D pam matrix to 2-D
 */
void initpam2 (struct pstruct *ppst)
{
   int i, j, k, nsq, pam_xx, pam_xm;
   int sa_x, sa_t, tmp;

   nsq = ppst->nsq;
   sa_x = pascii['X'];
   sa_t = pascii['*'];

   ppst->pam2[0][0][0] = -BIGNUM;
   ppst->pam_h = -1; ppst->pam_l = 1;

   k = 0;
   for (i = 1; i <= nsq; i++) {
     ppst->pam2[0][0][i] = ppst->pam2[0][i][0] = -BIGNUM;
     for (j = 1; j <= i; j++) {
       ppst->pam2[0][j][i] = ppst->pam2[0][i][j] = pam[k++] - ppst->pamoff;
       if (ppst->pam_l > ppst->pam2[0][i][j]) ppst->pam_l =ppst->pam2[0][i][j];
       if (ppst->pam_h < ppst->pam2[0][i][j]) ppst->pam_h =ppst->pam2[0][i][j];
     }
   }

   ppst->nt_align = (ppst->dnaseq== SEQT_DNA || ppst->dnaseq == SEQT_RNA);

   if (ppst->dnaseq == SEQT_RNA) {
     tmp = ppst->pam2[0][nascii['G']][nascii['G']] - 1;
     ppst->pam2[0][nascii['A']][nascii['G']] = 
       ppst->pam2[0][nascii['C']][nascii['T']] = 
       ppst->pam2[0][nascii['C']][nascii['U']] = tmp;
   }

   if (ppst->pam_x_set) {
     for (i=1; i<=nsq; i++) {
       ppst->pam2[0][sa_x][i] = ppst->pam2[0][i][sa_x]=ppst->pam_xm;
       ppst->pam2[0][sa_t][i] = ppst->pam2[0][i][sa_t]=ppst->pam_xm;
     }
     ppst->pam2[0][sa_x][sa_x]=ppst->pam_xx;
     ppst->pam2[0][sa_t][sa_t]=ppst->pam_xm;
   }
   else {
     ppst->pam_xx = ppst->pam2[0][sa_x][sa_x];
     ppst->pam_xm = ppst->pam2[0][1][sa_x];
   }
}

double
get_lambda(int **pam2p, int n0, int nsq, char *aa0) {
  double lambda, H;
  double *pr, tot, sum;
  int aa0i;
  int i, ioff, j, min, max;

  /* get min and max scores */
  min = BIGNUM;
  max = -BIGNUM;
  if(pam2p[0][1] == -BIGNUM) {
    ioff = 1;
    n0++;
  } else {
    ioff = 0;
  }

  for (i = ioff ; i < n0 ; i++) {
    for (j = 1; j <= nsq ; j++) {
      if (min > pam2p[i][j])
	min = pam2p[i][j];
      if (max < pam2p[i][j])
	max = pam2p[i][j];
    }
  }

  fprintf(stderr, "min: %d\tmax:%d\n", min, max);
  
  if ((pr = (double *) calloc(max - min + 1, sizeof(double))) == NULL) {
    fprintf(stderr, "Couldn't allocate memory for score probabilities: %d\n", max - min + 1);
    exit(1);
  }

  tot = (double) rrtotal * (double) rrtotal * (double) n0;
  for (i = ioff ; i < n0 ; i++) {
    if (aa0[i] < 'A') {aa0i = aa0[i];}
    else {aa0i = aascii[aa0[i]];}
    for (j = 1; j <= nsq ; j++) {
      pr[pam2p[i][j] - min] +=
	(double) ((double) rrcounts[aa0i] * (double) rrcounts[j]) / tot;
    }
  }

  sum = 0.0;
  for(i = 0 ; i <= max-min ; i++) { 
    sum += pr[i];
    /*    fprintf(stderr, "%3d: %g %g\n", i+min, pr[i], sum); */
  }
  /* fprintf(stderr, "pr[] sum: %g\n", sum); */

  for(i = 0 ; i <= max-min ; i++) { pr[i] /= sum; }

  if (!karlin(min, max, pr, &lambda, &H)) {
    fprintf(stderr, "Karlin lambda estimation failed\n");
  }

  /*   fprintf(stderr, "lambda: %g\n", lambda); */
  free(pr);

  return lambda;
}

int **
alloc_pam2p(int len, int nsq) {
  int i;
  int **pam2p;

  if ((pam2p = (int **)calloc(len,sizeof(int *)))==NULL) {
    fprintf(stderr," Cannot allocate pam2p: %d\n",len);
    return NULL;
  }

  if((pam2p[0] = (int *)calloc((nsq+1)*len,sizeof(int)))==NULL) {
    fprintf(stderr, "Cannot allocate pam2p[0]: %d\n", (nsq+1)*len);
    free(pam2p);
    return NULL;
  }

  for (i=1; i<len; i++) {
    pam2p[i] = pam2p[0] + (i*(nsq+1));
  }

  return pam2p;
}

void free_pam2p(int **pam2p) {
  if (pam2p) {
    free(pam2p[0]);
    free(pam2p);
  }
}

int
parse_pssm_asn_fa(FILE *afd, int *n_rows, int *n_cols,
		  unsigned char **query, double ***wfreqs,double ***freqs, int ***iscores,
		  char *matrix, int *gap_open, int *gap_extend,
		  double *lambda);

/* the ASN.1 pssm includes information about the scoring matrix used
   (though not the gap penalty in the current version PSSM:2) The PSSM
   scoring matrix and gap penalties should become the default if they
   have not been set explicitly.
*/

/* read the PSSM from an open FILE *fp - but nothing has been read
   from *fp */

int
read_asn_pssm(unsigned char *aa0, int n0, int nsq,
	      double pamscale, FILE *fp, struct pstruct *ppst) {

  int i, j, len, k, itmp;
  int qi, rj;	/* qi - index query; rj - index residues (1-20) */
  int **pam2p;
  int first, too_high;
  char *query, ctmp;
  char dline[512];
  char matrix[MAX_SSTR];
  double psi2_lambda;
  double freq, **wfreq2d, **freq2d, lambda, new_lambda;
  double scale, scale_high, scale_low;
  int **iscores2d;
  int gap_open, gap_extend;
  int n_rows, n_cols;

  pam2p = ppst->pam2p[0];

  if (parse_pssm_asn_fa(fp, &n_rows, &n_cols, (unsigned char **)&query, &wfreq2d, &freq2d, &iscores2d,
			matrix, &gap_open, &gap_extend, &psi2_lambda)<=0) {
    return -1;
  }

  if (!query) { query = (char *)aa0;}

  if (gap_open) {
    if (gap_open > 0) {gap_open = -gap_open;}
    ppst->gdelval = gap_open;
  }
  else if (strncmp(matrix,"BLOSUM62",8)==0) {
    ppst->gdelval = -11;
  }

  if (gap_extend) {
    if (gap_extend > 0) {gap_extend = -gap_extend;}
    ppst->ggapval = gap_extend;
  }
  else if (strncmp(matrix,"BLOSUM62",8)==0) {
    ppst->ggapval = -1;
  }

  if (strncmp(matrix, "BLOSUM62", 8)== 0 && !ppst->pam_set) {
    strncpy(ppst->pamfile, "BP62", 120);
    strncpy(ppst->pamfile_save, ppst->pamfile, 120);
    standard_pam(ppst->pamfile,ppst,0, 0);
    if (!ppst->have_pam2) {
     alloc_pam (MAXSQ, MAXSQ, ppst);
    }
    initpam2(ppst);
    ppst->pam_set = 1;
  }

  if (n_cols < n0) { 
    fprintf(stderr, " query length: %d != n_cols: %d\n",n0, n_cols);
    exit(1);
  }

  printf(" rrtotal: %ld; pamscale=%0.5f\n", rrtotal, pamscale);

  printf("        ");
  for (rj = 1; rj <= N_EFFECT; rj++ ) {
    printf(" %c    ",NCBIstdaa[rj]);
  }
  printf("\n");
  for (qi = 0 ; qi < n0 ; qi++) {
    printf("%3d %c",qi+1, NCBIstdaa[aa0[qi]]);
    for (rj = 0 ; rj < N_EFFECT ; rj++) {
      freq = freq2d[qi][rj];
/*    printf(" %8.7g",freq*10.0); */

      if (freq > 1e-12) {
	freq = log(freq /((double) (rrcounts[rj+1])/(double) rrtotal));
	freq /= pamscale; /* this gets us close to originial pam scores */
	freq2d[qi][rj] = freq;
      }
      printf(" %5.2f",freq);
    }
    printf("\n");
  }

  /* now figure out the right scale */
  scale = 1.0;
  lambda = get_lambda(ppst->pam2[0], 20, 20, "\0ARNDCQEGHILKMFPSTWYV");

  /* should be near 1.0 because of our initial scaling by ppst->pamscale */
  fprintf(stderr, "real_lambda: %g\n", lambda);

  /* get initial high/low scale values: */
  first = 1;
  while (1) {
    fill_pam(pam2p, n0, 20, freq2d, scale);
    new_lambda = get_lambda(pam2p, n0, 20, query); 
    fprintf(stderr, " new_lambda: %g; scale: %g\n",new_lambda,scale);

    if (new_lambda > lambda) {
      if (first) {
	first = 0;
	scale = scale_high = 1.0 + 0.05;
	scale_low = 1.0;
	too_high = 1;
      } else {
	if (!too_high) break;
	scale = (scale_high += scale_high - 1.0);
      }
    } else if (new_lambda > 0) {
      if (first) {
	first = 0;
	scale_high = 1.0;
	scale = scale_low = 1.0 - 0.05;
	too_high = 0;
      } else {
	if (too_high) break;
	scale = (scale_low += scale_low - 1.0);
      }
    } else {
      fprintf(stderr, "new_lambda (%g) <= 0; matrix has positive average score", new_lambda);
      exit(1);
    }
  }

  /* now do binary search between low and high */
  for (i = 0 ; i < 10 ; i++) {
    scale = 0.5 * (scale_high + scale_low);
    fill_pam(pam2p, n0, 20, freq2d, scale);
    new_lambda = get_lambda(pam2p, n0, 20, query);
    fprintf(stderr, "it: %d - new_lambda: %g; scale: %g\n",i,new_lambda,scale);
    
    if (new_lambda > lambda) scale_low = scale;
    else scale_high = scale;
  }

  scale = 0.5 * (scale_high + scale_low);
  fill_pam(pam2p, n0, 20, freq2d, scale);

  fprintf(stderr, "final scale: %g\n", scale);

  fprintf(stderr,"        ");
  for (rj = 1; rj <= N_EFFECT+3; rj++ ) {
    fprintf(stderr,"  %c",NCBIstdaa[rj]);
  }
  fprintf(stderr,"\n");

  for (qi = 0 ; qi < n0 ; qi++) {
    fprintf(stderr, "%4d %c: ", qi+1, NCBIstdaa[aa0[qi]]);
    for (rj = 1 ; rj <= N_EFFECT+3 ; rj++) {
      fprintf(stderr, "%3d", pam2p[qi][rj]);
    }
    fprintf(stderr, "\n");
  }

  if (iscores2d != NULL) {
    fprintf(stderr,"         ");
    for (rj = 1; rj <= N_EFFECT+3; rj++ ) {
      fprintf(stderr," %c ",NCBIstdaa[rj]);
    }
    fprintf(stderr,"\n");
    for (qi = 0 ; qi < n0 ; qi++) {
      fprintf(stderr, "%4d %c: ", qi+1, NCBIstdaa[aa0[qi]]);
      for (rj = 1 ; rj <= N_EFFECT+3 ; rj++) {
	itmp = iscores2d[qi][rj];
	if (itmp < -256) itmp=0;
	fprintf(stderr, "%3d", itmp );
      }
      fprintf(stderr, "\n");
    }
    free(iscores2d[0]);
    free(iscores2d);
  }

  if (wfreq2d != NULL) {
    free(wfreq2d[0]);
    free(wfreq2d);
  }

  if (freq2d != NULL) {
    free(freq2d[0]);
    free(freq2d);
  }

  free(query);
  return 1;
}

void
init_altpam(struct pstruct *ppst) {
  int ix_i, ix_l, ix_j, p_i, p_j, i;

  /* add values for 'J' (I/L) value, which are not present in 1-D matrices */
    ix_i = pascii['I'];
    ix_l = pascii['L'];
    ix_j = pascii['J'];
    if (strchr(pam_sq,'J')==NULL) {
      ppst->pam2[0][ix_j][0] = ppst->pam2[0][0][ix_j] = -BIGNUM;
      /* get the identities */
      ppst->pam2[0][ix_j][ix_j] =
	max(ppst->pam2[0][ix_i][ix_i],ppst->pam2[0][ix_l][ix_l]);
      for (i=1; i < pam_sq_n; i++) {
	p_i = pascii[pam_sq[i]];
	/* do not assume symmetric matrices */
	ppst->pam2[0][ix_j][p_i] =
	  max(ppst->pam2[0][ix_i][p_i],ppst->pam2[0][ix_l][p_i]);
	ppst->pam2[0][p_i][ix_j] =
	  max(ppst->pam2[0][p_i][ix_i],ppst->pam2[0][p_i][ix_l]);
      }
    }
    /* add values for 'O' (K) value, which are not present in 1-D matrices */
    ix_i = pascii['K'];
    ix_j = pascii['O'];  
      if (ix_j < ppst->nsq) {	/* is it in the NCBIstdaa alphabet ? */
      ppst->pam2[0][ix_j][0] = ppst->pam2[0][0][ix_j] = -BIGNUM;
      /* get the identity */
      ppst->pam2[0][ix_j][ix_j] = ppst->pam2[0][ix_i][ix_i];
      /* do not assume symmetric matrices */
      for (i=1; i < pam_sq_n; i++) {
	p_i = pascii[pam_sq[i]];
	ppst->pam2[0][ix_j][p_i] = ppst->pam2[0][ix_i][p_i];
	ppst->pam2[0][p_i][ix_j] = ppst->pam2[0][p_i][ix_i];
      }
    }
    else {
      pascii['O'] = pascii['K'];
      pascii['o'] = pascii['k'];
    }

    /* add values for 'U' (C) value, which are not present in 1-D matrices */
    ix_i = pascii['C'];
    ix_j = pascii['U'];  
      if (ix_j < ppst->nsq) {	/* is it in the NCBIstdaa alphabet */
      ppst->pam2[0][ix_j][0] = ppst->pam2[0][0][ix_j] = -BIGNUM;
      /* get the identity */
      ppst->pam2[0][ix_j][ix_j] = ppst->pam2[0][ix_i][ix_i];
      /* do not assume symmetric matrices */
      for (i=1; i < pam_sq_n; i++) {
	p_i = pascii[pam_sq[i]];
	ppst->pam2[0][ix_j][p_i] = ppst->pam2[0][ix_i][p_i];
	ppst->pam2[0][p_i][ix_j] = ppst->pam2[0][p_i][ix_i];
      }
    }
    else {
      pascii['U'] = pascii['C'];
      pascii['u'] = pascii['c'];
    }
}

