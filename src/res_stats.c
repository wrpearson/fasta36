/* $Id: res_stats.c 1227 2013-09-26 19:19:28Z wrp $  */

/* copyright (C) 2005, 2014 by William R. Pearson and The Rector &
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

/* calculate stats from results file using scalesws.c */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <limits.h>
#include <math.h>

#define MAX_LLEN 200

#define LN_FACT 10.0

#include "defs.h"
#include "structs.h"
#include "param.h"

struct beststr {
  int score;	/* smith-waterman score */
  int sscore;	/* duplicate for compatibility with fasta */
  double comp;
  double H;
  double zscore;
  double escore;
  int n1;
#ifndef USE_FTELLO
  long lseek;	/* position in library file */
#else
  off_t lseek;
#endif
  int cont;	/* offset into sequence */
  int frame;
  int lib;
  char libstr[13];
} *bbp, *bestptr, **bptr, *best;

struct stat_str {
  int score;
  int n1;
  double comp;
  double H;
};

static struct db_str qtt = {0l, 0l, 0};

char info_gstring2[MAX_STR];                  /* string for label */
char info_gstring3[MAX_STR];
char info_hstring1[MAX_STR];

FILE *outfd;

int nbest;	/* number of sequences better than bestcut in best */
int bestcut=1; 	/* cut off for getting into MAX_BEST */
int bestfull;

int dohist = 0;
int zsflag = 1;
int outtty=1;
int llen=40;

/* statistics functions */
extern void
process_hist(struct stat_str *sptr, int nstat, struct pstruct pst,
	     struct hist_str *hist, void **);
extern void addhistz(double, struct hist_str *); /* scaleswn.c */
void selectbestz(struct beststr **, int, int );

extern double zs_to_E(double, int, int, long, struct db_str);
extern double zs_to_Ec(double zs, long entries);

extern double find_z(int score, int length, double comp, void *);

void prhist(FILE *, struct mngmsg, struct pstruct, struct hist_str, 
	    int, struct db_str, char *);

int nshow=20, mshow=50, ashow= -1;
double e_cut=10.0;

main(argc, argv)
     int argc; char **argv;
{
  FILE *fin;
  char line[512];
  int max, icol, iarg, i, qsfnum, lsfnum, n0, n1, s[3], frame;
  double comp, H;
  int idup, ndup, max_s;
  char libstr[MAX_UID], *bp;
  char bin_file[80];
  FILE *bout=NULL;
  struct mngmsg m_msg;		/* Message from host to manager */
  struct pstruct pst;
  struct stat_str *stats;
  int nstats;
  double zscor, mu, var;

#if defined(UNIX)
  outtty = isatty(1);
#else
  outtty = 1;
#endif

  if (argc < 2 ) {
    fprintf(stderr," useage - res_stats -c col -r bin_file file\n");
    exit(1);
  }

  m_msg.db.length = qtt.length = 0l;
  m_msg.db.entries = m_msg.db.carry = qtt.entries = qtt.carry = 0;
  m_msg.pstat_void = NULL;
  m_msg.hist.hist_a = NULL;
  m_msg.nohist = 0;
  m_msg.markx = 0;

  pst.n0 = 200;		/* sensible dummy value */
  pst.zsflag = 1;
  pst.dnaseq = 0;
  pst.histint = 2;

  bin_file[0]='\0';
  icol = 1;
  iarg = 1;
  ndup = 1;
  while (1) {
    if (argv[iarg][0]=='-' && argv[iarg][1]=='c') {
      sscanf(argv[iarg+1],"%d",&icol);
      iarg += 2;
    }
    else if (argv[iarg][0]=='-' && argv[iarg][1]=='r') {
      strncpy(bin_file,argv[iarg+1],sizeof(bin_file));
      iarg += 2;
    }
    else if (argv[iarg][0]=='-' && argv[iarg][1]=='z') {
      sscanf(argv[iarg+1],"%d",&pst.zsflag);
      iarg += 2;
    }
    else if (argv[iarg][0]=='-' && argv[iarg][1]=='n') {
      pst.dnaseq = 1;
      iarg += 1;
    }
    else if (argv[iarg][0]=='-' && argv[iarg][1]=='s') {
      sscanf(argv[iarg+1],"%d",&ndup);
      iarg += 2;
    }
    else if (argv[iarg][0]=='-' && argv[iarg][1]=='q') {
      outtty = 0;
      iarg += 1;
    }
    else break;
  }

  icol--;

  if ((fin=fopen(argv[iarg],"r"))==NULL) {
    fprintf(stderr," cannot open %s\n",argv[1]);
    exit(1);
  }

  if (bin_file[0]!='\0' && ((bout=fopen(bin_file,"w"))==NULL)) {
    fprintf(stderr,"cannot open %s for output\n",bin_file);
  }

  if ((stats =
       (struct stat_str *)malloc((MAX_STATS)*sizeof(struct stat_str)))==NULL)
    s_abort ("Cannot allocate stats struct","");
  nstats = 0;

  initbest(MAX_BEST+1);	/* +1 required for select() */

  for (nbest=0; nbest<MAX_BEST+1; nbest++)
    bptr[nbest] = &best[nbest];
  bptr++; best++;
  best[-1].score= BIGNUM;
  
  nbest = 0;

  pst.Lambda=0.232;
  pst.K = 0.11;
  pst.H = 0.34;

  /* read the best scores from the results file */

  max_s = -1;
  idup = 0;

  /* get first line with sequence length */
  fgets(line,sizeof(line),fin);
  sscanf(line,"%d",&n0);
  if (n0 > 0) pst.n0 = n0;

  while (fgets(line,sizeof(line),fin)!=NULL) {
    if (line[0]=='/' && line[1]=='*') {
      fputs(line,stdout);
      strncpy(info_gstring2,line,sizeof(info_gstring2));
      if ((bp=strchr(info_gstring2,'\n'))!=NULL) *bp = '\0';
      break;
    }
    if (line[0]==';') {
      if ((bp=strchr(line,'|'))!=NULL) qsfnum = atoi(bp+1);
      else continue;
      if ((bp=strchr(line,'('))!=NULL) {
	n0 = atoi(bp+1);
	pst.n0 = n0;
      }
      else {
	fprintf(stderr, "cannot find n0:\n %s\n",line);
	continue;
      }
    }
    else {
	sscanf(line,"%s %d %d %d %lf %lf %d %d %d",
	       libstr,&lsfnum,&n1,&frame,&comp, &H, &s[0],&s[1],&s[2]);
	if (lsfnum==0 && n1==0) {
	  fputs(line,stderr);
	  continue;
	}
	if (n1 < 10 || s[icol]<=0) fputs(line,stderr);
	idup++;

	if (s[icol] > max_s) max_s = s[icol];
	if (idup < ndup) continue;

	m_msg.db.entries++;
	m_msg.db.length += n1;

	if (dohist) addhistz(zscor=find_z(max_s,n1,comp,m_msg.pstat_void),
			     &m_msg.hist);
	else zscor = (double)max_s;

	if (nstats < MAX_STATS) {
	  stats[nstats].n1 = n1;
	  stats[nstats].comp = comp;
	  stats[nstats].H = H;
	  stats[nstats++].score = max_s;
	}

	else if (!dohist) {
	  /*	  do_bout(bout,stats,nstats); */
	  process_hist(stats,nstats,pst,&m_msg.hist, &m_msg.pstat_void);
	  for (i=0; i<nbest; i++)
	    bptr[i]->zscore = 
	      find_z(bptr[i]->score,bptr[i]->n1,bptr[i]->comp,
			 m_msg.pstat_void);
	  dohist = 1;
	}

	if (dohist) {
	  zscor =find_z(max_s,n1,comp,m_msg.pstat_void);
	  addhistz(zscor,&m_msg.hist);
	}
	else zscor = (double)max_s;

	if (nbest >= MAX_BEST) {
	  bestfull = nbest-MAX_BEST/4;
	  selectz(bestfull-1,nbest);
	  bestcut = (int)(bptr[bestfull-1]->zscore+0.5);
	  nbest = bestfull;
	}
	bestptr = bptr[nbest];
	bestptr->score = max_s;
	bestptr->sscore = max_s;
	bestptr->n1 = n1;
	bestptr->comp = comp;
	bestptr->H = H;
	bestptr->lib = lsfnum;
	bestptr->zscore = zscor;
	strncpy(bestptr->libstr,libstr,12);
	bestptr->libstr[12]='\0';
	nbest++;

	max_s = -1;
	idup = 0;
    }
  }	/* done with reading results */

  if (!dohist) {
    if (nbest < 20) {
      zsflag = 0;
    }
    else {
      /*      do_bout(bout,stats,nstats); */
      process_hist(stats,nstats,pst,&m_msg.hist,&m_msg.pstat_void);
      for (i=0; i<nbest; i++)
	bptr[i]->zscore = 
	  find_z(bptr[i]->score,bptr[i]->n1,bptr[i]->comp,m_msg.pstat_void);
      dohist = 1;
    }
  }
  
  printf(" using n0: %d\n",pst.n0);

  /* print histogram, statistics */

  m_msg.nbr_seq = m_msg.db.entries;
  pst.zdb_size = m_msg.db.entries;
  /* get_param(&pst, info_gstring2,info_gstring3); */

  prhist(stdout,m_msg,pst,m_msg.hist,nstats,m_msg.db,info_gstring2);

  if (!zsflag) sortbest();
  else {
    sortbestz(bptr,nbest);
    for (i=0; i<nbest; i++)
      bptr[i]->escore = zs_to_E(bptr[i]->zscore,bptr[i]->n1,pst.dnaseq,
				pst.zdb_size, m_msg.db);
  }
  
  outfd = stdout;
  showbest(m_msg.db);	/* display best matches */
}

initbest(nbest)		/* allocate arrays for best sort */
     int nbest;
{

  if ((best=(struct beststr *)calloc((size_t)nbest,sizeof(struct beststr)))
      == NULL) {fprintf(stderr,"cannot allocate best struct\n"); exit(1);}
  if ((bptr=(struct beststr **)calloc((size_t)nbest,sizeof(struct beststr *)))
      == NULL) {fprintf(stderr,"cannot allocate bptr\n"); exit(1);}
}

void
prhist(FILE *fd, struct mngmsg m_msg,
       struct pstruct pst, 
       struct hist_str hist, 
       int nstats,
       struct db_str ntt,
       char *info_gstring2)
{
  int i,j,hl,hll, el, ell, ev;
  char hline[80], pch, *bp;
  int mh1, mht;
  int maxval, maxvalt, dotsiz, ddotsiz,doinset;
  double cur_e, prev_e, f_int;
  double max_dev, x_tmp;
  double db_tt;
  int n_chi_sq, cum_hl, max_i;


  fprintf(fd,"\n");
  
  if (pst.zsflag < 0 || nstats <= 10) {
    fprintf(fd, "%7ld residues in %5ld sequences\n", ntt.length,ntt.entries);
    fprintf(fd,"\n%s\n",info_gstring2);
    return;
  }

  max_dev = 0.0;
  mh1 = hist.maxh-1;
  mht = (3*hist.maxh-3)/4 - 1;

  if (!m_msg.nohist && mh1 > 0) {
    for (i=0,maxval=0,maxvalt=0; i<hist.maxh; i++) {
      if (hist.hist_a[i] > maxval) maxval = hist.hist_a[i];
      if (i >= mht &&  hist.hist_a[i]>maxvalt) maxvalt = hist.hist_a[i];
    }
    n_chi_sq = 0;
    cum_hl = -hist.hist_a[0];
    dotsiz = (maxval-1)/60+1;
    ddotsiz = (maxvalt-1)/50+1;
    doinset = (ddotsiz < dotsiz && dotsiz > 2);

    if (pst.zsflag>=0)
      fprintf(fd,"       opt      E()\n");
    else 
      fprintf(fd,"     opt\n");

    prev_e =  zs_to_Ec((double)(hist.min_hist-hist.histint/2),hist.entries);
    for (i=0; i<=mh1; i++) {
      pch = (i==mh1) ? '>' : ' ';
      pch = (i==0) ? '<' : pch;
      hll = hl = hist.hist_a[i];
      if (pst.zsflag>=0) {
	cum_hl += hl;
	f_int = (double)(i*hist.histint+hist.min_hist)+(double)hist.histint/2.0;
	cur_e = (double)zs_to_Ec(f_int,hist.entries);
	ev = el = ell = (int)(cur_e - prev_e + 0.5);
	if (hl > 0  && i > 5 && i < (90-hist.min_hist)/hist.histint) {
	  x_tmp  = fabs(cum_hl - cur_e);
	  if ( x_tmp > max_dev) {
	    max_dev = x_tmp;
	    max_i = i;
	  }
	  n_chi_sq++;
	}
	if ((el=(el+dotsiz-1)/dotsiz) > 60) el = 60;
	if ((ell=(ell+ddotsiz-1)/ddotsiz) > 40) ell = 40;
	fprintf(fd,"%c%3d %5d %5d:",
		pch,(i<mh1)?(i)*hist.histint+hist.min_hist :
		mh1*hist.histint+hist.min_hist,hl,ev);
      }
      else fprintf(fd,"%c%3d %5d :",
		   pch,(i<mh1)?(i)*hist.histint+hist.min_hist :
		   mh1*hist.histint+hist.min_hist,hl);

      if ((hl=(hl+dotsiz-1)/dotsiz) > 60) hl = 60;
      if ((hll=(hll+ddotsiz-1)/ddotsiz) > 40) hll = 40;
      for (j=0; j<hl; j++) hline[j]='='; 
      if (pst.zsflag>=0) {
	if (el <= hl ) {
	  if (el > 0) hline[el-1]='*';
	  hline[hl]='\0';
	}
	else {
	  for (j = hl; j < el; j++) hline[j]=' ';
	  hline[el-1]='*';
	  hline[hl=el]='\0';
	}
      }
      else hline[hl] = 0;
      if (i==1) {
	for (j=hl; j<10; j++) hline[j]=' ';
	sprintf(&hline[10]," one = represents %d library sequences",dotsiz);
      }
      if (doinset && i == mht-2) {
	for (j = hl; j < 10; j++) hline[j]=' ';
	sprintf(&hline[10]," inset = represents %d library sequences",ddotsiz);
      }
      if (i >= mht&& doinset ) {
	for (j = hl; j < 10; j++) hline[j]=' ';
	hline[10]=':';
	for (j = 11; j<11+hll; j++) hline[j]='=';
	hline[11+hll]='\0';
	if (pst.zsflag>=0) {
	  if (ell <= hll) hline[10+ell]='*';
	  else {
	    for (j = 11+hll; j < 10+ell; j++) hline[j]=' ';
	    hline[10+ell] = '*';
	    hline[11+ell] = '\0';
	  }
	}
      }

      fprintf(fd,"%s\n",hline);
      prev_e = cur_e;
    }
  }

  if (ntt.carry==0) {
    fprintf(fd, "%7ld residues in %5ld sequences\n", ntt.length, ntt.entries);
  }
  else {
    db_tt = (double)ntt.carry*(double)LONG_MAX + (double)ntt.length;
    fprintf(fd, "%.0f residues in %5ld library sequences\n", db_tt, ntt.entries);
  }

  if (pst.zsflag>=0) {
    if (MAX_STATS < hist.entries)
      fprintf(fd," statistics extrapolated from %d to %ld sequences\n",
	      MAX_STATS,hist.entries);
    /*    summ_stats(stat_info); */
    fprintf(fd," %s\n",hist.stat_info);
    if (!m_msg.nohist && cum_hl > 0)
      fprintf(fd," Kolmogorov-Smirnov  statistic: %6.4f (N=%d) at %3d\n",
	      max_dev/(double)cum_hl, n_chi_sq,max_i*hist.histint+hist.min_hist);
    if (m_msg.markx & MX_M10FORM) {
      while ((bp=strchr(hist.stat_info,'\n'))!=NULL) *bp=' ';
      if (cum_hl <= 0) cum_hl = -1;
      sprintf(info_hstring1,"; mp_extrap: %d %ld\n; mp_stats: %s\n; mp_KS: %6.4f (N=%d) at %3d\n",
	      MAX_STATS,hist.entries,hist.stat_info,max_dev/(double)cum_hl, n_chi_sq,max_i*hist.histint+hist.min_hist);
    }
  }
  fprintf(fd,"\n%s\n",info_gstring2);
  fflush(fd);
}

showbest(struct db_str ntt)
  {
    int ib, istart, istop;
    char bline[200], fmt[40], pad[200];
    char rline[20];
    int ntmp;
    int lcont, ccont, loff;
    int hcutoff;

    sprintf(fmt,"%%-%ds (%%3d)",llen-10);

    nshow = min(20,nbest);
    mshow = min(20,nbest);

    if (outtty) {
      printf(" How many scores would you like to see? [%d] ",nshow);
      fflush(stdout);
      if (fgets(rline,sizeof(rline),stdin)==NULL) exit(0);
      if (rline[0]!='\n' && rline[0]!=0) sscanf(rline,"%d",&nshow);
      if (nshow<=0) nshow = min(20,nbest);
    }
    else nshow=mshow;

    memset(pad,' ',llen-10);
    pad[llen-31]='\0';
    if (zsflag)
      fprintf(outfd,"The best scores are:%s s-w Z-score E(%ld)\n",pad,ntt.entries);
    else
      fprintf(outfd,"The best scores are:%s s-w\n",pad);

    if (outfd != stdout)
      if (zsflag)
	fprintf(stdout,"The best scores are:%s s-w Z-score E(%ld)\n",pad,ntt.entries);
      else
	fprintf(stdout,"The best scores are:%s s-w\n",pad);

    istart = 0;
  l1:	istop = min(nbest,nshow);
  for (ib=istart; ib<istop; ib++) {
    bbp = bptr[ib];

    if (!outtty && zsflag && bbp->escore > e_cut) {
      nshow = ib;
      goto done;
    }

    sprintf(bline,"%-12s %d",bbp->libstr,bbp->lib);
    bline[13]='\0';

    fprintf(outfd,fmt,bline,bbp->n1);

    if (zsflag)
      fprintf(outfd,"%4d %4.1f %6.2g\n",
	      bbp->score,bbp->zscore,
	      bbp->escore);
    else 
      fprintf(outfd,"%4d\n",bbp->score);

    if (outfd!=stdout) {
      fprintf(stdout,fmt,bline,bbp->n1);
      if (zsflag)
	printf("%4d %4.1f %6.2g\n",
	       bbp->score,bbp->zscore,
	       bbp->escore);
      else 
	printf("%4d\n",bbp->score);
    }
  }

  fflush(outfd); if (outfd!=stdout) fflush(stdout);

  if (outtty) {
    printf(" More scores? [0] ");
    fflush(stdout);
    if (fgets(rline,sizeof(rline),stdin)==NULL) exit(0);
    ntmp = 0;
    if (rline[0]!='\n' && rline[0]!=0) sscanf(rline,"%d",&ntmp);
    if (ntmp<=0) ntmp = 0;
    if (ntmp>0) {
      istart = istop;
      nshow += ntmp;
      mshow += ntmp;
      goto l1;
    }
  }
  else if (zsflag && bbp->escore < e_cut) {
    istart=istop;
    nshow += 10;
    goto l1;
  }

  done:
  if (outfd!=stdout) fprintf(outfd,"\n");
}

selectz(k,n)	/* k is rank in array */
     int k,n;
{
  int t, i, j, l, r;
  double v;
  struct beststr *tmptr;

  l=0; r=n-1;

  while ( r > l ) {
    i = l-1;
    j = r;
    v = bptr[r]->zscore;
    do {
      while (bptr[++i]->zscore > v ) ;
      while (bptr[--j]->zscore < v ) ;
      tmptr = bptr[i]; bptr[i]=bptr[j]; bptr[j]=tmptr;
    } while (j > i);
    bptr[j]=bptr[i]; bptr[i]=bptr[r]; bptr[r]=tmptr;
    if (i>=k) r = i-1;
    if (i<=k) l = i+1;
  }
}

sortbest()
{
  int cmps(), cmp1(), cmpa(), cmpz();
  ksort(bptr,nbest,cmps);
}

sortbeste()
{
  int cmpe();
  ksort(bptr,nbest,cmpe);
}

sortbestz()
{
  int cmpz();
  ksort(bptr,nbest,cmpz);
}

cmps(ptr1,ptr2)
     struct beststr *ptr1, *ptr2;
{
  if (ptr1->score < ptr2->score) return (1);
  else if (ptr1->score > ptr2->score) return (-1);
  else return (0);
}

cmpe(ptr1,ptr2)
     struct beststr *ptr1, *ptr2;
{
  if (ptr1->escore < ptr2->escore) return (-1);
  else if (ptr1->escore > ptr2->escore) return (1);
  else return (0);
}

cmpz(ptr1,ptr2)
     struct beststr *ptr1, *ptr2;
{
  if (ptr1->zscore < ptr2->zscore) return (1);
  else if (ptr1->zscore > ptr2->zscore) return (-1);
  else return (0);
}

ksort(v,n,comp)
     char *v[]; int n, (*comp)();
{
  int gap, i, j;
  char *tmp;
	
  for (gap=n/2; gap>0; gap/=2)
    for (i=gap; i<n; i++)
      for (j=i-gap; j>=0; j -= gap) {
	if ((*comp)(v[j],v[j+gap]) <=0)
	  break;
	tmp = v[j]; v[j]=v[j+gap]; v[j+gap]=tmp;
      }
}

/*
do_bout(FILE *bout,struct stat_str **bptr, int nbest)
{
  int i, min_hist, max_hist;
  double mu, var;

  if (bout==NULL) return;

  inithist();
  for (i = 0; i<nbest; i++)
    addhist(bptr[i]->score,bptr[i]->n1);

  for (i=0; i<MAX_LLEN; i++)
    if (llen_hist[i]>0) {
      min_hist=i;
      break;
    }

  for (i=MAX_LLEN-1; i>=0; i--)
    if (llen_hist[i]>0) {
      max_hist=i;
      break;
    }

  for (i=min_hist; i<=max_hist; i++) {
    mu=(double)score_sums[i]/(double)llen_hist[i];
    if (llen_hist[i]>1) {
      var = ((double)score2_sums[i]-(double)llen_hist[i]*mu*mu)/
	(double)(llen_hist[i]-1);

      fprintf(bout,"%d\t%d\t%.1f\t%.1f\t%.1f\t%.4f\t%.4f\n",
	      i,llen_hist[i],exp(((double)(i))/LN_FACT),
	      score_sums[i],score2_sums[i],mu,var);
    }
  }
  free_hist();
  fclose(bout);
}
*/

s_abort()
{
  exit(1);
}
