
/*****************************************************************************/
/***  (pseg.c)                                                             ***/
/*** #include precomputed ln(fact(n)) array in "lnfac.h"                   ***/
/*****************************************************************************/

/*--------------------------------------------------------------(includes)---*/

#include "pgenwin.h"
#include "lnfac.h"

/*---------------------------------------------------------------(defines)---*/

#define LENCORRLIM 120
#define MIN(a,b)	((a) <= (b) ? (a) : (b))

/*---------------------------------------------------------------(structs)---*/

struct Segment
  {
   int begin;
   int end;
   struct Segment *next;
  };

/*---------------------------------------------------------------(globals)---*/

extern void entropy_init(int);

/* forward declarations */
void segseq(struct Sequence *, struct Segment **, int);
void mergesegs(struct Sequence *, struct Segment *);
void appendseg(struct Segment *, struct Segment *);

struct Segment *per_mergesegs(struct Sequence *, struct Segment **);
void singreport(struct Sequence	*, struct Segment *);
void per_singreport(struct Sequence *, struct Segment **, struct Segment *);
void prettyreport(struct Sequence *, struct Segment *);
void per_prettyreport(struct Sequence *, struct Segment **, struct Segment *);
void pretreereport(struct Sequence *, struct Segment *); 
void per_pretreereport(struct Sequence *, struct Segment **, struct Segment *);
void report(struct Sequence *, struct Segment *);
void freesegs(struct Segment *);
struct Sequence *seq_phase(struct Sequence *, int, int);
void getparams(int, char **);
void per_seqprep(struct Sequence *, struct Segment **);
void usage();
int findlo(int, int, double *);
int findhi(int, int, double *);
void trim(struct Sequence *, int *, int *);
int hasdash(struct Sequence *);
int min(int, int);
int max(int, int);
double getprob(int *, int );
double lnperm(int *, int);
double lnass(int *);
void seqout(struct Sequence *, int , int , int) ;
void space(int);

int window = 12;
int downset, upset;
double locut = 2.2;
double hicut = 2.5;
int period = 0;

int hilenmin = 0;
int overlaps = FALSE;
int hionly = FALSE;
int loonly = FALSE;
int entinfo = FALSE;
int filterseq = FALSE;
int singleseq = FALSE;
int prettyseq = FALSE;
int prettytree = FALSE;
int charline = 60;
int maxtrim = 100;


/*------------------------------------------------------------------(main)---*/

int main(int argc, char **argv) {

  struct Database *db;
   struct Sequence *seq, *perseq;
   struct Segment *segs, **persegs;
   int ctime;
   int i;

   genwininit();
/* readlencorr(); */                        /* #include lencorr file */
   getparams(argc, argv);
   entropy_init(window);

   if ((db=opendbase(argv[1]))==NULL)
     {
      fprintf(stderr, "Error opening file %s\n", argv[1]);
      exit(1);
     }

   if (period<=0) for (seq=firstseq(db); seq!=NULL; seq=nextseq(db))
     {
      segs = (struct Segment *) NULL;
      segseq(seq, &segs, 0);
      mergesegs(seq, segs);

      if (singleseq) singreport(seq, segs);
      else if (prettyseq) prettyreport(seq, segs);
      else if (prettytree) pretreereport(seq, segs);
      else report(seq, segs);

      freesegs(segs);
      closeseq(seq);
     }
   else for (seq=firstseq(db); seq!=NULL; seq=nextseq(db))
     {
      persegs = (struct Segment **) malloc(period*sizeof(double));
      for (i=0; i<period; i++)
        {
         perseq = seq_phase(seq, i, period);
         persegs[i] = (struct Segment *) NULL;
         segseq(perseq, persegs+i, 0);
         mergesegs(perseq, persegs[i]);
         closeseq(perseq);                                                  
        }

      segs = per_mergesegs(seq, persegs);
      per_seqprep(seq, persegs);

      if (singleseq) per_singreport(seq, persegs, segs);
      else if (prettyseq) per_prettyreport(seq, persegs, segs);
      else if (prettytree) per_pretreereport(seq, persegs, segs);
      else report(seq, segs);

      freesegs(segs);
      for (i=0; i<period; i++) freesegs(persegs[i]);
      free(persegs);
      closeseq(seq);
     }

   closedbase(db);
   exit(0);
  }

/*-------------------------------------------------------------(seq_phase)---*/

struct Sequence *seq_phase(struct Sequence *seq, int phase, int period) {

   struct Sequence *perseq;
   int len, i, j;

   perseq = seqnew();

   len = ((seq->length)/period) + 1;
   perseq->seq = (char *) malloc((len+1)*sizeof(char));

   perseq->length = 0;
   for (i=0, j=phase; j<seq->length; i++, j+=period)
     {
      perseq->seq[i] = seq->seq[j];
      perseq->length++;
     }
   perseq->seq[i] = '\0';

   return(perseq);
  }

/*-------------------------------------------------------------(getparams)---*/

void getparams(int argc, char **argv) {

  int i;
   int nargc;
   char **nargv;
   extern char *optarg;
   extern int optind;
   int c;

   if (argc<2)
     {
      usage();
      exit(1);
     }

   for (i=2; argc>i && argv[i][0]!='-'; i++)
     {
      if (i==2)
        {
         window = atoi(argv[2]);
        }
      else if (i==3)
        {
         locut = atof(argv[3]);
         hicut = locut + 0.3;
        }
      else if (i==4)
        {
         hicut = atof(argv[4]);
        }
     }

   if (locut>hicut)
     {
      fprintf(stderr, "Warning: trigger complexity greater than extension\n");
      hicut = locut;
     }

   downset = (window+1)/2 - 1;
   upset = window - downset;

   if (i==argc) return;

   nargc = argc-i+1;
   nargv = argv+(i-1);
   while ((c=getopt(nargc, nargv, "m:olhaxpqc:nt:z:"))!=-1)
     {
      switch(c)
        {
         case 'z':
            period = atoi(optarg);
            break;
         case 'm':
            hilenmin = atoi(optarg);
            break;
         case 'o':
            overlaps = TRUE;
            hilenmin = 0;
            break;
         case 'l':
            loonly = TRUE;
            break;
         case 'h':
            hionly = TRUE;
            break;
         case 'a':
            hionly = FALSE;
            loonly = FALSE;
            singleseq = FALSE;
            prettyseq = FALSE;
            prettytree = FALSE;
            break;
         case 'x':
            singleseq = TRUE;
            filterseq = TRUE;
            prettyseq = FALSE;
            prettytree = FALSE;
            hionly = TRUE;
            loonly = FALSE;
            break;
         case 'p':
            prettytree = TRUE;
            prettyseq = FALSE;
            singleseq = FALSE;
            filterseq = FALSE;
            break;
         case 'q':
            prettyseq = TRUE;
            prettytree = FALSE;
            singleseq = FALSE;
            filterseq = FALSE;
            break;
         case 'c':
            charline = atoi(optarg);
            break;
         case 'n':
            entinfo = FALSE;
            break;
         case 't':
            maxtrim = atoi(optarg);
            break;
         case '?':
            fprintf(stderr, "Unknown option.\n");
            usage();
            exit(1);
            break;
        }
     }   

   return;
  }

/*----------------------------------------------------------------(segseq)---*/

void segseq(struct Sequence *seq, struct Segment **segs, int offset) {

  struct Segment *seg, *leftsegs;
   struct Sequence *leftseq;
   int first, last, lowlim;
   int loi, hii, i;
   int leftend, rightend, lend, rend;
   double *H, *seqent(struct Sequence *);

   H = seqent(seq);
   if (H==NULL) return;

   first = downset;
   last = seq->length - upset;
   lowlim = first;

   for (i=first; i<=last; i++) {
      if (H[i]<=locut && H[i]!=-1) {
         loi = findlo(i, lowlim, H);
         hii = findhi(i, last, H);

         leftend = loi - downset;
         rightend = hii + upset - 1;

         trim(openwin(seq, leftend, rightend-leftend+1), &leftend, &rightend);

         if (i+upset-1<leftend)   /* check for trigger window in left trim */
           {
            lend = loi - downset;
            rend = leftend - 1;

            leftseq = openwin(seq, lend, rend-lend+1);
            leftsegs = (struct Segment *) NULL;
            segseq(leftseq, &leftsegs, offset+lend);
            if (leftsegs!=NULL)
              {
               if (*segs==NULL) *segs = leftsegs;
               else appendseg(*segs, leftsegs);
              }
            closewin(leftseq);

/*          trim(openwin(seq, lend, rend-lend+1), &lend, &rend);
            seg = (struct Segment *) malloc(sizeof(struct Segment));
            seg->begin = lend;
            seg->end = rend;
            seg->next = (struct Segment *) NULL;
            if (segs==NULL) segs = seg;
            else appendseg(segs, seg);  */
           }

         seg = (struct Segment *) malloc(sizeof(struct Segment));
         seg->begin = leftend + offset;
         seg->end = rightend + offset;
         seg->next = (struct Segment *) NULL;

         if (*segs==NULL) *segs = seg;
         else appendseg(*segs, seg);

         i = MIN(hii, rightend+downset);
         lowlim = i + 1;
/*       i = hii;     this ignores the trimmed residues... */
        }
     }

   free(H);
   return;
  }

/*----------------------------------------------------------------(seqent)---*/

double *seqent(struct Sequence *seq) {

   struct Sequence *win;
   double *H;
   int i, first, last;

   if (window>seq->length)
     {
      return((double *) NULL);
     }

   H = (double *) malloc(seq->length*sizeof(double));

   for (i=0; i<seq->length; i++)
     {
      H[i] = -1.;
     }

   win = openwin(seq, 0, window);
   enton(win);

   first = downset;
   last = seq->length - upset;

   for (i=first; i<=last; i++)
     {
      if (seq->punctuation && hasdash(win))
        {H[i] = -1;
         shiftwin1(win);
         continue;}
      H[i] = win->entropy;
      shiftwin1(win);
     }

   closewin(win);
   return(H);
  }

/*---------------------------------------------------------------(hasdash)---*/

int hasdash(struct Sequence *win) 
{
	register char	*seq, *seqmax;

	seq = win->seq;
	seqmax = seq + win->length;

	while (seq < seqmax) {
		if (*seq++ == '-')
			return TRUE;
	}
	return FALSE;
}

/*----------------------------------------------------------------(findlo)---*/

int findlo(int i, int limit, double *H) {

   int j;

   for (j=i; j>=limit; j--)
     {
      if (H[j]==-1) break;
      if (H[j]>hicut) break;
     }

   return(j+1);
  }

/*----------------------------------------------------------------(findhi)---*/

int findhi(int i, int limit, double *H) {

   int j;

   for (j=i; j<=limit; j++) {
      if (H[j]==-1) break;
      if (H[j]>hicut) break;
     }

   return(j-1);
  }

/*------------------------------------------------------------------(trim)---*/

void trim(struct Sequence *seq, int *leftend, int *rightend) {

   struct Sequence *win;
   double prob, minprob, realprob;
   int shift, len, i;
   int lend, rend;
   int minlen;

/* fprintf(stderr, "%d %d\n", *leftend, *rightend);  */

   lend = 0;
   rend = seq->length - 1;
   minlen = 1;
   if ((seq->length-maxtrim)>minlen) minlen = seq->length-maxtrim;

   minprob = 1.;
/*   fprintf(stderr, "\n");                                         /*-*/
   for (len=seq->length; len>minlen; len--)
     {
/*      fprintf(stderr, "%5d ", len);                               /*-*/
      win = openwin(seq, 0, len);
      i = 0;

      shift = TRUE;
      while (shift)
        {
         prob = getprob(win->state, len);
/*         realprob = exp(prob);                  /*-(for tracing the trim)-*/
/*         fprintf(stderr, "%2e ", realprob);                        /*-*/
         if (prob<minprob)
           {
            minprob = prob;
            lend = i;
            rend = len + i - 1;
           }
         shift = shiftwin1(win);
         i++;
        }
      closewin(win);
/*      fprintf(stderr, "\n");                                     /*-*/
     }

/*   fprintf(stderr, "%d-%d ", *leftend, *rightend);               /*-*/

   *leftend = *leftend + lend;
   *rightend = *rightend - (seq->length - rend - 1);

/*   fprintf(stderr, "%d-%d\n", *leftend, *rightend);              /*-*/

   closewin(seq);
   return;
  }

/*---------------------------------------------------------------(getprob)---*/

double getprob(int *sv, int total) {

  double ans, totseq;

#define LN20	2.9957322735539909
#define LN4	1.3862943611198906

   totseq = ((double) total) * LN20;

   ans = lnass(sv) + lnperm(sv, total) - totseq;

   return(ans);
  }

/*----------------------------------------------------------------(lnperm)---*/

double lnperm(int *sv, int tot) {

   double ans;
   int i;

   ans = lnfac[tot];

   for (i=0; sv[i]!=0; i++) 
     {
      ans -= lnfac[sv[i]];
     }

   return(ans);
  }

/*-----------------------------------------------------------------(lnass)---*/

double lnass(int *sv) {
	double	ans;
	register int	svi, svim1;
	register int	class, total;
	register int    i;

	ans = lnfac[ALPHA_SIZE];
	if (sv[0] == 0)
		return ans;

	total = ALPHA_SIZE;
	class = 1;
	svim1 = sv[0];
	for (i=0;; svim1 = svi) {
	        if (++i==ALPHA_SIZE) {
		        ans -= lnfac[class];
			break;
		      }
		else if ((svi = *++sv) == svim1) {
			class++;
			continue;
		}
		else {
			total -= class;
			ans -= lnfac[class];
			if (svi == 0) {
				ans -= lnfac[total];
				break;
			}
			else {
				class = 1;
				continue;
			}
		}
	}

	return ans;
}

/*-------------------------------------------------------------(mergesegs)---*/

void mergesegs(struct Sequence *seq, struct Segment *segs) {

  struct Segment *seg, *nextseg;
   int len;

   if (overlaps) return;
   if (segs==NULL) return;

   if (segs->begin<hilenmin) segs->begin = 0;

   seg = segs;
   nextseg = seg->next;

   while (nextseg!=NULL)
     {
      if (seg->end>=nextseg->begin && seg->end>=nextseg->end)
        {
         seg->next = nextseg->next;
         free(nextseg);
         nextseg = seg->next;
         continue;
        }
      if (seg->end>=nextseg->begin)               /* overlapping segments */
        {
         seg->end = nextseg->end;
         seg->next = nextseg->next;
         free(nextseg);
         nextseg = seg->next;
         continue;
        }
      len = nextseg->begin - seg->end - 1;
      if (len<hilenmin)                            /* short hient segment */
        {
         seg->end = nextseg->end;
         seg->next = nextseg->next;
         free(nextseg);
         nextseg = seg->next;
         continue;
        }
      seg = nextseg;
      nextseg = seg->next;
     }

   len = seq->length - seg->end - 1;
   if (len<hilenmin) seg->end = seq->length - 1;

   return;
  }

/*---------------------------------------------------------(per_mergesegs)---*/

struct Segment *per_mergesegs(struct Sequence *seq, struct Segment **persegs) 

  {struct Segment **localsegs;
   struct Segment *firstseg, *segs, *seg;
   int first;
   int phase, savephase;

   segs = (struct Segment *) NULL;
   if (persegs==NULL) return(segs);

   localsegs = (struct Segment **) malloc(period*sizeof(double));
   for (phase=0; phase<period; phase++) localsegs[phase] = persegs[phase];

   while (1)
     {
      firstseg = (struct Segment *) NULL;
      first = -1;
      savephase = -1;

      for (phase=0; phase<period; phase++)
        {
         if (localsegs[phase]==NULL) continue;
         if (first==-1 || localsegs[phase]->begin<first)
           {
            savephase = phase;
            firstseg = localsegs[phase];
            first = firstseg->begin;
           }
        }

      if (firstseg==NULL) break;

      seg = (struct Segment *) malloc(sizeof(struct Segment));
      seg->begin = ((firstseg->begin)*period)+savephase;
      seg->end = ((firstseg->end)*period)+savephase;
      seg->next = (struct Segment *) NULL;
      if (segs==NULL) segs = seg; else appendseg(segs, seg);

      localsegs[savephase] = localsegs[savephase]->next;
     }

   free(localsegs);
   mergesegs(seq, segs);
   return(segs);
  }

/*-----------------------------------------------------------(per_seqprep)---*/

void per_seqprep(struct Sequence *seq, struct Segment **persegs) {

   char *proseq, *proseqmax;
   struct Segment *segs, *seg;
   int begin, end, i, phase, pos;

   proseq = seq->seq;
   proseqmax = proseq +seq->length;
   upper(proseq, seq->length);

   for (phase=0; phase<period; phase++)
     {
      segs = persegs[phase];
      for (seg=segs; seg!=NULL; seg=seg->next) 
        {
         begin = seg->begin;
         end = seg->end; 
         for (i=begin; i<=end; i++)
           {
            pos = phase + period*i;
            if (isalpha(proseq[pos])) proseq[pos] = tolower(proseq[pos]);
           }
        }
     }

   if (hionly && singleseq && !prettyseq && !prettytree)
     {
      for (i=0; i<seq->length; i++)
        {
         if (islower(proseq[i])) proseq[i] = 'x';
        }
     }

   if (loonly && prettyseq)
     {
      for (i=0; i<seq->length; i++)
        {
         if (isupper(proseq[i])) proseq[i] = '.';
        }
     }

   if (hionly && prettyseq)
     {
      for (i=0; i<seq->length; i++)
        {
         if (islower(proseq[i])) proseq[i] = '.';
        }
     }

   return;
  }

/*----------------------------------------------------------------(report)---*/

void report(struct Sequence *seq, struct Segment *segs) {

  struct Sequence *subseq;
   struct Segment *seg, *nextseg;
   static int hi = 1;
   static int lo = 0;

   if (segs==NULL)
     {
      enton(seq);
      seqout(seq, hi, 1, seq->length);
/*    fputc('\n', stdout);   -- for spacing after each sequence */
      return;
     }

   if (segs->begin>0)
     {
      subseq = openwin(seq, 0, segs->begin);
      enton(subseq);
      seqout(subseq, hi, 1, segs->begin);
      closewin(subseq);
     }

   for (seg=segs; seg!=NULL; seg=seg->next)
     {
      subseq = openwin(seq, seg->begin, seg->end-seg->begin+1);
      enton(subseq);
      seqout(subseq, lo, seg->begin+1, seg->end+1);
      closewin(subseq);

      if (seg->next==NULL)
        {
         break;
        }

      nextseg = seg->next;
      
      if (nextseg->begin<=seg->end)
        {
         fprintf(stderr, "Overlapping segments: %s\n", seq->id);
         continue;
        }

      if (nextseg->begin==seg->end+1)
        {
         continue;
        }

      subseq = openwin(seq, seg->end+1, nextseg->begin-seg->end-1);
      enton(subseq);
      seqout(subseq, hi, seg->end+2, nextseg->begin);
      closewin(subseq);
     }

   if (seg->end+1==seq->length)
     {
/*    fputc('\n', stdout);   -- for spacing after each sequence */
      return;
     }

   subseq = openwin(seq, seg->end+1, seq->length-seg->end-1);
   enton(subseq);
   seqout(subseq, hi, seg->end+2, seq->length);
   closewin(subseq);

/* fputc('\n', stdout);   -- for spacing after each sequence */
   return;
  }

/*------------------------------------------------------------(singreport)---*/

void singreport(struct Sequence	*seq, struct Segment *segs)
{
	char	*proseq, *proseqmax;
	struct Segment	*seg;
	int	begin, end, i, ctr;

	proseq = seq->seq;
	proseqmax = proseq + seq->length;
	upper(proseq, seq->length);

	for (seg=segs; seg!=NULL; seg=seg->next) {
		begin = seg->begin;
		end = seg->end;
		memset(proseq + begin, 'x', end - begin +1);
	}

	fprintf(stdout, "%s\n", seq->header);

	for (i=0, ctr=0; proseq < proseqmax; ++i, ++ctr, ++proseq) {
		if (ctr==charline) {
			putc('\n', stdout);
			ctr = 0;
		}
		putc(*proseq, stdout);
	}

	putc('\n', stdout);
	if (putc('\n', stdout) == EOF) {
		fprintf(stderr, "premature EOF on write\n");
		exit(2);
	}
}

/*--------------------------------------------------------(per_singreport)---*/

void per_singreport(struct Sequence *seq, struct Segment **persegs, struct Segment *psegs) {

   char *proseq, *proseqmax;
   int i, ctr;

   proseq = seq->seq;
   proseqmax = proseq +seq->length;
   fprintf(stdout, "%s\n", seq->header);

   for (i=0, ctr=0; proseq < proseqmax; ++i, ++ctr, ++proseq) 
     {
      if (ctr==charline) 
        {
         putc('\n', stdout);
         ctr = 0;
        }
      putc(*proseq, stdout);
     }

   putc('\n', stdout);
   if (putc('\n', stdout) == EOF) 
     {
      fprintf(stderr, "premature EOF on write\n");
      exit(2);
     }
}

/*----------------------------------------------------------(prettyreport)---*/

void prettyreport(struct Sequence *seq, struct Segment *segs)

{
	char	*proseq, *proseqmax;
	char	format[10];
	struct Segment	*seg;
	int	begin, end, i, ctr;
	int	leftspace;

	leftspace = (int) ceil(log10((double)seq->length));
	sprintf(format, "%%%dd ", leftspace);

	upper(proseq = seq->seq, seq->length);

	for (seg=segs; seg!=NULL; seg=seg->next) {
		begin = seg->begin;
		end = seg->end;
		lower(proseq + begin, end - begin + 1);
	}

	fprintf(stdout, "%s\n", seq->header);

	space(leftspace+1);
	for (i=0, ctr=1; i<charline; ++i, ++ctr) {
		if (ctr==10) {
			putc('.', stdout);
			ctr = 0;
		}
		else
			putc(' ', stdout);
	}
	putc('\n', stdout);
	fprintf(stdout, format, 1);

	proseqmax = proseq + seq->length;
	for (i=0, ctr=0; proseq < proseqmax; ++i, ++ctr, ++proseq) {
		if (ctr==charline) {
			putc('\n', stdout);
			fprintf(stdout, format, i+1);
			ctr = 0;
		}
		putc(*proseq, stdout);
	}

	fprintf(stdout, " %d\n", seq->length);
	if (putc('\n', stdout) == EOF) {
		fprintf(stderr, "premature EOF on write\n");
		exit(2);
	}
}

/*------------------------------------------------------(per_prettyreport)---*/

void per_prettyreport(struct Sequence *seq, struct Segment **persegs, struct Segment *psegs)   {
   char *proseq, *proseqmax;
   int i, ctr;

   proseq = seq->seq;
   proseqmax = proseq +seq->length;
   fprintf(stdout, "%s\n", seq->header);

   for (i=0, ctr=0; proseq < proseqmax; ++i, ++ctr, ++proseq) 
     {
      if (ctr==charline) 
        {
         putc('\n', stdout);
         ctr = 0;
        }
      putc(*proseq, stdout);
     }

   putc('\n', stdout);
   if (putc('\n', stdout) == EOF) 
     {
      fprintf(stderr, "premature EOF on write\n");
      exit(2);
     }
}

/*---------------------------------------------------------(pretreereport)---*/

void pretreereport(struct Sequence *seq, struct Segment *segs) {

	struct Sequence	*win;
	struct Segment	*seg;
	char	buffer[100], leftfmt[20], rightfmt[20];
	char	*curseq;
	int	i, left, right, len;
	int	current, nextloent;
	int	cline;

	cline = charline / 2;
   
	fprintf(stdout, "%s\n\n", seq->header);
	sprintf(leftfmt, "%%%ds", cline);
	sprintf(rightfmt, "%%-%ds", cline);

	current = 0;

	for (seg=segs; ; seg=seg->next) {
		if (seg==NULL)
			nextloent = seq->length;
		else
			nextloent = seg->begin;

		if (current < nextloent) {
			left = current;
			right = nextloent - 1;
			len = right - left + 1;
			win = openwin(seq, left, len);
			upper(curseq = win->seq, win->length);

			space(cline);
			fprintf(stdout, " %4d-%-4d ", left+1, right+1);

			while (len>0) {
				if (len<=cline) {
					fwrite(curseq, len, 1, stdout);
					putc('\n', stdout);
					break;
				}
				else {
					fwrite(curseq, cline, 1, stdout);
					putc('\n', stdout);
					space(cline+11);
					curseq += cline;
					len -= cline;
				}
			}
			closewin(win);
		}

		if (seg==NULL) break;

      left = seg->begin;
      right = seg->end;
      len = right - left + 1;
      win = openwin(seq, left, len);
      lower(curseq = win->seq, win->length);

		i = MIN(cline, len);
		if (i < cline)
			fprintf(stdout, "%*s", cline-i, "");
		fwrite(curseq, i, 1, stdout);
		fprintf(stdout, " %4d-%-4d ", left+1, right+1);
		putc('\n', stdout);

		len -= cline;
		if (len>0)
			curseq += cline;

		while (len>0) {
			i = MIN(cline, len);
			if (i < cline)
				space(cline - i);
			fwrite(curseq, i, 1, stdout);
			putc('\n', stdout);
			len -= i;
			if (len>0)
				curseq += i;
		}

		closewin(win);
		current = right + 1;
	}

	putc('\n', stdout);
}

/*-----------------------------------------------------(per_pretreereport)---*/

void per_pretreereport(struct Sequence *seq, struct Segment **persegs, struct Segment *segs) 

  {
	struct Sequence	*win;
	struct Segment	*seg;
	char	buffer[100], leftfmt[20], rightfmt[20];
	char	*curseq;
	int	i, left, right, len;
	int	current, nextloent;
	int	cline;

	cline = charline / 2;
   
	fprintf(stdout, "%s\n\n", seq->header);
	sprintf(leftfmt, "%%%ds", cline);
	sprintf(rightfmt, "%%-%ds", cline);

	current = 0;

	for (seg=segs; ; seg=seg->next) {
		if (seg==NULL)
			nextloent = seq->length;
		else
			nextloent = seg->begin;

		if (current < nextloent) {
			left = current;
			right = nextloent - 1;
			len = right - left + 1;
			win = openwin(seq, left, len);
                        curseq = win->seq;

			space(cline);
			fprintf(stdout, " %4d-%-4d ", left+1, right+1);

			while (len>0) {
				if (len<=cline) {
					fwrite(curseq, len, 1, stdout);
					putc('\n', stdout);
					break;
				}
				else {
					fwrite(curseq, cline, 1, stdout);
					putc('\n', stdout);
					space(cline+11);
					curseq += cline;
					len -= cline;
				}
			}
			closewin(win);
		}

		if (seg==NULL) break;

      left = seg->begin;
      right = seg->end;
      len = right - left + 1;
      win = openwin(seq, left, len);
      curseq = win->seq;

		i = MIN(cline, len);
		if (i < cline)
			fprintf(stdout, "%*s", cline-i, "");
		fwrite(curseq, i, 1, stdout);
		fprintf(stdout, " %4d-%-4d ", left+1, right+1);
		putc('\n', stdout);

		len -= cline;
		if (len>0)
			curseq += cline;

		while (len>0) {
			i = MIN(cline, len);
			if (i < cline)
				space(cline - i);
			fwrite(curseq, i, 1, stdout);
			putc('\n', stdout);
			len -= i;
			if (len>0)
				curseq += i;
		}

		closewin(win);
		current = right + 1;
	}

	putc('\n', stdout);

  }

/*-----------------------------------------------------------------(space)---*/

void space(int len)
{
	static char	spaces[] =
		{' ',' ',' ',' ',' ',' ',' ',' ',' ',' ',' ',' ',' ',' ',' ',' ',' ',
		' ',' ',' ',' ',' ',' ',' ',' ',' ',' ',' ',' ',' ',' ',' ',' ',' ',
		' ',' ',' ',' ',' ',' ',' ',' ',' ',' ',' ',' ',' ',' ',' ',' ',' ',
		' ',' ',' ',' ',' ',' ',' ',' ',' ',' ',' ',' ',' ',' ',' ',' ',' ',
		' ',' ',' ',' ',' ',' ',' ',' ',' ',' ',' ',' ',' ',' ',' ',' ',' '
		};
	register int i;

	while (len > 0) {
		i = MIN(len, sizeof(spaces)/sizeof(spaces[0]));
		fwrite(spaces, i, 1, stdout);
		len -= i;
	}
}

/*----------------------------------------------------------------(seqout)---*/

#define HDRLEN 60
void seqout(struct Sequence *seq, int hilo, int begin, int end) {

   char *proseq, *proseqmax, *id, *header;
   char outbuf[HDRLEN+1];
   static int hi = 1;
   static int lo = 0;
   int i, ctr, iend;

   if (hionly && hilo==lo) return;
   if (loonly && hilo==hi) return;

   proseq = seq->seq;
   proseqmax = proseq + seq->length;
   id = seq->id;
   if (id==NULL) id = seq->parent->id;
   header = seq->header;
   if (header==NULL) header = seq->parent->header;

   iend = findchar(header, ' ');
   if (iend!=-1) header = header+iend;

   if (entinfo)
     {
      fprintf(stdout, ">%s(%d-%d)", id, begin, end);
/*    if (iend!=-1 && strlen(header)<=HDRLEN) fprintf(stdout, "%s", header);
      else if (iend!=-1) for (i=0; i<HDRLEN; i++) putc(header[i], stdout); */
      fprintf(stdout, " complexity=%4.2f (%d/%4.2f/%4.2f)\n",
           seq->entropy, window, locut, hicut);
     }
   else
     {
      fprintf(stdout, ">%s(%d-%d)", id, begin, end);
      if (iend!=-1)   /* fprintf(stdout, "%s\n", header); */
        {
		 i = MIN(HDRLEN, strlen(header));
		 fwrite(header, i, 1, stdout);
         putc('\n', stdout);
        }
      else putc('\n', stdout);
     }
   
/*   if (hilo==lo)
/*     {
/*      lower(proseq, seq->length);
/*     }
/*   else if (hilo==hi && seq->length>=hilenmin)
/*     {
/*      upper(proseq, seq->length);
/*     }
/*   else
/*     {
/*      lower(proseq, seq->length);
/*     }                                               */

   for (; proseq < proseqmax; proseq+=i) {
		i = MIN(charline, proseqmax - proseq);
		fwrite(proseq, i, 1, stdout);
		putc('\n', stdout);
	}

	if (putc('\n', stdout) == EOF) {
		fprintf(stderr, "premature EOF on write\n");
		exit(2);
	}
}

/*-------------------------------------------------------------(appendseg)---*/

void appendseg(struct Segment *segs, struct Segment *seg) {

   struct Segment *temp;

   temp = segs;
   while (1)
     {
      if (temp->next==NULL)
        {
         temp->next = seg;
         break;
        }
      else
        {
         temp = temp->next;
        }
     }

   return;
  }

/*--------------------------------------------------------------(freesegs)---*/

void freesegs(struct Segment *segs) {

  struct Segment *temp;

   while (segs!=NULL)
     {
      temp = segs->next;
      free(segs);
      segs = temp;
     }
  }

/*-----------------------------------------------------------------(usage)---*/

void usage()

  {
   fprintf(stderr, "\
Usage: pseg <file> <window> <locut> <hicut> <options>\n\
         <file>   - filename containing fasta-formatted sequence(s) \n\
         <window> - OPTIONAL window size (default 12) \n\
         <locut>  - OPTIONAL low (trigger) complexity (default 2.2) \n\
         <hicut>  - OPTIONAL high (extension) complexity (default locut + 0.3) \n\
	 <options> \n\
            -x  each input sequence is represented by a single output \n\
                sequence with low-complexity regions replaced by \n\
                strings of 'x' characters \n\
            -c <chars> number of sequence characters/line (default 60)\n\
            -m <size> minimum length for a high-complexity segment \n\
                (default 0).  Shorter segments are merged with adjacent \n\
                low-complexity segments \n\
            -l  show only low-complexity segments (fasta format) \n\
            -h  show only high-complexity segments (fasta format) \n\
            -a  show all segments (fasta format) \n\
            -n  do not add complexity information to the header line \n\
            -o  show overlapping low-complexity segments (default merge) \n\
            -t <maxtrim> maximum trimming of raw segment (default 100) \n\
            -p  prettyprint each segmented sequence (tree format) \n\
            -q  prettyprint each segmented sequence (block format) \n\
            -z <period>\n");
   exit(1);
  }

