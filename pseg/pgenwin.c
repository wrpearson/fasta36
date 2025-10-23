
/*****************************************************************************/
/***   (pgenwin.c)                                                         ***/
/*****************************************************************************/

/*--------------------------------------------------------------(includes)---*/

#include "pgenwin.h"

/*---------------------------------------------------------------(defines)---*/

#define STRSIZE 100

/*----------------------------------------------------------------(protos)---*/

struct Sequence *readentry(struct Database *);
/*---------------------------------------------------------------(globals)---*/

char *blastdbs[] =
  {"bba", "bbn", "embl", "gbupdate", "genbank", "genpept", "gpupdate",
   "nr", "nrdb", "nrdb.shuf", "pir", "pseq", "swissprot", "tfdaa"};

int nblastdbs = 14;

int blastdb(char *);

#ifndef MIN
#define MIN(a,b) ((a) <= (b) ? (a) : (b))
#endif

char *blastdir = "/net/cruncher/usr/ncbi/db/fasta/";
char *indexdir = "/net/cruncher/usr/ncbi/db/index/";

int nabets;
struct Alphabet **abets;
int ntvecs;
struct TransVector **tvecs;
int nsvecs;
struct ScoreVector **svecs;
int nsmats;
struct ScoreMatrix **smats;

int aaindex[128];
unsigned char	aaflag[128];
char aachar[ALPHA_SIZE];

struct strlist
  {
   char string[STRSIZE];
   struct strlist *next;
  } *str, *curstr;

/*---------------------------------------------------------------(tmalloc)---*/

#define TESTMAX 1000
void *tmalloc(size_t);
size_t record_ptrs[TESTMAX] = {0,0,0,0};
int rptr = 0;

/* forward declarations */
static void readseq(struct Sequence *);
static void decrementsv(int *sv, int class), incrementsv(int *sv, int class);
static void skipline(FILE *);
extern int findchar(char *, char);
static int readhdr(struct Sequence *);
static int min(int, int);

/*------------------------------------------------------------(genwininit)---*/

void genwininit()
{
	char	*cp, *cp0;
	int		i;
	char	c;

	for (i = 0; i < sizeof(aaindex)/sizeof(aaindex[0]); ++i) {
		aaindex[i] = ALPHA_SIZE;
		aaflag[i] = TRUE;
	}
		
	for (cp = cp0 = "ACDEFGHIKLMNPQRSTVWY"; (c = *cp) != '\0'; ++cp) {
		i = cp - cp0;
		aaindex[c] = i;
		aaindex[tolower(c)] = i;
		aachar[i] = tolower(c);
		aaflag[c] = FALSE;
		aaflag[tolower(c)] = FALSE;
	}
	return;
}
        
/*-------------------------------------------------------------(opendbase)---*/

extern struct Database *opendbase(char *name) {
 
  struct Database *dbase;

   dbase = (struct Database *) malloc(sizeof(struct Database));

   if (blastdb(name))
     {
      dbase->filename = (char *) malloc(strlen(blastdir)+strlen(name)+1);
      dbase->indexname = (char *) malloc(strlen(indexdir)+strlen(name)+1);
      strcpy(dbase->filename, blastdir);
      strcat(dbase->filename, name);
      strcpy(dbase->indexname, indexdir);
      strcat(dbase->indexname, name);
     }
   else
     {
      dbase->filename = (char *) malloc(strlen(name)+1);
      dbase->indexname = (char *) malloc(strlen(name)+1);
      strcpy(dbase->filename, name);
      strcpy(dbase->indexname, name);
     }

   if (strcmp(dbase->filename, "-")==0)
     {
      dbase->fp = stdin;
     }
   else if ((dbase->fp=fopen(dbase->filename, "r"))==NULL)
     {
      free(dbase->filename);
      free(dbase->indexname);
      free(dbase);
      return((struct Database *) NULL);
     }

   dbase->filepos = 0L;

   return(dbase);
  }

/*---------------------------------------------------------------(blastdb)---*/

int blastdb(char *name) {

   int i;

   for (i=0; i<nblastdbs; i++)
     {
      if (strcmp(name, blastdbs[i])==0) {return(TRUE);}
     }

   return(FALSE);
  }

/*------------------------------------------------------------(closedbase)---*/

extern void closedbase(struct Database *dbase) {

   fclose(dbase->fp);
   free(dbase->filename);
   free(dbase->indexname);
   free(dbase);

   return;
}

/*--------------------------------------------------------------(firstseq)---*/

extern struct Sequence *firstseq(struct Database *dbase) {

   if (dbase->filepos!=0L)
     {
      dbase->filepos = 0L;
      if (fseek(dbase->fp, dbase->filepos, 0)!=0)
        {fprintf(stderr, "Error positioning file %s for firstseq.\n",
                           dbase->filename);
         exit(1);}
     }

   return(readentry(dbase));
  }

/*---------------------------------------------------------------(nextseq)---*/

extern struct Sequence *nextseq(struct Database *dbase) {

   return(readentry(dbase));
}

/*--------------------------------------------------------------(closeseq)---*/

extern void closeseq(struct Sequence *seq) {

   if (seq==NULL) return;

   if (seq->id!=NULL)          free(seq->id);
   if (seq->name!=NULL)        free(seq->name);
   if (seq->organism!=NULL)    free(seq->organism);
   if (seq->header!=NULL)      free(seq->header);
   if (seq->state!=NULL)       free(seq->state);
   if (seq->composition!=NULL) free(seq->composition);
   free(seq->seq);
   free(seq);
   return;
  }

/*---------------------------------------------------------------(openwin)---*/

extern struct Sequence *openwin(struct Sequence *parent, int start, int length) {

  struct Sequence *win;
  int i;

   if (start<0 || length<0 || start+length>parent->length)
     {
      return((struct Sequence *) NULL);
     }

   win = (struct Sequence *) malloc(sizeof(struct Sequence));

/*---                                          ---[set links, up and down]---*/

   win->parent = parent;
   if (parent->root==NULL)
     {win->root = parent;}
   else
     {win->root = parent->root;}
   win->children = (struct Sequence **) NULL;

/* parent->children = ***foo***                   ---[not yet implemented]---*/

   win->id = (char *) NULL;
   win->name = (char *) NULL;
   win->organism = (char *) NULL;
   win->header = (char *) NULL;

/*---                          ---[install the local copy of the sequence]---*/

   win->start = start;
   win->length = length;
#if 0
   win->seq = (char *) malloc(sizeof(char)*length + 1);
   memcpy(win->seq, (parent->seq)+start, length);
   win->seq[length] = '\0';
#else
	win->seq = parent->seq + start;
#endif

/*---                          ---[setup window implementation parameters]---*/

/*---                                                 ---[set local flags]---*/

	win->rubberwin = FALSE;
	win->floatwin = FALSE;
	win->punctuation = FALSE;

/*---                                   ---[initially unconfiguerd window]---*/

	win->entropy = -2.;
	win->state = (int *) NULL;
	win->composition = (int *) NULL;
	win->classvec = (char *) NULL;
	win->scorevec = (double *) NULL;

	stateon(win);

	return win;
}

/*---------------------------------------------------------------(nextwin)---*/

extern struct Sequence *nextwin(struct Sequence *win, int shift) {

   if ((win->start+shift)<0 ||
       (win->start+win->length+shift)>win->parent->length)
     {
      return((struct Sequence *) NULL);
     }
   else
     {
      return(openwin(win->parent, win->start+shift, win->length));
     }
  }

/*--------------------------------------------------------------(shiftwin1)---*/

int shiftwin1(struct Sequence *win) {
	register int	j, length;
	register int	*comp;

	length = win->length;
	comp = win->composition;

	if ((++win->start + length) > win->parent->length) {
		--win->start;
		return FALSE;
	}

	if (!aaflag[j = win->seq[0]])
		decrementsv(win->state, comp[aaindex[j]]--);

	j = win->seq[length];
	++win->seq;

	if (!aaflag[j])
		incrementsv(win->state, comp[aaindex[j]]++);

	if (win->entropy > -2.)
		win->entropy = entropy(win->state);

	return TRUE;
}

/*--------------------------------------------------------------(closewin)---*/

extern void closewin(struct Sequence *win)  {

   if (win==NULL) return;

   if (win->state!=NULL)       free(win->state);
   if (win->composition!=NULL) free(win->composition);
   if (win->classvec!=NULL)    free(win->classvec);
   if (win->scorevec!=NULL)    free(win->scorevec);

   free(win);
   return;
  }

/*----------------------------------------------------------------(compon)---*/

extern void compon(struct Sequence *win) {
	register int	*comp;
	register int	aa;
	register char	*seq, *seqmax;

	win->composition = comp = (int *) calloc(ALPHA_SIZE*sizeof(*comp), 1);
	seq = win->seq;
	seqmax = seq + win->length;

	while (seq < seqmax) {
		aa = *seq++;
		if (!aaflag[aa])
			comp[aaindex[aa]]++;
	}

	return;
}

/*---------------------------------------------------------------(stateon)---*/

extern void stateon(struct Sequence *win) {

  int aa, i, j;

   if (win->composition==NULL) {compon(win);}

   win->state = (int *) malloc(21*sizeof(int));

   for (aa=0; aa<21; aa++)
     {
      win->state[aa] = 0;
     }

   for (aa=0; aa<20; aa++)
     {
      if (win->composition[aa]==0) {continue;}

      for (i=0; i<20; i++)
        {
         if (win->composition[aa]<win->state[i]) {continue;}

         for (j=19; j>i; j--)
           {
            win->state[j] = win->state[j-1];
           }

         win->state[i] = win->composition[aa];
         break;
        }
     }

   return;
  }

/*-----------------------------------------------------------------(enton)---*/

extern void enton(struct Sequence *win) {
   if (win->state==NULL) {stateon(win);}

   win->entropy = entropy(win->state);

   return;
  }

/*---------------------------------------------------------------(entropy)---*/
static int		thewindow;
static double	*entray;

#define LN2	0.69314718055994530941723212145818

void entropy_init(int window) {
	int		i;
	double	x, xw;

	entray = (double *)malloc((window+1) * sizeof(*entray));
	xw = window;
	for (i = 1; i <= window; ++i) {
		x = i / xw;
		entray[i] = -x * log(x) / LN2;
	}

	thewindow = window;
}

extern double entropy(int *sv) {

  int	*sv0 = sv;
  register double	ent;
  register int	i, total;
  register int	*svmax;
  register double	xtotrecip, xsv;

  for (total = 0; (i = *sv) != 0; ++sv) {
    total += i;
  }
  svmax = sv;
  ent = 0.0;
  if (total == thewindow) {
    for (sv = sv0; sv < svmax; ) {
      ent += entray[*sv++];
    }
    return ent;
  }
  if (total == 0)
    return 0.;

  xtotrecip = 1./(double)total;
  for (sv = sv0; sv < svmax; ) {
    xsv = *sv++;
    ent += xsv * log(xsv * xtotrecip);
  }
  return -ent * xtotrecip / LN2;
}

/*-----------------------------------------------------------(decrementsv)---*/

static void
decrementsv(int *sv, int class) {
  register int	svi;

  while ((svi = *sv++) != 0) {
    if (svi == class && *sv < class) {
      sv[-1] = svi - 1;
      break;
    }
  }
}

/*-----------------------------------------------------------(incrementsv)---*/

static void incrementsv(int *sv, int class) {

  for (;;) {
    if (*sv++ == class) {
      sv[-1]++;
      break;
    }
  }
}

/*-------------------------------------------------------------(readentry)---*/

struct Sequence *readentry(struct Database *dbase) {

  struct Sequence *seq;
  int	c;

  seq = seqnew();

  seq->db = dbase;

  /*---                                    ---[backpointers null at the top]---*/

  seq->parent = (struct Sequence *) NULL;
  seq->root = (struct Sequence *) NULL;
  seq->children = (struct Sequence **) NULL;

  /*---                                                       ---[set flags]---*/

  seq->rubberwin = FALSE;
  seq->floatwin = FALSE;

  /*---                                                  ---[read from file]---*/

  if (!readhdr(seq))
    {
      return((struct Sequence *) NULL);
    }
  while (1)  /*---[skip multiple headers]---*/
    {
      c = getc(dbase->fp);
      if (c == EOF)
	break;
      if (c != '>') {
	ungetc(c, dbase->fp);
	break;
      }
      while ((c=getc(dbase->fp)) != EOF && c !='\n')
	;
      if (c == EOF)
	break;
    }
  readseq(seq);

  /*---                                   ---[set implementation parameters]---*/

  /*---                                          ---[initially unconfigured]---*/

  seq->entropy = -2.;
  seq->state = (int *) NULL;
  seq->composition = (int *) NULL;
  seq->classvec = (char *) NULL;
  seq->scorevec = (double *) NULL;

  return(seq);
}

/*---------------------------------------------------------------(readhdr)---*/

int readhdr(struct Sequence *seq) {

  FILE *fp;
  char *bptr, *curpos;
  int	c, i, itotal;
  int idend, namend, orgend;

  fp = seq->db->fp;

  if ((c=getc(fp)) == EOF) {
      free(seq);
      return(FALSE);
    }
   
  while (c != EOF && isspace(c)) {
      c = getc(fp);
    }

  if (c!='>')
    {fprintf(stderr, "Error reading fasta format - '>' not found.\n");
      exit(1);}
  ungetc(c, fp);
  /*                                               ---[read the header line]---*/
  str = (struct strlist *) malloc (sizeof(struct strlist));
  str->next = NULL;
  curstr = str;

  for (i=0,itotal=0,c=getc(fp); c != EOF; c=getc(fp)) {
      if (c=='\n') break;

      if (i==STRSIZE-1)
        {curstr->string[i] = '\0';
	  curstr->next = (struct strlist *) malloc (sizeof(struct strlist));
	  curstr = curstr->next;
	  curstr->next = NULL;
	  i = 0;}

      curstr->string[i] = c;
      itotal++;
      i++;
    }

  curstr->string[i] = '\0';
  seq->header = (char *) malloc (itotal+2);
  seq->header[0] = '\0';

  for (curstr=str, curpos=seq->header; curstr!=NULL;)
    {
      if (curstr->next==NULL)
        {memccpy(curpos, curstr->string, '\0', STRSIZE);}
      else
        {memccpy(curpos, curstr->string, '\0', STRSIZE-1);}

      str = curstr;
      curstr = curstr->next;
      free (str);

      if (curstr!=NULL) {curpos = curpos+STRSIZE-1;}
    }

  bptr = (seq->header)+1;
  seq->name = (char *) NULL;
  seq->organism = (char *) NULL;
  /*                                                   ---[parse out the id]---*/
  idend = findchar(bptr, ' ');
  if (idend==-1) {idend = findchar(bptr, '\n');}
  if (idend==-1) {idend = findchar(bptr, '\0');}
  if (idend==-1)
    {fprintf(stderr, "Error parsing header line - id.\n");
      fputs(seq->header, fp);
      exit(1);}   

  seq->id = (char *) malloc((idend+1)*sizeof(char));
  memcpy(seq->id, bptr, idend);
  seq->id[idend] = '\0';

  if (bptr[idend]=='\n' || bptr[idend]=='\0') {return(TRUE);}

  /*                                         ---[parse out the protein name]---*/
  bptr = bptr + idend + 1;
  while (bptr[0]==' ') {bptr++;}

  namend = findchar(bptr, '-');
  if (namend==-1) {namend = findchar(bptr, '\n');}
  if (namend==-1) {namend = findchar(bptr, '\0');}
  if (namend==-1)
    {fprintf(stderr, "Error parsing header line - name.\n");
      fputs(seq->header, fp);
      return(TRUE);}

  seq->name = (char *) malloc((namend+1)*sizeof(char));
  memcpy(seq->name, bptr, namend);
  seq->name[namend] = '\0';

  if (bptr[namend]=='\n' || bptr[namend]=='\0') {return(TRUE);}

  /*                                                 ---[parse out organism]---*/
  bptr = bptr + namend + 1;
  while (bptr[0]==' ') {bptr++;}

  orgend = findchar(bptr, '|');
  if (orgend==-1) {orgend = findchar(bptr, '#');}
  if (orgend==-1) {orgend = findchar(bptr, '\n');}
  if (orgend==-1) {orgend = findchar(bptr, '\0');}
  if (orgend==-1)
    {fprintf(stderr, "Error parsing header line - organism.\n");
      fputs(seq->header, fp);
      return(TRUE);}

  seq->organism = (char *) malloc((orgend+1)*sizeof(char));
  memcpy(seq->organism, bptr, orgend);
  seq->organism[orgend] = '\0';

  /*                                    ---[skip over multiple header lines]---*/
  while (TRUE)
    {
      c = getc(fp);
      if (c == EOF)
	return(TRUE);
      if (c=='>')
        {
	  skipline(fp);
        }
      else
        {
	  ungetc(c,fp);
	  break;
        }
    }

  return(TRUE);
}

/*--------------------------------------------------------------(skipline)---*/

void skipline(FILE *fp) {

  int c;
  while ((c=getc(fp))!='\n' && c!=EOF)
    ;

  return;
}

/*--------------------------------------------------------------(findchar)---*/

extern int findchar(char *str, char chr) {

  int i;

  for (i=0; ; i++)
    {
      if (str[i]==chr)
        {
	  return(i);
        }
      if (str[i]=='\0')
        {
	  return(-1);
        }
    }
}

/*---------------------------------------------------------------(readseq)---*/

void readseq(struct Sequence *seq) {

  FILE *fp;
  int i, itotal;
  int	c;
  char *curpos;

  fp = seq->db->fp;

  seq->punctuation = FALSE;   

  str = (struct strlist *) malloc (sizeof(struct strlist));
  str->next = NULL;
  curstr = str;

  for (i = 0, itotal = 0, c = getc(fp); c != EOF; c = getc(fp)) {
    if (!aaflag[c]) {
    Keep:
      if (i < STRSIZE-1) {
	curstr->string[i++] = c;
	continue;
      }
      itotal += STRSIZE-1;
      curstr->string[STRSIZE-1] = '\0';
      curstr->next = (struct strlist *) malloc(sizeof(*curstr));
      curstr = curstr->next;
      curstr->next = NULL;
      curstr->string[0] = c;
      i = 1;
      continue;
    }

    switch (c) {
    case '>':
      ungetc(c, fp);
      goto EndLoop;
    case '*': case '-':
      seq->punctuation = TRUE;
      goto Keep;
    case 'x': case 'X':
      goto Keep;

    }

    if (isalpha(c)) {
      c = 'X';
      goto Keep;
    } else continue;
  }
 EndLoop:
  itotal += i;

  curstr->string[i] = '\0';
  seq->seq = (char *) malloc (itotal+2);
  seq->seq[0] = '\0';

  for (curstr = str, curpos = seq->seq; curstr != NULL;) {
    if (curstr->next == NULL)
      memccpy(curpos, curstr->string, '\0', STRSIZE);
    else
      memccpy(curpos, curstr->string, '\0', STRSIZE-1);

    str = curstr;
    curstr = curstr->next;
    free(str);

    if (curstr != NULL)
      curpos = curpos+STRSIZE-1;
  }

  seq->length = strlen(seq->seq);

  return;
}
/*-----------------------------------------------------------------(upper)---*/

extern void upper(char *string, size_t len)
{
  register char	*stringmax, c;

  for (stringmax = string + len; string < stringmax; ++string)
    if (islower(c = *string))
      *string = toupper(c);
}

/*-----------------------------------------------------------------(lower)---*/

extern void lower(char *string, size_t len)
{
  register char	*stringmax, c;

  for (stringmax = string + len; string < stringmax; ++string) {
    if (isupper(c = *string)) {
      *string = tolower(c);
    }
  }
}

/*-------------------------------------------------------------------(min)---*/

int min(int a, int b) {

   if (a<b) {return(a);}
   else {return(b);}
}

/*-------------------------------------------------------------------(max)---*/

int max(int a, int b) {

  if (a<b) {return(b);}
  else {return(a);}

}

/*---------------------------------------------------------------------------*/

/*---------------------------------------------------------------(tmalloc)---*/

void *tmalloc(size_t size) {

  void *ptr;

  ptr = (void *) malloc(size);

  if (rptr>TESTMAX)
    {
      exit(2);
    }

  record_ptrs[rptr] = (size_t) ptr;
  rptr++;

  return(ptr);
}

/*-----------------------------------------------------------------(tfree)---*/

void tfree(void *ptr) {

  int i;

  for (i=0; i<rptr; i++) {
    if (record_ptrs[i]==(size_t)ptr) {
      record_ptrs[i] = 0;
      break;
    }
  }

  free(ptr);
}

/*---------------------------------------------------------------(seqnew)----*/

struct Sequence *seqnew() {

  struct Sequence *seq;

  seq = (struct Sequence *) malloc (sizeof (struct Sequence));

  seq->db = (struct Database *) NULL;
  seq->id = (char *) NULL;
  seq->name = (char *) NULL;
  seq->organism = (char *) NULL;
  seq->header = (char *) NULL;
  seq->seq = (char *) NULL;
  seq->start = 0;
  seq->length = 0;
  seq->alphabet = (struct Alphabet *) NULL;
  seq->parent = (struct Sequence *) NULL;
  seq->root = (struct Sequence *) NULL;
  seq->children = (struct Sequence **) NULL;
  seq->bogus = FALSE;
  seq->punctuation = FALSE;
  seq->rubberwin = FALSE;
  seq->floatwin = FALSE;
  seq->seedwin = FALSE;
  seq->state = (int *) NULL;
  seq->entropy = 0.0;
  seq->composition = (int *) NULL;
  seq->comptst = FALSE;
  seq->compindex = (int *) NULL;
  seq->classvec = (char *) NULL;
  seq->clalphabet = (struct Alphabet *) NULL;
  seq->scorevec = (double *) NULL;
  seq->score = 0.0;
  seq->foo = (struct PerScoreVec *) NULL;
  seq->bar = 0;

  return(seq);
}

/*---------------------------------------------------------------------------*/
