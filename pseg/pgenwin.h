
/*****************************************************************************/
/***   (genwin.h)                                                          ***/
/*****************************************************************************/

/*--------------------------------------------------------------(includes)---*/

#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <string.h>
#include <fcntl.h>
#include <ctype.h>
#include <math.h>

/*--------------------------------------------------------------(typedefs)---*/

typedef int Bool;

/*---------------------------------------------------------------(defines)---*/
#ifndef TRUE
#define TRUE 1
#define FALSE 0
#endif

#define ALPHA_SIZE 20
#define ALPHA_SIZE_PLUS 21

/*---------------------------------------------------------------(structs)---*/

struct Database
  {
   char *filename;
   FILE *fp;
   char *indexname;
   long int filepos;
  };

struct Sequence
  {
   struct Database *db;

   char *id;
   char *name;
   char *organism;
   char *header;

   char *seq;
   int start;                     /* for windows */
   int length;
   struct Alphabet *alphabet;

   struct Sequence *parent;       /* for windows */
   struct Sequence *root;         /* for subwindows */
   struct Sequence **children;    /* only the floaters? */

   Bool bogus;
   Bool punctuation;
   Bool rubberwin;             /* for windows */
   Bool floatwin;              /* for subwindows */
   Bool seedwin;

   int *state;
   double entropy;
   int *composition;
   Bool	*comptst;
   int	*compindex;

   char *classvec;               /* from ClaVec[aa] */
   struct Alphabet *clalphabet;
   double *scorevec;             /* from ScoVec[aa] or ScoFun(pos) */
   double score;                 /* from ScoFun(win) */
                          /* union, for integer scores? */
   struct PerScoreVec *foo;
   int *bar;
  };

struct Configuration
  {
   char *iseq;
   int ilength;

   int printper;
  };

struct Matrix
  {
   struct Sequence *parent;

   char **seq;
   int start;
   int period;
   int length;

   int **rowcomp;
   int **colcomp;
  };

/*------------------------------------------------------(alphabet structs)---*/

struct Alphabet
  {
   char *name;
   int size;

   int *index[128];
   char *chars;       /*  [size]  */
   char **charnames;  /*  [size]  */

   Bool caseinvariant;
  };

struct TransVector
  {
   struct Alphabet *from;
   struct Alphabet *to;

   int *index;     /*  [from->size]  */
  };

struct ScoreVector
  {
   struct Alphabet *from;

   double *score;   /*  [from->size]  */
   char *units;
  };

struct ScoreMatrix
  {
   struct Alphabet *from1;
   struct Alphabet *from2;

   double **score;  /*  [from1->size][from2->size]  */
   char *units;
  };

/*---------------------------------------------------------(bogus structs)---*/

struct PerScoreVec
  {
   int hits;
   int tot;

   double pct;
   double std;
   double prob;
  };

/*----------------------------------------------------------------(protos)---*/

extern struct Database *opendbase(char *);
extern void closedbase(struct Database *);

extern struct Sequence *openseq(struct Sequence *, int, int); 
extern struct Sequence *firstseq(struct Database *);
extern struct Sequence *nextseq(struct Database *);
extern struct Sequence *seqnew();

extern void closeseq(struct Sequence *);

extern void  genwininit();
extern struct Sequence *openwin(struct Sequence *parent, int start, int length);
extern struct Sequence *nextwin(struct Sequence *, int);
extern void  closewin(struct Sequence *);

extern int shiftwin(struct Sequence *);
extern int shiftwin1(struct Sequence *);

extern void compon(struct Sequence *), stateon(struct Sequence *), enton(struct Sequence *);
extern double entropy(int *);

/* extern struct Matrix *openmat(); */
/* extern closemat(); */

extern void upper(char *, size_t), lower(char *, size_t);
int findchar(char *, char);

/*----------------------------------------------------------------(macros)---*/

/***   #define bogus(A) (aaindex[A]>=ALPHA_SIZE)   ***/

/*---------------------------------------------------------------(globals)---*/

extern int aaindex[];
extern char aachar[];

/*---------------------------------------------------------------------------*/
