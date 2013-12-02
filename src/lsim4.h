
/* $Id: lsim4.h 625 2011-03-23 17:21:38Z wrp $ */
/* $Revision: 625 $  */

/* global definitions, #defines, for lsim4.c */

typedef int bool;

void *ckalloc(size_t amount, size_t size);

#define maxi(x, y)  (((x) > (y)) ? x: y)

typedef struct ONE {int COL ; struct ONE  *NEXT ;} pair, *pair_ptr;

#define PAIRNULL (pair_ptr)NULL

typedef struct NODE
{ int  SCORE;
  int  STARI, STARJ;
  int  ENDI,  ENDJ;
  int  TOP,   BOT;
  int  LEFT,  RIGHT; 
  struct NODE *next;
}  vertex, *vertex_p;
		
struct lrr_str {
  int CC, DD;			/* saving row matrix scores */
  int RR, SS, EE, FF;		/* saving row start-points */
};
  
struct lcc_str {
  int HH, WW;			/* saving col matrix scores HH=CC, */
  int II, JJ, XX, YY; 		/* saving col start-points , II=RR, JJ=EE */
};

typedef struct spa {
  int CC, RR, EE, DD, SS, FF; 
} space;

typedef struct spa *space_ptr;

struct vert_str {
  int numnode;
  vertex_p LIST, most;
};

struct l_struct {
  space_ptr CCC;
  struct lrr_str *r_ss;
  struct lcc_str *c_ss;
  pair_ptr *row;		/* for saving used aligned pairs */
  int  m1, mm, n1, nn;	/* boundaries of recomputed area */
  int  rl, cl;		/* left and top boundaries */
  int I, J;		/* current positions of A ,B - used by diff() */
};

static void big_pass(const unsigned char *A,
		     const unsigned char *B,
		     int M, int N,
		     int mini_score,
		     int **pam2, int Q, int R,
		     int nseq,
		     struct vert_str *v_ptr,
		     struct l_struct *l_ptr);

static void locate(const unsigned char *A,
		   const unsigned char *B,
		   int mini_score,
		   int **pam2, int Q, int R,
		   int nseq,
		   int *flag_p,
		   struct vert_str *v_ptr,
		   struct l_struct *l_ptr);

static void small_pass(const unsigned char *A,
		       const unsigned char *B,
		       int mini_score,
		       int **pam2, int Q, int R,
		       int nseq,
		       struct vert_str *v_ptr,
		       struct l_struct *l_ptr);

static void addnode(int c, int ci, int cj, int i, int j,
		    struct vert_str *v_ptr);

static bool no_cross(int *flag_p, vertex_p LIST, struct l_struct *l_ptr);

static int diff(const unsigned char *A,
		const unsigned char *B,
		int M, int N, int tb, int te,
		int two_seq,
		int **pam2, int q, int r,
		int **sapp, int *last,
		struct l_struct *l_ptr);

static int CHECK_SCORE(const unsigned char *A, const unsigned char *B,
		       int M, int N,
		       int *S, int **W, int G, int H, int *nres);

static vertex_p findmax(struct vert_str *v_ptr);

/* DIAG() assigns value to x if (ii,jj) is never used before */
#define DIAG(ii, jj, x, value)	\
{ for ( z = l_ptr->row[(ii)]; z != 0 && z->COL != (jj); z = z->NEXT ) ; \
  if ( !z )   x = ( value );	\
  }

/* replace (ss1, xx1, yy1) by (ss2, xx2, yy2) if the latter is large */
#define ORDER(ss1, xx1, yy1, ss2, xx2, yy2) \
{ if ( ss1 < ss2 ) \
  { ss1 = ss2; xx1 = xx2; yy1 = yy2; } \
else \
if ( ss1 == ss2 ) \
  { if ( xx1 < xx2 ) { xx1 = xx2; yy1 = yy2; } \
    else \
       if ( xx1 == xx2 && yy1 < yy2 )  yy1 = yy2; \
  } \
}

#define ORDER1(ss1, xx1, yy1, ss2, xx2, yy2)  \
{ if (ss1 <= ss2) {  \
    if (ss1 == ss2) {  \
	if (xx1 < xx2) {  \
           xx1 = xx2; yy1 = yy2; \
    } else {                           \
      if (xx1 == xx2 && yy1 < yy2)  \
         yy1 = yy2; \
    } \
    } else {  \
      ss1 = ss2; xx1 = xx2; yy1 = yy2;  \
    } \
  } \
}

#define ORDER2(ss1, xx1, ss2, xx2) \
{  \
   if (ss1 <= ss2) { \
      if (ss1 == ss2) {  \
         if (xx1 < xx2) xx1 = xx2; \
      } else { \
      ss1 = ss2; xx1 = xx2; \
      } \
 } \
}

