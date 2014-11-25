/* Concurrent read version */
/*	ascii.gbl	ascii translation to amino acids */
/*	modified 10-Mar-1987 for B, Z	*/

/* $Id: uascii.h 989 2012-07-24 19:37:38Z wrp $ */

#define NA 123
#define NANN 60	/* changed 24-July-2012 because NCBIstdaa_ext_n = 56 */
#define ESS 59	/* code for ',' in FASTS,FASTF, FASTM */
#define EL 125
#define ES 126
#define AAMASK 127

#ifndef XTERNAL
/*       0  1  2  3  4  5  6  7  8  9  10 11 12 13 14 15	*/
/* 32	    !  "  #  $  %  &  '  (  )  *  +  ,  -  .  / 	*/
/* 48	 0  1  2  3  4  5  6  7  8  9  :  ;  <  =  >  ? 	*/
/* 64	 @  A  B  C  D  E  F  G  H  I  J  K  L  M  N  O		*/ 
/* 80	 P  Q  R  S  T  U  V  W  X  Y  Z  [  \  ]  ^  _		*/ 
/* 96	 `  a  b  c  d  e  f  g  h  i  j  k  l  m  n  o 	*/
/*112	 p  q  r  s  t  u  v  w  x  y  z  {  |  }  ~  ^?	*/ 

int aascii[128]={
	EL,NA,NA,NA,NA,NA,NA,NA,NA,NA,EL,NA,NA,EL,NA,NA,	/* 15 */
	NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,	/* 31 */
	NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,24,NA,NA,NA,NA,NA,	/* 47 */
	NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,	/* 63 */
	NA, 1,21, 5, 4, 7,14, 8, 9,10,25,12,11,13, 3,26,	/* 79 */
	15, 6, 2,16,17,27,20,18,23,19,22,NA,NA,NA,NA,NA,	/* 95 */
	NA, 1,21, 5, 4, 7,14, 8, 9,10,25,12,11,13, 3,26,	/*111 */
	15, 6, 2,16,17,27,20,18,23,19,22,NA,NA,NA,NA,NA};	/*127 */

int nascii[128]={
/*	 0  1  2  3  5  6  7  8  9 10 11 12 13 14 15 15
	 @  A  B  C  D  E  F  G  H  I  J  K  L  M  N  O
	 P  Q  R  S  T  U  V  W  X  Y  Z		*/
	EL,NA,NA,NA,NA,NA,NA,NA,NA,NA,EL,NA,NA,EL,NA,NA,
	NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,
	NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,ES,NA,NA,16,NA,NA,
	NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,ES,NA,NA,ES,NA,
	NA, 1,15, 2,12,NA,NA, 3,13,NA,NA,11,NA, 8,16,NA,
	 6, 7, 6,10, 4, 5,14, 9,17, 7,NA,NA,NA,NA,NA,NA,
	NA, 1,15, 2,12,NA,NA, 3,13,NA,NA,11,NA, 8,16,NA,
	 6, 7, 6,10, 4, 5,14, 9,17, 7,NA,NA,NA,NA,NA,NA};

int *pascii;
int qascii[128];
int lascii[128];
int l_ann_ascii[128];
#else
#define AAMASK 127
extern int aascii[128];
extern int nascii[128];

extern int *pascii;
extern int qascii[128];
extern int lascii[128];
extern int l_ann_ascii[128];
#endif
