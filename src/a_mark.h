/* a_mark.h - symbols used to indicate match/mismatch alignment code */

/* copyright (c) 2003 William R. Pearson and the U. of Virginia */

/* $Id: a_mark.h 1024 2012-08-07 18:08:45Z wrp $ */
/* $Revision: 1024 $  */

/* character types in aln_code strings */
#define M_BLANK 0
#define M_NEG 1
#define M_ZERO 2
#define M_POS 3
#define M_IDENT 4
#define M_DEL 5

#define MX_A0 0
#define MX_A1 1
#define MX_A2 2
#ifdef M10_CONS_L
#define MX_A10 3
#else
#define MX_A10 4
#endif
#define MX_ACC 5
#define MX_ABLAST 6

static char *
aln_map_sym[] = {"  ..: ",	/* 0 */
		 " Xxx  ",	/* 1 */
		 "    . ",	/* 2 */
		 " mzp=-",	/* 3: 10a */
		 "  ..:-",	/* 4: 10b */
		 " <z>=-",	/* 5: calc_code */
		 "   += "	/* 6: MX_MBLAST blast */
		  };



