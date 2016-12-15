/*	initfa.c	*/
/*  $Id: initfa.c 1274 2014-08-07 18:30:56Z wrp $ */

/* copyright (c) 1996, 1997, 1998, 2014 by William R. Pearson and The
   Rector & Visitors of the University of Virginia */

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

/* init??.c files provide function specific initializations */

/* h_init()	- called from comp_lib.c, comp_thr.c to initialize pstruct ppst
   		  which includes the alphabet, and pam matrix

   alloc_pam()	- allocate pam matrix space
   init_pam2()	- convert from 1D to 2D pam

   init_pamx()	- extend pam2 for 'X', 'N', or lower case characters

   f_initenv()	- set up mngmsg and pstruct defaults
   f_getopt()	- read fasta specific command line options
   f_getarg()	- read ktup

   resetp()	- reset the parameters, scoring matrix for DNA-DNA/DNA-prot

   query_parm()	- ask for ktup
   last_init()	- some things must be done last

   f_initpam()	- set some parameters based on the pam matrix
*/

#include <stdio.h>
#include <stdlib.h>
#include <ctype.h>
#include <string.h>
#include <math.h>

#ifdef UNIX
#include <sys/types.h>
#include <sys/stat.h>
#endif

#include "defs.h"
#include "structs.h"
#include "param.h"
#include "best_stats.h"

#define XTERNAL
#include "upam.h"
#include "uascii.h"
#undef XTERNAL

#define MAXWINDOW 32

int initpam(char *, struct pstruct *);
void init_pam2 (struct pstruct *ppst);
void init_altpam(struct pstruct *ppst);
void init_pamx (struct pstruct *ppst);
void extend_pssm(unsigned char *aa0, int n0, struct pstruct *ppst);
void build_xascii(int *qascii, char *save_str);
void add_ascii_ann(int *qascii, unsigned char *ann_arr);
void re_ascii(int *qascii, int *pascii, int max_ann_arr);
extern int my_nrand(int, void *);

/*  at some point, all the defaults should be driven from this table */
/*
#pgm	q_seq	l_seq	p_seq	matrix	g_open	g_ext	fr_shft	e_cut	ktup  E_band_opt
#	-n/-p		-s	-e	-f	-h/-j	-E	argv[3]
fasta	prot(0)	prot(0)	prot(0)	bl50	-10	-2	-	10.0	2	0.02
fasta	dna(1)	dna(1)	dna(1)	+5/-4	-12	-4	-	2.0	6	0.01
ssearch	prot(0)	prot(0)	prot(0)	bl50	-10	-2	-	10.0	-	-
ssearch	dna(1)	dna(1)	dna(1)	+5/-4	-16	-4	-	2.0	-	-
fastx	dna(1)	prot(0)	prot(0)	BL50	-12	-2	-20	5.0	2	0.02
fasty	dna(1)	prot(0)	prot(0)	BL50	-12	-2	-20/-24	5.0	2	0.02
tfastx	dna(1)	prot(0)	prot(0)	BL50	-14	-2	-20	5.0	2	0.01
tfasty	dna(1)	prot(0)	prot(0)	BL50	-14	-2	-20/-24	5.0	2	0.01
fasts	prot(0)	prot(0)	prot(0)	MD20-MS	-	-	-	5.0	-	-
fasts	dna(1)	dna(1)	dna(1)	+2/-4	-	-	-	5.0	1	-
tfasts	prot(0)	dna(1)	prot(0)	MD10-MS	-	-	-	2.0	1	-
fastf	prot(0)	prot(0)	prot(0)	MD20	-	-	-	2.0	1	-
tfastf	prot(0)	dna(1)	prot(0)	MD10	-	-	-	1.0	1	-
fastm	prot(0)	prot(0)	prot(0)	MD20	-	-	-	5.0	1	-
fastm	dna(1)	dna(1)	dna(1)	+2/-4	-	-	-	2.0	1	-
tfastm	prot(0)	dna(1)	prot(0)	MD10	-	-	-	2.0	1	-
lalign	prot(0) prot(0) prot(0) BL50    -12     -2	-       1.0	-	-
lalign	dna(1)  dna(1)  dna(1)  +5/-4   -12     -4	-	0.1	-	-
*/

void show_help(char *, int );

char *ref_str_a[]={
/* 0 */ "W.R. Pearson & D.J. Lipman PNAS (1988) 85:2444-2448\n",
/* 1 */ "T. F. Smith and M. S. Waterman, (1981) J. Mol. Biol. 147:195-197; \n W.R. Pearson (1991) Genomics 11:635-650\n",
/* 2 */ "Pearson et al, Genomics (1997) 46:24-36\n",
/* 3 */ "Mackey et al. Mol. Cell. Proteomics  (2002) 1:139-147\n",
/* 4 */ "W.R. Pearson (1996) Meth. Enzymol. 266:227-258\n",
/* 5 */ "X. Huang and W. Miller (1991) Adv. Appl. Math. 12:373-381\n",
  ""
};

#define FA_PID	1
#define SS_PID	2
#define FX_PID	3
#define FY_PID	4
#define FS_PID	5
#define FF_PID	6
#define FM_PID	7
#define RSS_PID	8
#define RFX_PID	9
#define SSS_PID 10	/* old (slow) non-PG Smith-Waterman */
#define TFA_PID	FA_PID+10
#define TFX_PID	FX_PID+10
#define TFY_PID	FY_PID+10
#define TFS_PID	FS_PID+10
#define TFF_PID	FF_PID+10
#define TFM_PID FM_PID+10
#define LAL_PID 18
#define LNW_PID 19
#define GNW_PID 20

struct pgm_def_str {
  int pgm_id;
  char *prog_func;
  char *info_pgm_abbr;
  char *iprompt0;
  char *ref_str;
  int PgmDID;
  char *smstr;
  int g_open_mod;
  int g_ext_mod;
  int gshift;
  int hshift;
  double e_cut;
  int ktup;
  double E_band_opt;
  int can_pre_align;
};

static struct pgm_def_str
pgm_def_arr[21] = {
  {0, "", "", "", NULL, 400, "", 0, 0, 0, 0, 1.0, 0, 0 },  /* 0 */
  {FA_PID, "FASTA", "fa",
   "FASTA searches a protein or DNA sequence data bank",
   NULL, 401, "BL50", 0, 0, 0, 0, 10.0, 2, 0.2, 1}, /* 1 - FASTA */
  {SS_PID, "SSEARCH","gsw","SSEARCH performs a Smith-Waterman search",
   NULL, 404, "BL50", 0, 0, 0, 0, 10.0, 0, 0.0, 1}, /* 2 - SSEARCH */
  {FX_PID, "FASTX","fx",
   "FASTX compares a DNA sequence to a protein sequence data bank",
   NULL, 405, "BL50", -2, 0, -20, 0, 5.0, 2, 0.10, 1}, /* 3 - FASTX */
  {FY_PID, "FASTY", "fy",
   "FASTY compares a DNA sequence to a protein sequence data bank",
   NULL, 405, "BL50", -2, 0, -20, -24, 5.0, 2, 0.10, 1}, /* 4 - FASTY */
  {FS_PID, "FASTS", "fs",
   "FASTS compares linked peptides to a protein data bank",
   NULL, 400, "MD20-MS", 0, 0, 0, 0, 5.0, 1, 0.0, 0}, /* 5 - FASTS */
  {FF_PID, "FASTF", "ff",
   "FASTF compares mixed peptides to a protein databank",
   NULL, 400, "MD20", 0, 0, 0, 0, 2.0, 1, 0.0, 0 }, /* 6 - FASTF */
  {FM_PID, "FASTM", "fm",
   "FASTM compares ordered peptides to a protein data bank",
   NULL, 400, "MD20", 0, 0, 0, 0, 5.0, 1, 0.0, 0 }, /* 7 - FASTM */
  {RSS_PID, "PRSS", "rss",
   "PRSS evaluates statistical signficance using Smith-Waterman",
   NULL, 401, "BL50", 0, 0, 0, 0, 1000.0, 0, 0.0, 1 }, /* 8 - PRSS */
  {RFX_PID,"PRFX", "rfx",
   "PRFX evaluates statistical signficance using FASTX",
   NULL, 401, "BL50", -2, 0, -20, -24, 1000.0, 2, 0.2, 1 }, /* 9 - PRFX */
  {SSS_PID, "OSEARCH","ssw","OSEARCH searches a sequence data bank",
   NULL, 404, "BL50", 0, 0, 0, 0, 10.0, 0, 0.0, 1}, /* 2 - OSEARCH */
  {TFA_PID, "TFASTA", "tfa",
   "TFASTA compares a protein  to a translated DNA data bank",
   NULL, 402, "BL50", -2, 0, 0, 0, 5.0, 2, 0.1, 1},
  {0, "", "", "", NULL, 400, "", 0, 0, 0, 0, 1.0, 0, 0.0 },  /* 0 */
  {TFX_PID, "TFASTX", "tfx",
   "TFASTX compares a protein to a translated DNA data bank",
   NULL, 406, "BL50", -2, 0, -20, 0, 2.0, 2, 0.10, 1},
  {TFY_PID, "TFASTY", "tfy",
   "TFASTY compares a protein to a translated DNA data bank",
   NULL, 406, "BL50", -2, 0, -20, -24, 2.0, 2, 0.10, 1},
  {TFS_PID, "TFASTS", "tfs",
   "TFASTS compares linked peptides to a translated DNA data bank",
   NULL, 400, "MD10-MS", 0, 0, 0, 0, 2.0, 2, 0.0, 0 },
  {TFF_PID, "TFASTF", "tff",
   "TFASTF compares mixed peptides to a protein databank",
   NULL, 400, "MD10", 0, 0, 0, 0, 1.0, 1, 0.0, 0 },
  {TFM_PID, "TFASTM", "tfm",
   "TFASTM compares ordered peptides to a translated DNA databank",
   NULL, 400, "MD10", 0, 0, 0, 0, 1.0, 1, 0.0, 0 },
  {LAL_PID, "LALIGN", "lal",
   "LALIGN finds non-overlapping local alignments",
   NULL, 404, "BL50", -2, 0, 0, 0, 1.0, 0, 0.0, 1}, 	/* 18 - LALIGN */
  {LNW_PID, "GLSEARCH", "lnw",
   "GLSEARCH performs a global-query/local-library search",
   NULL, 404, "BL50", -2, 0, 0, 0, 10.0, 0, 0.0, 1}, 	/* 19 - GLSEARCH */
  {GNW_PID, "GGSEARCH", "gnw",
   "GGSEARCH performs a global/global database searches",
   NULL, 404, "BL50", 0, 0, 0, 0, 10.0, 0, 0.0, 1}, 	/* 20 - GGSEARCH */
};

struct msg_def_str {
  int pgm_id;
  int q_seqt;
  int l_seqt;
  int p_seqt;
  int sw_flag;
  int stages;
  int qframe;
  int nframe;
  int nrelv, srelv, arelv;
  char *f_id0, *f_id1, *label, *alabel;
};

/* align_label must be < MAX_SSTR (32) */
char *align_label[]={
  "Smith-Waterman",		/* 0 */
  "banded Smith-Waterman",	/* 1 */
  "Waterman-Eggert",		/* 2 */
  "trans. Smith-Waterman",	/* 3 */
  "global/local",		/* 4 */
  "trans. global/local",	/* 5 */
  "global/global (N-W)"	/* 6 */
};

/* pgm_id    q_seqt     l_seqt   p_seqt sw_f st qf nf nrv srv arv s_ix */
static struct msg_def_str
msg_def_arr[21] = {
  {0, 0, 0,  0, 0, 0, 0, 0, 0, 0, 0, "", "", ""},	/* ID=0 */
  {FA_PID, SEQT_UNK, SEQT_PROT, SEQT_PROT, 1, 1, 1, -1, 3, 1, 3,
   "fa","sw", "opt"},
  {SS_PID, SEQT_UNK, SEQT_PROT, SEQT_PROT, 1, 1, 1, -1, 1, 1, 1,
   "sw","sw", "s-w"},
  {FX_PID, SEQT_DNA, SEQT_PROT, SEQT_PROT, 1, 1, 2, -1, 3, 1, 3,
   "fx","sx", "opt"},
  {FY_PID, SEQT_DNA, SEQT_PROT, SEQT_PROT, 1, 1, 2, -1, 3, 1, 3,
   "fy","sy", "opt"},
  {FS_PID, SEQT_UNK, SEQT_PROT, SEQT_PROT, 1, 1, 1, -1, 3, 2, 3,
   "fs","fs", "initn init1"},
  {FF_PID, SEQT_PROT,SEQT_PROT, SEQT_PROT, 1, 1, 1, -1, 3, 2, 3,
   "ff","ff", "initn init1"},
  {FM_PID, SEQT_UNK,SEQT_PROT, SEQT_PROT, 1, 1, 1, -1, 3, 2, 3,
   "fm","fm","initn init1"},
  {RSS_PID, SEQT_UNK,SEQT_PROT, SEQT_PROT, 0, 1, 1, -1, 1, 1, 1,
   "rss","sw","s-w"},
  {RFX_PID, SEQT_DNA,SEQT_PROT, SEQT_PROT, 0, 1, 2, -1, 3, 1, 3,
   "rfx","sx","opt"},
  {SSS_PID, SEQT_UNK,SEQT_PROT, SEQT_PROT, 1, 1, 1, -1, 1, 1, 1,
   "sw","sw", "s-w"},
  {TFA_PID, SEQT_PROT,SEQT_DNA, SEQT_PROT, 0, 1, 1, 6, 3, 1, 3,
   "tfa","fa","initn init1"},
  {0, 0, 0,  0, 0, 0, 0, 0, 0, 0, 0, "", "", ""},	/* ID=12 */
  {TFX_PID, SEQT_PROT,SEQT_DNA, SEQT_PROT, 1, 1, 1, 2, 3, 2, 3,
   "tfx","sx","initn opt"},
  {TFY_PID, SEQT_PROT,SEQT_DNA, SEQT_PROT, 1, 1, 1, 2, 3, 2, 3,
   "tfy","sy","initn opt"},
  {TFS_PID, SEQT_PROT,SEQT_DNA, SEQT_PROT, 1, 1, 1, 6, 3, 2, 3,
   "tfs","fs","initn init1"},
  {TFF_PID, SEQT_PROT,SEQT_DNA, SEQT_PROT, 1, 1, 1, 6, 3, 2, 3,
   "tff","ff","initn init1"},
  {TFM_PID, SEQT_PROT,SEQT_DNA, SEQT_PROT, 1, 1, 1, 6, 3, 2, 3,
   "tfm","fm","initn init1"},
  {LAL_PID, SEQT_UNK, SEQT_PROT, SEQT_PROT, 1, 1, 1, -1, 1, 1, 1,
   "lsw","lsw", "ls-w"},
  {LNW_PID, SEQT_UNK, SEQT_PROT, SEQT_PROT, 1, 1, 1, -1, 1, 1, 1,
   "gnw","gnw", "n-w"},
  {GNW_PID, SEQT_UNK, SEQT_PROT, SEQT_PROT, 1, 1, 1, -1, 1, 1, 1,
   "gnw","gnw", "n-w"},
};

int
get_pgm_id() {

  int rval=0;

#ifdef FASTA
#ifndef TFAST
  pgm_def_arr[FA_PID].ref_str = ref_str_a[0];
  msg_def_arr[FA_PID].alabel = align_label[0];
  rval=FA_PID;
#else
  pgm_def_arr[TFA_PID].ref_str = ref_str_a[0];
  msg_def_arr[TFA_PID].alabel = align_label[2];
  rval=TFA_PID;
#endif
#endif

#ifdef FASTX
#ifndef TFAST
#ifndef PRSS
  pgm_def_arr[FX_PID].ref_str = ref_str_a[2];
  msg_def_arr[FX_PID].alabel = align_label[3];
  rval=FX_PID;
#else
  pgm_def_arr[RFX_PID].ref_str = ref_str_a[2];
  msg_def_arr[FX_PID].alabel = align_label[3];
  rval=RFX_PID;
#endif
#else
  pgm_def_arr[TFX_PID].ref_str = ref_str_a[2];
  msg_def_arr[TFX_PID].alabel = align_label[3];
  rval=TFX_PID;
#endif
#endif

#ifdef FASTY
#ifndef TFAST
  pgm_def_arr[FY_PID].ref_str = ref_str_a[2];
  msg_def_arr[FY_PID].alabel = align_label[3];
  rval=FY_PID;
#else
  pgm_def_arr[TFY_PID].ref_str = ref_str_a[2];
  msg_def_arr[TFY_PID].alabel = align_label[3];
  rval=TFY_PID;
#endif
#endif

#ifdef FASTS
#ifndef TFAST
  pgm_def_arr[FS_PID].ref_str = ref_str_a[3];
  msg_def_arr[FS_PID].alabel = align_label[4];
  rval=FS_PID;
#else
  pgm_def_arr[TFS_PID].ref_str = ref_str_a[3];
  msg_def_arr[TFS_PID].alabel = align_label[5];
  rval=TFS_PID;
#endif
#endif

#ifdef FASTF
#ifndef TFAST
  pgm_def_arr[FF_PID].ref_str = ref_str_a[3];
  msg_def_arr[FF_PID].alabel = align_label[4];
  rval=FF_PID;
#else
  pgm_def_arr[TFF_PID].ref_str = ref_str_a[3];
  msg_def_arr[TFF_PID].alabel = align_label[5];
  rval=TFF_PID;
#endif
#endif

#ifdef FASTM
#ifndef TFAST
  pgm_def_arr[FM_PID].ref_str = ref_str_a[3];
  msg_def_arr[FM_PID].alabel = align_label[4];
  rval=FM_PID;
#else
  pgm_def_arr[TFM_PID].ref_str = ref_str_a[3];
  msg_def_arr[TFM_PID].alabel = align_label[5];
  rval=TFM_PID;
#endif
#endif

#ifdef SSEARCH
#define CAN_PSSM
  pgm_def_arr[SS_PID].ref_str = ref_str_a[1];
  msg_def_arr[SS_PID].alabel = align_label[0];
  rval=SS_PID;
#endif

#ifdef OSEARCH
  pgm_def_arr[SSS_PID].ref_str = ref_str_a[1];
  msg_def_arr[SSS_PID].alabel = align_label[0];
  rval=SSS_PID;
#endif

#ifdef LALIGN
#define CAN_PSSM
  pgm_def_arr[LAL_PID].ref_str = ref_str_a[5];
  msg_def_arr[LAL_PID].alabel = align_label[2];
  rval=LAL_PID;
#endif

#ifdef GLSEARCH
#define CAN_PSSM
  pgm_def_arr[LNW_PID].ref_str = ref_str_a[6];
  msg_def_arr[LNW_PID].alabel = align_label[4];
  rval=LNW_PID;
#endif

#ifdef GGSEARCH
#define CAN_PSSM
  pgm_def_arr[GNW_PID].ref_str = ref_str_a[6];
  msg_def_arr[GNW_PID].alabel = align_label[6];
  rval=GNW_PID;
#endif

  return rval;
}

extern struct opt_def_str g_options[];
extern void set_opt_disp_defs(char opt_char, struct opt_def_str *options, int type,
			      int i_param1, int i_param2,
			      double d_param1, double d_param2, char *s_param);

static char z_opt_descr[] = "Statistics estimation method:\n      1 - regression; -1 - no stats.; 0 - no scaling; 2 - Maximum Likelihood Est.;\n      3 - Altschul/Gish; 4 - iter. regress.; 5 - regress w/variance;\n      6 - MLE with comp. adj.;\n     11 - 16 - estimates from shuffled library sequences;\n     21 - 26 - E2()-stats from shuffled high-scoring sequences;";

static char s_opt_descr[] = "Scoring matrix: (protein)\n      BL50, BP62 (sets -f -11 -g -1); P250, OPT5, VT200,\n      VT160, P120, VT120, BL80, VT80, MD40, VT40, MD20, VT20, MD10, VT10;\n      scoring matrix file name; -s ?BL50 adjusts matrix for short queries;";


struct opt_def_str f_options[] = {
  {'3', 0, "norevcomp", "compare forward strand only", NULL, 0, 0, 0, 0, 0.0, 0.0, NULL},
#if defined(FASTA) || defined(SSEARCH)
  {'a', 0, "show_all", "show complete Query/Sbjct sequences in alignment", NULL, 0, 0, 0, 0, 0.0, 0.0, NULL},
  {'W', 1, "context", "alignment context length (surrounding unaligned sequence)", NULL, 0, 0, 0, 0, 0.0, 0.0, NULL},
#endif
#if defined(FASTA)
  {'A', 0, "sw_align", "Smith-Waterman for final DNA alignment, band alignment for protein\n      default is band-alignment for DNA, Smith-Waterman for protein", NULL, 0, 0, 0, 0, 0.0, 0.0, NULL},
#endif
  {'b', 1, "num_descriptions", "high scores reported (limited by -E by default)", 
   "high scores reported (limited by -E by default);\n      =<int> forces <int> results;", 0, 0, 0, 0, 0.0, 0.0, NULL},
  {'d', 1, "num_alignments", "number of alignments shown (limited by -E by default)", NULL, 0, 0, 0, 0, 0.0, 0.0, NULL},
#if defined(FASTA) || defined(FASTX) || defined(FASTY)
  {'c', 1, "opt_join", "expected fraction for band-optimization, joining", NULL, 0, 0, 0, 0, 0.0, 0.0, NULL},
#endif
  {'E', 1, "evalue", "E()-value threshold", "E()-value,E()-repeat threshold", 0, 0, 0, 0, 0.0, 0.0, NULL},
  {'f', 1, "gapopen", "gap-open penalty", NULL, 0, 0, 0, 0, 0.0, 0.0, NULL},
  {'g', 1, "gapext", "gap-extension penalty", NULL, 0, 0, 0, 0, 0.0, 0.0, NULL},
#ifdef SHOW_HELP
  {'h', 0, "help", "help - show options, arguments", NULL, 0, 0, 0, 0, 0.0, 0.0, NULL},
#ifdef FASTY
  {'j', 1, "frame_subs", "frame-shift, codon substitution penalty", NULL, 0, 0, 0, 0, 0.0, 0.0, NULL},
#else
#ifdef FASTX
  {'j', 1, "frame_shift", "frame-shift penalty", NULL, 0, 0, 0, 0, 0.0, 0.0, NULL},
#endif
#endif
#else
#ifndef LALIGN
  {'h', 1, "frame", "frameshift penalty", NULL, 0, 0, 0, 0, 0.0, 0.0, NULL},
  {'j', 1, "codon_subs", "codon substitution", NULL, 0, 0, 0, 0, 0.0, 0.0, NULL},
#endif
#endif
#if defined(LALIGN)
  {'J', 0, "show_ident", "show identity alignment", NULL, 0, 0, 0, 0, 0.0, 0.0, NULL},
  {'K', 1, "max_repeat", "maximum number of non-intersecting alignments", NULL, 0, 0, 0, 0, 0.0, 0.0, NULL},
#endif
  {'k', 1, "nshuffle", "number of shuffles", NULL, 0, 0, 0, 0, 0.0, 0.0, NULL},
  {'M', 1, "range", "filter on library sequence length", NULL, 0, 0, 0, 0, 0.0, 0.0, NULL},
  {'n', 0, "dna", "DNA/RNA query", NULL, 0, 0, 0, 0, 0.0, 0.0, NULL},
  {'p', 0, "prot", "protein query", NULL, 0, 0, 0, 0, 0.0, 0.0, NULL},
#if defined(CAN_PSSM)
  {'P', 1, "pssm", "PSSM file", NULL, 0, 0, 0, 0, 0.0, 0.0, NULL},
#endif
  {'r', 1, "dna_ratio", " +match/-mismatch for DNA/RNA", NULL, 0, 0, 0, 0, 0.0, 0.0, NULL},
  {'s', 1, "matrix", "scoring matrix", &s_opt_descr[0], 0, 0, 0, 0, 0.0, 0.0, NULL},
#if !defined(LALIGN) && !defined(FASTS) && !defined(FASTM)
  {'S', 0, "seg", "filter lowercase (seg) residues", NULL, 0, 0, 0, 0, 0.0, 0.0, NULL},
#endif
#if defined(TFAST) || defined(FASTX) || defined(FASTY)
  {'t', 1, "gencode", "translation genetic code", NULL, 0, 0, 0, 0, 0.0, 0.0, NULL},
#endif
  {'U', 0, "rna", "RNA query", NULL, 0, 0, 0, 0, 0.0, 0.0, NULL},
  {'X', 1, "ext_opts", "Extended options", NULL, 0, 0, 0, 0, 0.0, 0.0, NULL},
  {'z', 1, "stats", "Statistics estimation method", &z_opt_descr[0], 0, 0, 0, 0, 0.0, 0.0, NULL},
  {'\0', 0, "", "", NULL, 0, 0, 0, 0, 0.0, 0.0, NULL}
};

void f_init_opts(int pgm_id, struct mngmsg *m_msp, struct pstruct *ppst) {
#if defined(FASTA) || defined(FASTX) || defined(FASTY)
  set_opt_disp_defs('c', f_options, 4, 0, 0, ppst->param_u.fa.E_band_opt, max(ppst->param_u.fa.E_band_opt*5.0,1.0), NULL);
#endif
#ifndef LALIGN
  set_opt_disp_defs('E', f_options, 4, 0, 0, ppst->e_cut, ppst->e_cut_r, NULL);
#else
  set_opt_disp_defs('E', f_options, 3, 0, 0, ppst->e_cut, 0.0, NULL);
#endif
  set_opt_disp_defs('s', f_options, 5, 0, 0, 0.0, 0.0, ppst->pamfile);
  set_opt_disp_defs('f', f_options, 1, ppst->gdelval, 0, 0.0, 0.0, NULL);
  set_opt_disp_defs('g', f_options, 1, ppst->ggapval, 0, 0.0, 0.0, NULL);
#ifdef FASTY
  set_opt_disp_defs('j', f_options, 2, ppst->gshift, ppst->gsubs, 0.0, 0.0, NULL);
#endif
#ifdef FASTX
  set_opt_disp_defs('j', f_options, 1, ppst->gshift, 0, 0.0, 0.0, NULL);
#endif
#ifdef LALIGN
  set_opt_disp_defs('K', f_options, 1, ppst->max_repeat, 0, 0.0, 0.0, NULL);
#endif
  set_opt_disp_defs('k', f_options, 1, m_msp->shuff_max, 0, 0.0, 0.0, NULL);
  set_opt_disp_defs('r', f_options, 2, ppst->p_d_mat, ppst->p_d_mis, 0.0, 0.0, NULL);
  set_opt_disp_defs('t', f_options, 1, ppst->tr_type, 0, 0.0, 0.0, NULL);
  set_opt_disp_defs('z', f_options, 1, ppst->zsflag, 0, 0.0, 0.0, NULL);
#if defined(FASTA) || defined(SSEARCH)
  set_opt_disp_defs('W', f_options, 1, m_msp->aln.llcntx, 0, 0.0, 0.0, NULL);
#endif
}

char *iprompt1=" test sequence file name: ";
char *iprompt2=" database file name: ";

#ifdef PCOMPLIB
char *verstr="36.3.8e Dec, 2016 MPI";
#else
char *verstr="36.3.8e Dec, 2016";
#endif

static int mktup=3;
static int ktup_set = 0;
static int gap_set=0;
static int del_set=0;
static int mshuff_set = 0;
static int prot2dna = 0;
void
parse_ext_opts(char *opt_arg, int pgm_id, struct mngmsg *m_msp, struct pstruct *ppst);

extern int fa_max_workers;

extern void s_abort(char *, char *);
extern void init_ascii0(int *sascii, char *sq_map, int sq_map_n, struct pstruct *ppst);
extern void init_ascii(int ext_sq, int *sascii, int nsq, int dnaseq);
extern void validate_novel_aa(int *sascii, int p_nsq, int dnaseq);
extern int standard_pam(char *smstr, struct pstruct *ppst,
			int del_set, int gap_set);
extern int min_pam_bits(int n0_eff, double bit_thresh, struct pstruct *ppst,
			int del_set, int gap_set);
extern void mk_n_pam(int *arr,int siz, int mat, int mis);
extern int karlin(int , int, double *, double *, double *);
extern void init_karlin_a(struct pstruct *, double *, double **);
extern int do_karlin_a(int **, const struct pstruct *, double *,
		       double *, double *, double *, double *);

#if defined(TFAST) || defined(FASTX) || defined(FASTY)
extern void aainit(int tr_type, int debug);
#endif

char *iprompt0, *prog_func, *refstr;

/* Sets defaults assuming a protein sequence */
void h_init (struct pstruct *ppst, struct mngmsg *m_msp, char *info_pgm_abbr)
{
  struct pgm_def_str pgm_def;
  int i, pgm_id;

  ppst->pgm_id  = pgm_id =  get_pgm_id();
  pgm_def = pgm_def_arr[pgm_id];

  /* check that pgm_def_arr[] is valid */
  if (pgm_def.pgm_id != pgm_id) {
    fprintf(stderr,
	    "**pgm_def integrity failure: def.pgm_id %d != pgm_id %d**\n",
	    pgm_def.pgm_id, pgm_id);
    exit(1);
  }

  /* check that msg_def_arr[] is valid */
  if (msg_def_arr[pgm_id].pgm_id != pgm_id) {
    fprintf(stderr,
	    "**msg_def integrity failure: def.pgm_id %d != pgm_id %d**\n",
	    msg_def_arr[pgm_id].pgm_id, pgm_id);
    exit(1);
  }

  SAFE_STRNCPY(info_pgm_abbr,pgm_def.info_pgm_abbr,MAX_SSTR);
  iprompt0 = pgm_def.iprompt0;
  refstr = pgm_def.ref_str;
  prog_func = pgm_def.prog_func;

  /* used to be MAXTOT = MAXTST+MAXLIB, but now fixed at MAXLIB for
     pre-loaded libraries */
  m_msp->max_tot = MAXLIB;

  init_ascii0(aascii, NCBIstdaa, NCBIstdaa_n, ppst);
  pascii = aascii;

  /* set up DNA query sequence if required*/
  if (msg_def_arr[pgm_id].q_seqt == SEQT_DNA) {
    memcpy(qascii,nascii,sizeof(qascii));
    m_msp->qdnaseq = SEQT_DNA;
  }
  else { 	/* when SEQT_UNK, start with protein */
    memcpy(qascii,aascii,sizeof(qascii));
    m_msp->qdnaseq = msg_def_arr[pgm_id].q_seqt;
  }

#if defined(FASTF) || defined(FASTS) || defined(FASTM)
  qascii[','] = ESS;
  /* also initialize aascii, nascii for databases */
  qascii['*'] = NA;
  ppst->pam_ms = 1;
  ppst->do_rep=0;		/* disable multiple alignments */
  ppst->pseudocts = 200;
#else
  ppst->pam_ms = 0;
  ppst->do_rep=1;		/* enable multiple alignments */
#endif

  /* initialize a pam matrix */
  SAFE_STRNCPY(ppst->pamfile,pgm_def.smstr,MAX_FN);
  standard_pam(ppst->pamfile,ppst,del_set,gap_set);
  ppst->have_pam2 = 0;

  /* specify pre-alignment */
  ppst->can_pre_align = pgm_def.can_pre_align;

  ppst->p_d_mat = 5;
  ppst->p_d_mis = -4;

  /* this is always protein by default */
  ppst->nsq = NCBIstdaa_n;
  ppst->nsqx = NCBIstdaa_ext_n;
  /* we need to populate ppst->sq to nsqx for direct lc mapping */
  for (i=0; i<ppst->nsqx; i++) {
    ppst->sq[i] = NCBIstdaa_l[i];
    ppst->hsq[i] = h_NCBIstdaa[i];
  }
  for (i=0; i<ppst->nsqx; i++) {
    ppst->sqx[i]=NCBIstdaa_ext[i];
    ppst->hsqx[i]=h_NCBIstdaa_ext[i];
  }
  ppst->sq[ppst->nsq] = ppst->sqx[ppst->nsqx] = '\0';

  /* set up the c_nt[] mapping */

#if defined(FASTS) || defined(FASTF) || defined(FASTM)
  ppst->c_nt[ESS] = ESS;
#endif
  ppst->c_nt[0]=0;
  for (i=1; i<nnt; i++) {
    ppst->c_nt[i]=gc_nt[i];
    ppst->c_nt[i+nnt]=gc_nt[i]+nnt;
  }

#ifdef CAN_PSSM
  ppst->pam2p[0] = NULL;
  ppst->pam2p[1] = NULL;
#endif
}

/*
 * alloc_pam(): allocates memory for the 2D pam matrix as well
 * as for the integer array used to transmit the pam matrix
 */
void
alloc_pam (int d1, int d2, struct pstruct *ppst)
{
  int     i, *d2p;
  char err_str[128];

  if ((ppst->pam2[0] = (int **) malloc (d1 * sizeof (int *))) == NULL) {
     sprintf(err_str,"Cannot allocate 2D pam matrix: %d",d1);
     s_abort (err_str,"");
  }

  if ((ppst->pam2[1] = (int **) malloc (d1 * sizeof (int *))) == NULL) {
     sprintf(err_str,"Cannot allocate 2D pam matrix: %d",d1);
     s_abort (err_str,"");
  }

  if ((d2p = (int *) calloc (d1 * d2, sizeof (int))) == NULL) {
     sprintf(err_str,"Cannot allocate 2D pam matrix: %d",d1);
     s_abort (err_str,"");
   }

   for (i = 0; i < d1; i++, d2p += d2)
      ppst->pam2[0][i] = d2p;

   if ((d2p= (int *) malloc (d1 * d2 * sizeof (int))) == NULL) {
     sprintf(err_str,"Cannot allocate 2d pam matrix: %d",d2);
     s_abort (err_str,"");
   }

   for (i = 0;  i < d1; i++, d2p += d2)
      ppst->pam2[1][i] = d2p;

   ppst->have_pam2 = 1;
}

/*
 *  init_pam2(struct pstruct pst): Converts 1-D pam matrix to 2-D
 *  currently, this function is very protein centric
 */
void
init_pam2 (struct pstruct *ppst) {
  int     i, j, k, nsq, sa_t;
  int ix_j, ix_l, ix_i, p_i, p_j;

  nsq = ppst->nsq;

  ppst->pam2[0][0][0] = -BIGNUM;
  ppst->pam_h = -1; ppst->pam_l = 1;

  k = 0;

  if (ppst->dnaseq == 0) { /* not DNA */
    sa_t = aascii['*'];	/* this is the last character for which pam[] is available */
    pam_sq = apam_sq;
    pam_sq_n = apam_sq_n;
  }
  else {	/* have DNA, no '*' */
    sa_t = nascii['X'];
    pam_sq = npam_sq;
    pam_sq_n = npam_sq_n;
  }

  /* we use sa_t here because that is the last position in the 1-D
     matrix */
  for (i = 1; i < sa_t; i++) {
    p_i = pascii[pam_sq[i]];
    ppst->pam2[0][0][p_i] = ppst->pam2[0][p_i][0] = -BIGNUM;
    for (j = 1; j <= i; j++) {
      /* here is where the pam file is actually set */
      p_j = pascii[pam_sq[j]];
      ppst->pam2[0][p_j][p_i] = ppst->pam2[0][p_i][p_j] = pam[k++] - ppst->pamoff;
      if (ppst->pam_l > ppst->pam2[0][p_i][p_j]) ppst->pam_l = ppst->pam2[0][p_i][p_j];
      if (ppst->pam_h < ppst->pam2[0][p_i][p_j]) ppst->pam_h = ppst->pam2[0][p_i][p_j];
    }
  }

  /* need to do the same thing for characters > sa_t */
  for (i = sa_t+1; i < pam_sq_n; i++) {
    p_i = pascii[pam_sq[i]];
    ppst->pam2[0][0][p_i] = ppst->pam2[0][p_i][0] = -BIGNUM;
  }  

  if (ppst->dnaseq == 0) {
    init_altpam(ppst);
  }
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

/* extend the standard pam matrix for special residues
   (a) 'X' for protein and 'N' for DNA, 'G' and 'U' for RNA, 
   (b) lower-case characters for ext_sq_set

   must be called after init_pam2()
*/
void
init_pamx (struct pstruct *ppst) {
  int     i, j, k, nsq;
  int sa_x, sa_t, tmp;

  nsq = ppst->nsq;

  ppst->nt_align = (ppst->dnaseq== SEQT_DNA || ppst->dnaseq == SEQT_RNA);

  /* sa_x is the 'unknown' character, 'X' for proteins, 'N' for DNA/RNA */
  /* sa_t is the termination character -- not used for DNA */

  if (ppst->nt_align) {
    sa_x = pascii['N'];
    sa_t = sa_x;
  }
  else {
    sa_x = pascii['X'];
    sa_t = pascii['*'];
  }

  /* build an asymmetric matrix for RNA */
  if (ppst->dnaseq == SEQT_RNA && !ppst->pam_set) {
    tmp = ppst->pam2[0][nascii['G']][nascii['G']] - 3;
    ppst->pam2[0][nascii['A']][nascii['G']] = 
      ppst->pam2[0][nascii['C']][nascii['T']] = 
      ppst->pam2[0][nascii['C']][nascii['U']] = tmp;
  }

  if (ppst->pam_x_set) {
    for (i=1; i<nsq; i++) {
      if (sa_x < nsq) ppst->pam2[0][sa_x][i] = ppst->pam2[0][i][sa_x]=ppst->pam_xm;
      if (sa_t < nsq) ppst->pam2[0][sa_t][i] = ppst->pam2[0][i][sa_t]=ppst->pam_xm;
    }
    if (sa_x < nsq) ppst->pam2[0][sa_x][sa_x]=ppst->pam_xx;
    if (sa_t < nsq) ppst->pam2[0][sa_t][sa_t]=ppst->pam_xx;
  }
  else {
    ppst->pam_xx = ppst->pam2[0][sa_x][sa_x];
    ppst->pam_xm = ppst->pam2[0][1][sa_x];
  }

    /* fill in pam2[1] matrix */
    ppst->pam2[1][0][0] = -BIGNUM;
    /* fill in additional parts of the matrix */
    for (i = 0; i < nsq; i++) {
      /* -BIGNUM to all matches vs 0 */
      ppst->pam2[0][0][i+nsq] = ppst->pam2[0][i+nsq][0] = 
	ppst->pam2[1][0][i+nsq] = ppst->pam2[1][i+nsq][0] = 
	ppst->pam2[0][0][i] = ppst->pam2[0][i][0] = 
	ppst->pam2[1][0][i] = ppst->pam2[1][i][0] = -BIGNUM;

      for (j = 0; j < nsq; j++) {
	/* replicate pam2[0] to i+nsq, j+nsq, also initialize lowest of pam2[1]*/
	ppst->pam2[0][i+nsq][j] = ppst->pam2[0][i][j+nsq] =  ppst->pam2[0][i+nsq][j+nsq] = 
	  ppst->pam2[1][i][j] = ppst->pam2[0][i][j];

	/* set the high portion of pam2[1] to the corresponding value
	   of pam2[1][sa_x][j] */

	ppst->pam2[1][i+nsq][j] = ppst->pam2[1][i][j+nsq]=
	  ppst->pam2[1][i+nsq][j+nsq]=ppst->pam2[0][sa_x][j];
      }
    }

    /* set matches to the internal '-' to pam_xx */
    /* this needs to be adjusted for multiple internal '-----' */
    for (i=1; i< nsq; i++) {
      ppst->pam2[0][nsq][i] = ppst->pam2[0][i][nsq] = 
	ppst->pam2[0][nsq][i+nsq] = ppst->pam2[0][i+nsq][nsq] = 
	ppst->pam2[1][nsq][i] = ppst->pam2[1][i][nsq] = 
	ppst->pam2[1][nsq][i+nsq] = ppst->pam2[1][i+nsq][nsq] = ppst->pam_xm;
    }
    ppst->pam2[0][nsq][nsq] = ppst->pam2[1][nsq][nsq] = ppst->pam_xm;
}

/*  function specific initializations */
void
f_initenv (struct mngmsg *m_msp, struct pstruct *ppst, unsigned char **aa0) {
  struct msg_def_str m_msg_def;
  int pgm_id;

  pgm_id = ppst->pgm_id;
  m_msg_def = msg_def_arr[pgm_id];

  m_msp->last_calc_flg=0;

  SAFE_STRNCPY(m_msp->f_id0,m_msg_def.f_id0,sizeof(m_msp->f_id0));
  SAFE_STRNCPY(m_msp->f_id1,m_msg_def.f_id1,sizeof(m_msp->f_id1));
  SAFE_STRNCPY (m_msp->label, m_msg_def.label, sizeof(m_msp->label));
  SAFE_STRNCPY(m_msp->alabel, m_msg_def.alabel, sizeof(m_msp->alabel));

#if !defined(SSEARCH) && !defined(GGSEARCH) && !defined(GLSEARCH) && !defined(LALIGN)
  SAFE_STRNCPY (m_msp->alab[0],"initn",20);
  SAFE_STRNCPY (m_msp->alab[1],"init1",20);
  SAFE_STRNCPY (m_msp->alab[2],"opt",20);
#else
#if defined(SSEARCH) || defined(LALIGN)
  SAFE_STRNCPY (m_msp->alab[0],"s-w opt",20);
#else
  SAFE_STRNCPY (m_msp->alab[0],"n-w opt",20);
#endif
#endif

#if defined(GGSEARCH) || defined(GLSEARCH)
  m_msp->zsflag = ppst->zsflag = ppst->zsflag_f = 0;
#else
  m_msp->zsflag = ppst->zsflag = ppst->zsflag_f = 1;
  m_msp->zsflag2 = ppst->zsflag2 = 1;
#endif

  ppst->gdelval += pgm_def_arr[pgm_id].g_open_mod;
  ppst->ggapval += pgm_def_arr[pgm_id].g_ext_mod;
#if defined(FASTX) || defined(FASTY)
  ppst->gshift = pgm_def_arr[pgm_id].gshift;
  ppst->gsubs = pgm_def_arr[pgm_id].hshift;
#endif
  ppst->sw_flag = m_msg_def.sw_flag;
  ppst->e_cut = m_msp->e_cut=pgm_def_arr[pgm_id].e_cut;
#ifndef LALIGN
  ppst->e_cut_r = ppst->e_cut/10.0;	/* more significant */
#else
  ppst->e_cut_r = ppst->e_cut;		/* everything if local */
#endif

  ppst->score_ix = 0;
  ppst->histint = 2;
  m_msp->qframe = m_msg_def.qframe;
  m_msp->nframe = m_msg_def.nframe;
  m_msp->nrelv = m_msg_def.nrelv;
  m_msp->srelv = m_msg_def.srelv;
  m_msp->arelv = m_msg_def.arelv;
  m_msp->stages = m_msg_def.stages;
  m_msp->shuff_wid = 0;
#if defined(GGSEARCH)
  m_msp->shuff_max = 100;
#else
  m_msp->shuff_max = MAX_RSTATS;
#endif
  m_msp->shuff_max_save = m_msp->shuff_max;

  /* see param.h for the definition of all these */

  m_msp->qshuffle = 0;
  m_msp->nm0 = 1;
  m_msp->escore_flg = 0;

  /* pam information */
  ppst->pam_pssm = 0;
#if defined(FASTS) || defined(FASTF) || defined(FASTM)
   ppst->pam_xx = ppst->pam_xm = 0;
#else
  ppst->pam_xx = 1;  /* set >0 to use pam['X']['X'] value */
  ppst->pam_xm = -1;  /* set >0 to use pam['X']['A-Z'] value */
#endif
  ppst->pam_x_set = 0;
  ppst->pam_x_id_sim = 0;
  ppst->pam_set = ppst->pam_variable = 0;
  ppst->pam_pssm = 0;
  ppst->p_d_set = 0;
  ppst->pamoff = 0;
  ppst->ext_sq_set = 0;
  ppst->nsq_e = ppst->nsq;

  /* initial settings for protein */
  if (pgm_def_arr[ppst->pgm_id].ktup > 0) {
    mktup = 3;
    ppst->param_u.fa.bestscale = 300;
    ppst->param_u.fa.bestoff = 36;
    ppst->param_u.fa.bkfact = 6;
    ppst->param_u.fa.scfact = 3;
    ppst->param_u.fa.bktup = mktup;
    ppst->param_u.fa.ktup = 0;
    ppst->param_u.fa.bestmax = 50;
    ppst->param_u.fa.pamfact = 1;
    ppst->param_u.fa.altflag = 0;
    ppst->param_u.fa.optflag = 1;
    ppst->param_u.fa.iniflag = 0;
    ppst->param_u.fa.optcut = 0;
    ppst->param_u.fa.optcut_set = 0;
    ppst->param_u.fa.cgap = 0;
    ppst->param_u.fa.optwid = 16;
    ppst->param_u.fa.optwid_set = 0;
    ppst->param_u.fa.E_band_opt = pgm_def_arr[ppst->pgm_id].E_band_opt;
    ppst->param_u.fa.use_E_thresholds = 1; /* disable E-thresholds for now */
  }

  f_init_opts(pgm_id, m_msp, ppst);
}

/*  switches for fasta only */

static int shift_set=0;
static int subs_set=0;
static int sw_flag_set=0;
static int nframe_set=0;
static int E_thresh_set = 0;
static int E_cgap_set = 0;

void
f_getopt (char copt, char *optarg,
	  struct mngmsg *m_msg, struct pstruct *ppst)
{
  int pgm_id;
  double tmp_f, tmp_f1;
  double tmp_e_cut, tmp_e_rep;
  int dnaseq_save;
  char *bp;

  pgm_id = ppst->pgm_id;

  switch (copt) {
  case '3':
    nframe_set = 1;
    if (pgm_id == TFA_PID) {
      m_msg->nframe = 3; break;
    }
    else {
      m_msg->nframe = 1;	/* for TFASTXY */
      m_msg->qframe = 1;  /* for FASTA, FASTX */
    }
    break;
  case 'a': m_msg->aln.showall = 1; break;
  case 'A':
    if (ppst->sw_flag) ppst->sw_flag=0;
    else ppst->sw_flag= 1;
    sw_flag_set = 1;
    break;
  case 'b':
    if (optarg[0] == '$') {
      m_msg->mshow = -1;
      m_msg->e_cut = 10000000.0;
      break;
    }
    else if (optarg[0] == '=') {
      m_msg->e_cut = 10000000.0;
      m_msg->e_cut_set = 1;
      m_msg->mshow_min = 1;
      sscanf (optarg+1, "%d", &m_msg->mshow);
    }
    else if (optarg[0] == '>') {
      m_msg->mshow_min = 2;
      sscanf (optarg+1, "%d", &m_msg->mshow);
    }
    else {
      sscanf (optarg, "%d", &m_msg->mshow);
      m_msg->mshow_min = 0;
    }
    m_msg->mshow_set = 1;
    break;
  case 'c':
    tmp_f = tmp_f1 = 0.0;
    if (*optarg == 'O') {
      ppst->param_u.fa.use_E_thresholds = 0;
      optarg++;
    }
    if (*optarg != '\0' && pgm_def_arr[pgm_id].ktup > 0) {
      sscanf (optarg, "%lf %lf", &tmp_f, &tmp_f1);
      if (tmp_f > 1.0) {
	ppst->param_u.fa.optcut = (int)(tmp_f+0.1);
	ppst->param_u.fa.use_E_thresholds = 0;
	ppst->param_u.fa.optcut_set = 1;
      }
      else if (tmp_f <= 0.0) {
	ppst->param_u.fa.use_E_thresholds = 1;
      }
      else {	/* 0.0 < tmp_f <= 1.0 */
	ppst->param_u.fa.use_E_thresholds = 1;
	ppst->param_u.fa.E_band_opt = min(tmp_f,1.0);
	E_thresh_set = 1;
	if (tmp_f1 > 0.0) {
	  tmp_f1 = min(tmp_f1,1.0);	/* may want to do max(tmp_f1,tmp_f) */
	  ppst->param_u.fa.E_join = tmp_f1;
	  E_cgap_set=1;
	}
      }
    }
    break;
  case 'd': sscanf(optarg,"%d",&m_msg->ashow);
    if ((m_msg->mshow > 0) && (m_msg->ashow > m_msg->mshow)) m_msg->mshow=m_msg->ashow;
    m_msg->ashow_set = 1;
    break;
  case 'E':
    if (strchr(optarg,' ')) {	/* check for 1 or 2 values */
      sscanf(optarg,"%lf %lf",&tmp_e_cut, &tmp_e_rep);
      if (tmp_e_rep <= 0.0) {	/* two values, 2nd <= 0.0, no do_rep */
	ppst->do_rep = 0;
	ppst->e_cut_r = 1E-100;
	tmp_e_rep = -2.0;
      }
      else {ppst->do_rep = 1;}
    }
    else {	/* one value, do_rep; tmp_e_rep=10.0 */
      sscanf(optarg,"%lf",&tmp_e_cut);
#ifndef LALIGN
      tmp_e_rep = 10.0;
#else
      tmp_e_rep = 1.0;
#endif
      ppst->do_rep = 1;
    }
    if (!m_msg->e_cut_set && tmp_e_cut > 0.0 ) {
      ppst->e_cut = m_msg->e_cut = tmp_e_cut;
    }
    m_msg->e_cut_set = 1;

    if (tmp_e_rep > 0.0) {
      if (tmp_e_rep >= 1.0) { ppst->e_cut_r = ppst->e_cut/tmp_e_rep;}
      else { ppst->e_cut_r = tmp_e_rep;}
    }
    break;
  case 'f':
    sscanf (optarg, "%d", &ppst->gdelval);
    if (ppst->gdelval > 0) ppst->gdelval = -ppst->gdelval;
    del_set = 1;
    break;
  case 'g':
    sscanf (optarg, "%d", &ppst->ggapval);
    if (ppst->ggapval > 0) ppst->ggapval = -ppst->ggapval;
    gap_set = 1;
    break;
#ifndef SHOW_HELP
  case 'h':
    sscanf (optarg, "%d", &ppst->gshift);
    if (ppst->gshift > 0) ppst->gshift = -ppst->gshift;
    shift_set = 1;
    break;
  case 'j':
    sscanf (optarg, "%d", &ppst->gsubs);
    if (ppst->gsubs > 0) ppst->gsubs = -ppst->gsubs;
    subs_set = 1;
    break;
#else
  case 'h':
    show_help(m_msg->pgm_name, pgm_id);
    break;
  case 'j':
#ifdef FASTY
    if (strchr(optarg,' ')) {
      sscanf (optarg, "%d %d", &ppst->gshift, &ppst->gsubs);
      subs_set = 1;
      if (ppst->gsubs > 0) ppst->gsubs = -ppst->gsubs;
    }
    else if (strchr(optarg,',')) {
      sscanf (optarg, "%d,%d", &ppst->gshift, &ppst->gsubs);
      subs_set = 1;
      if (ppst->gsubs > 0) ppst->gsubs = -ppst->gsubs;
    }
    else {
      sscanf (optarg, "%d", &ppst->gshift);
    }
#else
#ifdef FASTX
    sscanf (optarg, "%d", &ppst->gshift);
#endif
#endif
    if (ppst->gshift > 0) ppst->gshift = -ppst->gshift;
    shift_set = 1;
    break;
#endif
  case 'J':
#ifdef LALIGN
    ppst->show_ident=1;
#else
    ppst->show_ident=0;
#endif
    break;


#ifdef LALIGN
  case 'K':
    sscanf(optarg,"%d", &ppst->max_repeat);
    break;
#endif
  case 'k':
    sscanf (optarg, "%d", &m_msg->shuff_max);
    m_msg->shuff_max_save = m_msg->shuff_max;
    mshuff_set = 1;
    break;
  case 'M':
    sscanf(optarg,"%d-%d",&m_msg->n1_low,&m_msg->n1_high);
    if (m_msg->n1_low < 0) {
      m_msg->n1_high = -m_msg->n1_low;
      m_msg->n1_low = 0;
    }
    if (m_msg->n1_high == 0) m_msg->n1_high = BIGNUM;
    if (m_msg->n1_low > m_msg->n1_high) {
      fprintf(stderr," low cutoff %d greater than high %d\n",
	      m_msg->n1_low, m_msg->n1_high);
      m_msg->n1_low = 0;
      m_msg->n1_high = BIGNUM;
    }
    ppst->n1_low = m_msg->n1_low;
    ppst->n1_high = m_msg->n1_high;
    break;
  case 'n':
    m_msg->qdnaseq = SEQT_DNA;
    re_ascii(qascii,nascii,strlen((char *)m_msg->ann_arr+1));
    SAFE_STRNCPY(m_msg->sqnam,"nt",4);
    prot2dna = 1;
    break;
  case 'o':
  case 'p':
    m_msg->qdnaseq = SEQT_PROT;
    ppst->dnaseq = SEQT_PROT;
    SAFE_STRNCPY(m_msg->sqnam,"aa",4);
    break;
  case 'P':
    SAFE_STRNCPY(ppst->pgpfile,optarg,MAX_FN);
    if ((bp=strchr(ppst->pgpfile,' '))!=NULL) {
      *bp='\0';
      ppst->pgpfile_type = atoi(bp+1);
    }
    else ppst->pgpfile_type = 0;
    ppst->pam_pssm = 1;
    break;
  case 'r':
    sscanf(optarg,"%d/%d",&ppst->p_d_mat,&ppst->p_d_mis);
    ppst->pam_set = 0;
    ppst->p_d_set = 1;

    SAFE_STRNCPY(ppst->pam_name, "DNA", 4);
    if (ppst->dnaseq != SEQT_RNA) ppst->dnaseq = SEQT_DNA;
    if (ppst->p_d_mat > 0 && ppst->p_d_mis < 0) {
      ppst->p_d_set = 1;
      SAFE_STRNCPY(ppst->pamfile,optarg,40);
    }
    break;
    /* modified Sept, 2011, to recognize that a scoring matrix
       specifies a sequence alphabet */
  case 's':
    if (*optarg == '?') {
      ppst->pam_variable = 1;
      optarg++;
    }
    if (*optarg == '\0') break;
    SAFE_STRNCPY (ppst->pamfile, optarg, MAX_FN);
    dnaseq_save = ppst->dnaseq;
    /* check for default abbreviation */
    if (!standard_pam(ppst->pamfile,ppst,del_set, gap_set)) {
      /* check/load matrix file */
      if (!initpam (ppst->pamfile, ppst)) {
	/* matrix file failed, use default matrix */
	SAFE_STRNCPY(ppst->pamfile,pgm_def_arr[pgm_id].smstr,MAX_FN);
      }
    }
    ppst->pam_set=1;
    /* check for changing alphabet here */
    if (ppst->dnaseq != dnaseq_save && ppst->dnaseq >= SEQT_DNA) {
      m_msg->qdnaseq = SEQT_DNA;
      re_ascii(qascii,nascii,strlen((char *)m_msg->ann_arr+1));
      SAFE_STRNCPY(m_msg->sqnam,"nt",4);
      prot2dna = 1;
    }
    break;
  case 'S':	/* turn on extended alphabet for seg */
    ppst->ext_sq_set = 1;
    ppst->nsq_e = ppst->nsqx;
    break;
  case 't':
    if (tolower(optarg[0])=='t') {
      m_msg->ldb_info.term_code = aascii['*'];
      optarg++;
    }
    if (*optarg) {sscanf (optarg, "%d", &ppst->tr_type);}
    break;
  case 'U':
    m_msg->qdnaseq = SEQT_RNA;
    memcpy(qascii,nascii,sizeof(qascii));
    SAFE_STRNCPY(m_msg->sqnam,"nt",4);
    nt[nascii['T']]='U';
    prot2dna=1;
    break;
  case 'W':
    sscanf (optarg,"%d",&m_msg->aln.llcntx);
    m_msg->aln.llcntx_set = 1;
    break;
  case 'X':
    parse_ext_opts(optarg, pgm_id, m_msg, ppst);
    break;
  case 'z':
    if (strchr(optarg,' ')!=NULL) {
      sscanf(optarg,"%d %d",&ppst->zsflag,&ppst->zsflag2);
      if (ppst->zsflag2 < 1 || ppst->zsflag2 > 6) ppst->zsflag2 = 2;
    }
    else if (strchr(optarg,',')!=NULL) {
      sscanf(optarg,"%d,%d",&ppst->zsflag,&ppst->zsflag2);
      if (ppst->zsflag2 < 1 || ppst->zsflag2 > 6) ppst->zsflag2 = 2;
    }
    else {
      sscanf(optarg,"%d",&ppst->zsflag);
      ppst->zsflag2 = (ppst->zsflag % 10);
    }
    break;
  }
}

static char my_opts[] = "1BIM:ox:y:N:";

void
parse_ext_opts(char *opt_arg, int pgm_id, struct mngmsg *m_msp, struct pstruct *ppst) {
  long l_arg;
  char c_arg, c_opt, *the_arg, *bp;

  c_opt = *opt_arg;
  if ((bp=strchr(my_opts, c_opt))==NULL) {
    return;
  }

  if (*(bp+1) == ':') the_arg = opt_arg+1;

  switch (c_opt) {
  case '1': 
    if (pgm_def_arr[pgm_id].ktup > 0) {
      ppst->param_u.fa.iniflag=1;
    }
    break;
  case 'B': m_msp->z_bits = 0; break;
  case 'I': 
    m_msp->tot_ident = 1;
    /*
    l_arg = 0;
    sscanf(the_arg,"%ld",&l_arg);
    if (l_arg > 0) m_msp->tot_ident = l_arg;
    */
    break;
  case 'M':
    c_arg = '\0';
    sscanf(the_arg,"%ld%c",&l_arg,&c_arg);
    if (l_arg < 0) m_msp->max_memK = BIGNUM;
    else {
      l_arg *= 1024;
      if (c_arg == 'G') l_arg *= 1024;
      m_msp->max_memK = l_arg;
    }
    break;
  case 'N':
  case 'X':
    ppst->pam_x_id_sim = 0;
    if (*the_arg == 'S' || *the_arg == '+') {
      ppst->pam_x_id_sim = 1;
    }
    else if (*the_arg == 'D' || *the_arg == '-') {
      ppst->pam_x_id_sim = -1;
    }
    break;
  case 'o':
    if (pgm_def_arr[pgm_id].ktup > 0) {
      ppst->param_u.fa.optflag = 0;
      msg_def_arr[pgm_id].nrelv = m_msp->nrelv = 2;
    }
    break;
  case 'x':
    if (strchr(the_arg,' ')!=NULL) {
      sscanf (the_arg,"%d %d",&ppst->pam_xx, &ppst->pam_xm);
    }
    else if (strchr(the_arg,',')!=NULL) {
      sscanf (the_arg,"%d,%d",&ppst->pam_xx, &ppst->pam_xm);
    }
    else {
      sscanf (the_arg,"%d",&ppst->pam_xx);
      ppst->pam_xm = ppst->pam_xx;
    }
    ppst->pam_x_set=1;
    break;
  case 'y':
    if (pgm_def_arr[pgm_id].ktup > 0) {
      sscanf (the_arg, "%d", &ppst->param_u.fa.optwid);
      ppst->param_u.fa.optwid_set = 1;
    }
    break;
  }
}

void
f_lastenv (struct mngmsg *m_msg, struct pstruct *ppst)
{
  char save_str[MAX_SSTR];

#if !defined(FASTM) && !defined(FASTS) && !defined(FASTF)
  SAFE_STRNCPY(save_str,"*",sizeof(save_str));
#else
  SAFE_STRNCPY(save_str,",",sizeof(save_str));
#endif

  if (m_msg->qdnaseq == SEQT_UNK) {
    build_xascii(qascii,save_str);
    if (m_msg->ann_flg) add_ascii_ann(qascii,m_msg->ann_arr);
  }  
/* this check allows lc DNA sequence queries with FASTX */
  else {
#if !defined(FASTS) && !defined(FASTM) && !defined(FASTF) && !defined(FASTX) && !defined(FASTY)
    init_ascii(ppst->ext_sq_set,qascii, ppst->nsq, m_msg->qdnaseq);
#endif
    validate_novel_aa(qascii, ppst->nsq, m_msg->qdnaseq);
  }
}

void
f_getarg (int argc, char **argv, int optind,
	  struct mngmsg *m_msg, struct pstruct *ppst)
{

  if (pgm_def_arr[ppst->pgm_id].ktup > 0) {
    if (argc - optind >= 4) {
      sscanf (argv[optind + 3], "%d", &ppst->param_u.fa.ktup);
      ktup_set = 1;
    }
    else {
      ppst->param_u.fa.ktup = -pgm_def_arr[ppst->pgm_id].ktup;
    }
  }
  
  if (ppst->pgm_id == RSS_PID && argc - optind > 3) {
    sscanf (argv[optind + 3], "%d", &m_msg->shuff_max);
  }

  if (ppst->pgm_id == RFX_PID && argc - optind > 4) {
    sscanf (argv[optind + 4], "%d", &m_msg->shuff_max);
  }
  m_msg->shuff_max_save = m_msg->shuff_max;
}

/* fills in the query ascii mapping from the parameter
   ascii mapping.
*/

void
re_ascii(int *qascii, int *pascii, int max_ann_arr) {
  int i;

  for (i=0; i < 128; i++) {
    if (qascii[i] > NANN+max_ann_arr || qascii[i] < ESS) {
      qascii[i] = pascii[i];
    }
  }
}


/* recode has become function specific to accommodate FASTS/M */
/* modified 28-Dec-2004 to ensure that all mapped characters
   are valid */
int
recode(unsigned char *seq, int n, int *qascii, int nsqx) {
  int i,j;
  char save_c;

#if defined(FASTS) || defined(FASTM)
  qascii[',']=ESS;
#endif

  for (i=0; i < n; i++) {
    save_c = seq[i];
    if (seq[i] > '@' || seq[i]=='*') seq[i] = qascii[seq[i]];
    if (seq[i] > nsqx && seq[i]!=ESS) {
      fprintf(stderr, "*** Warning - unrecognized residue at %d:%c - %2d\n",
	      i,save_c,save_c);
      seq[i] = qascii['X'];
    }
  }
  seq[i]=EOSEQ;
  return i;
}

/* here we have the query sequence, all the command line options,
   but we need to set various parameter options based on the type
   of the query sequence (m_msg->qdnaseq = 0:protein/1:DNA) and
   the function (FASTA/FASTX/TFASTA)

   29-Jun-2008 add code to ensure that weird ('O', 'U') amino-acids
   are read properly.

   15-Nov-2010 -- modify scoring matrix for very short query sequences
   (e.g. short read metagenomics)
*/

/* this resetp is for conventional a FASTA/TFASTXYZ search */
void
resetp (struct mngmsg *m_msg, struct pstruct *ppst) {
  int i, pgm_id;
  int n0_eff;

  pgm_id = ppst->pgm_id;

  /* check for alphabet conflict */

  ppst->shuffle_dna3 = 0;
#if defined(TFAST)
  if (m_msg->qdnaseq == SEQT_DNA || m_msg->qdnaseq == SEQT_RNA) {
    fprintf(stderr," %s compares a protein to a translated\n\
DNA sequence library.  Do not use a DNA query/scoring matrix.\n",prog_func);
    exit(1);
  }
  ppst->shuffle_dna3 = 1;
#else
#if (defined(FASTX) || defined(FASTY))
  if (!(m_msg->qdnaseq == SEQT_DNA || m_msg->qdnaseq == SEQT_RNA)) {
    fprintf(stderr," FASTX/Y compares a DNA sequence to a protein database\n");
    fprintf(stderr," Use a DNA query\n");
    exit(1);
  }
#endif
#endif

  /* **************************************************************** */
  /* adjust alphabets for prot:prot or DNA:DNA alignments */

  /* this code changes parameters for programs (FA_PID, SS_PID, FS_PID,
     RSS_PID) that can examine either protein (initial state) or DNA 
     Modified May, 2006 to reset e_cut for DNA comparisons.
  */
  /* **************************************************************** */

  if (msg_def_arr[pgm_id].q_seqt == SEQT_UNK) {
    if (m_msg->qdnaseq == SEQT_DNA || m_msg->qdnaseq == SEQT_RNA) {
      msg_def_arr[pgm_id].q_seqt = m_msg->qdnaseq;
      msg_def_arr[pgm_id].p_seqt = SEQT_DNA;
      msg_def_arr[pgm_id].l_seqt = SEQT_DNA;
      if (m_msg->qdnaseq == SEQT_DNA) msg_def_arr[pgm_id].qframe = 2;
      if (!m_msg->e_cut_set) {
	pgm_def_arr[pgm_id].e_cut /= 5.0;
	ppst->e_cut_r = 0.001;
      }
    }
    else {
      msg_def_arr[pgm_id].q_seqt = SEQT_PROT;
    }
  }

  /* set the comparison type (PROT/DNA) in ppst */
  ppst->dnaseq = msg_def_arr[pgm_id].p_seqt;

  if (!sw_flag_set) ppst->sw_flag = msg_def_arr[pgm_id].sw_flag;
  if (!m_msg->e_cut_set) {
    ppst->e_cut = m_msg->e_cut=pgm_def_arr[pgm_id].e_cut;
#ifdef LALIGN
    ppst->e_cut_r = ppst->e_cut;
#endif
  }

  if (ppst->dnaseq == SEQT_DNA && m_msg->qdnaseq==SEQT_RNA) {
    ppst->dnaseq = SEQT_RNA;
    ppst->nt_align = 1;
  }
  if (ppst->dnaseq==SEQT_DNA) pascii = &nascii[0];
  else if (ppst->dnaseq==SEQT_RNA) {
    pascii = &nascii[0];
    ppst->sq[nascii['T']] = 'U';
  }
  else pascii = &aascii[0];
  m_msg->ldb_info.ldnaseq = msg_def_arr[pgm_id].l_seqt;

  if (m_msg->ldb_info.ldnaseq & SEQT_DNA) {
    memcpy(lascii,nascii,sizeof(lascii));
#ifndef TFAST
#ifdef DNALIB_LC
    init_ascii(ppst->ext_sq_set,lascii, ppst->nsq, m_msg->ldb_info.ldnaseq);
#endif
#else
  /* no init_ascii() because we translate lower case library sequences */
#endif
    validate_novel_aa(lascii, ppst->nsq, m_msg->ldb_info.ldnaseq);
  }
  else {
    memcpy(lascii,aascii,sizeof(lascii));	/* initialize lib mapping */
    if (m_msg->ann_flg && strchr((char *)m_msg->ann_arr,'*')) {lascii['*'] = NA;}

#if defined(FASTF) || defined(FASTS) || defined(FASTM)
    lascii['*'] = NA;
#endif
    init_ascii(ppst->ext_sq_set,lascii, ppst->nsq, m_msg->ldb_info.ldnaseq);
    validate_novel_aa(lascii, ppst->nsq, m_msg->ldb_info.ldnaseq);
  }

  /* have lascii - initialize l_ann_ascii[] if necessary */
  if (m_msg->ann_flg) {
    memcpy(l_ann_ascii,lascii,sizeof(l_ann_ascii));
    /* make certain that '*' is treated correctly */
    if (strchr((char *)m_msg->ann_arr,'*')) {l_ann_ascii['*'] = NA;}
    add_ascii_ann(l_ann_ascii, m_msg->ann_arr);
  }

  /* **************************************************************** */
  /* adjust qframe/nframe if DNA/translated DNA search  */
  /* **************************************************************** */

  if (!nframe_set) {
    m_msg->qframe = msg_def_arr[pgm_id].qframe;
    m_msg->nframe = msg_def_arr[pgm_id].nframe;
  }

  /* the possibilities:
  	     -i  -3	qframe	revcomp
   FA_D/FX   -    -        2       0	
   FA_D/FX   +    -        2       1	
   FA_D/FX   -    +        1       0	
   FA_D/FX   +    +        2       1	
  */

  if (m_msg->qdnaseq == SEQT_DNA) {
    m_msg->nframe = 1;
    if (m_msg->qframe == 1 && m_msg->revcomp==1) {
      m_msg->qframe = m_msg->revcomp+1;
    }
  }
  else if (m_msg->qdnaseq == SEQT_RNA) {
    m_msg->qframe = m_msg->revcomp+1;
    m_msg->nframe = 1;
  }

  /* **************************************************************** */
  /* adjust FASTA heuristics for  DNA/translated DNA search  */
  /* **************************************************************** */

  if (ppst->dnaseq == SEQT_DNA || ppst->dnaseq == SEQT_RNA) {
    ppst->histint = 4;

    if (!del_set) {
#ifdef OLD_FASTA_GAP
      ppst->gdelval = -16;	/* def. del penalty */
#else
      ppst->gdelval = -12;	/* def. open penalty */
#endif
    }
    if (!gap_set) ppst->ggapval = -4;	/* def. gap penalty */

    ppst->nsq = nnt;
    ppst->nsqx = nntx;
    ppst->sq[ppst->nsqx+1] = ppst->sqx[ppst->nsqx+1] = '\0';

    /* reset parameters for DNA */
    if (pgm_def_arr[pgm_id].ktup > 0) {
      /* these parameters are used to scale optcut, and are being replaced
	 by statistically based parameters */
      /* largest ktup */
      pgm_def_arr[pgm_id].ktup = mktup = 6;
      if (!ppst->param_u.fa.optwid_set) ppst->param_u.fa.optwid = 16;
      ppst->param_u.fa.bestscale = 80;
      ppst->param_u.fa.bkfact = 5;
      ppst->param_u.fa.scfact = 1;
      ppst->param_u.fa.bktup = mktup;
      ppst->param_u.fa.bestmax = 80;
      ppst->param_u.fa.bestoff = 45;
      if (!E_thresh_set) ppst->param_u.fa.E_band_opt = 0.05;

      if (!sw_flag_set) {
	ppst->sw_flag = 0;
	SAFE_STRNCPY(m_msg->f_id1,"bs",sizeof(m_msg->f_id1));
	SAFE_STRNCPY(m_msg->alabel, align_label[1], sizeof(m_msg->alabel));
      }

      /* largest ktup */
      mktup = 6;

      if (ppst->param_u.fa.pamfact >= 0) ppst->param_u.fa.pamfact = 0;
      if (ppst->param_u.fa.ktup < 0)
	ppst->param_u.fa.ktup = -ppst->param_u.fa.bktup;
    }

    for (i=0; i<=ppst->nsqx; i++) {
      ppst->hsq[i] = hnt[i];
      ppst->sq[i] = nt[i];
      ppst->hsqx[i] = hntx[i];
      ppst->sqx[i] = ntx[i];
    }

    /* **************************************************************** */
    /* adjust scoring matrix for DNA:DNA search  */
    /* **************************************************************** */

    if (!ppst->pam_set) {
      if (ppst->p_d_set)
	mk_n_pam(npam,nnt,ppst->p_d_mat,ppst->p_d_mis);
#if !defined(FASTS) && !defined(FASTM)
      else if (ppst->pamfile[0]=='\0' || strncmp(ppst->pamfile,"BL50",4)==0) {
	SAFE_STRNCPY (ppst->pamfile, "+5/-4", sizeof(ppst->pamfile));
	SAFE_STRNCPY(ppst->pamfile_save, ppst->pamfile, sizeof(ppst->pamfile_save));
	SAFE_STRNCPY (ppst->pam_name, "+5/-4", sizeof(ppst->pamfile));
      }
#else
      else if (strncmp(ppst->pamfile,"MD20",4)==0) {
	SAFE_STRNCPY (ppst->pamfile, "+2/-2", sizeof(ppst->pamfile));
	SAFE_STRNCPY (ppst->pam_name, "+2/-2", sizeof(ppst->pam_name));
	SAFE_STRNCPY(ppst->pamfile_save, ppst->pamfile, sizeof(ppst->pamfile_save));
	ppst->p_d_mat = +2;
	ppst->p_d_mis = -2;
      	mk_n_pam(npam,nnt,ppst->p_d_mat,ppst->p_d_mis);
      }
#endif
      pam = npam;
    }

    SAFE_STRNCPY (m_msg->sqnam, "nt",sizeof(m_msg->sqnam));
    SAFE_STRNCPY (m_msg->sqtype, "DNA",sizeof(m_msg->sqtype));
  }	/* end DNA reset */

  else {  /* other parameters for protein comparison */
    if (pgm_def_arr[pgm_id].ktup > 0) {
      if (!ppst->param_u.fa.optwid_set) {
	if (ppst->param_u.fa.ktup==1) ppst->param_u.fa.optwid = 32;
	else ppst->param_u.fa.optwid = 16;
      }
    }
    if (!shift_set) {ppst->gshift = pgm_def_arr[pgm_id].gshift;}
    if (!subs_set) {ppst->gsubs = pgm_def_arr[pgm_id].hshift;}
  }

  SAFE_STRNCPY(ppst->pamfile_save, ppst->pamfile, 120);
}

/* query_parm()	this function asks for any additional parameters
	that have not been provided.  Could be null. */
void
query_parm (struct mngmsg *m_msp, struct pstruct *ppst)
{
   char    qline[40];

   if (pgm_def_arr[ppst->pgm_id].ktup > 0) {
     if (ppst->param_u.fa.ktup < 0)
       ppst->param_u.fa.ktup = -ppst->param_u.fa.ktup;

     if (ppst->param_u.fa.ktup == 0) {
       printf (" ktup? (1 to %d) [%d] ", mktup, pgm_def_arr[ppst->pgm_id].ktup);
       if (fgets (qline, sizeof(qline), stdin) == NULL) exit (0);
       else sscanf(qline,"%d",&ppst->param_u.fa.ktup);
     }
     if (ppst->param_u.fa.ktup == 0)
       ppst->param_u.fa.ktup = pgm_def_arr[ppst->pgm_id].ktup;
     else ktup_set = 1;
   }

#if defined(PRSS)
   if (m_msp->shuff_max < 10) m_msp->shuff_max = MAX_RSTATS;

   if (!mshuff_set) {
     printf(" number of shuffles [%d]? ",m_msp->shuff_max);
     fflush(stdout);
     if (fgets (qline, sizeof(qline), stdin) == NULL) exit (0);
     else sscanf(qline,"%d",&m_msp->shuff_max);
   }

   if (ppst->zs_win == 0) {
     printf (" local (window) (w) or uniform (u) shuffle [u]? ");
     if (fgets (qline, sizeof(qline), stdin) == NULL) exit (0);
     else if (qline[0]=='w' || qline[0]=='W') {
       m_msp->shuff_wid = 20;
       printf(" local shuffle window size [%d]? ",m_msp->shuff_wid);
       if (fgets (qline, sizeof(qline), stdin) == NULL) exit (0);
       else sscanf(qline,"%d",&m_msp->shuff_wid);
     }
   }
#endif
}

/* last_init() cannot look at aa0, n0, because it is only run once,
   it is not run before each new aa0 search */
void
last_init (struct mngmsg *m_msg, struct pstruct *ppst)
{
  int ix_l, ix_i, i, pgm_id;
  double *kar_p;
  double aa0_f[MAXSQ];

  m_msg->zsflag = ppst->zsflag;
  m_msg->zsflag2 = ppst->zsflag2;

  if (ppst->zsflag < 0) {
    ppst->do_rep = 0;
  }

  pgm_id = ppst->pgm_id;

#ifdef LALIGN
  m_msg->do_showbest = 1;
  m_msg->quiet = 1;
#endif

#if defined(FASTF) || defined(FASTS) || defined(FASTM)
  m_msg->nohist = 1;
  m_msg->shuff_max = 2000;
  ppst->shuff_node = m_msg->shuff_max/fa_max_workers;
#else
  m_msg->shuff_max = m_msg->shuff_max_save;
#endif

  if (m_msg->aln.llen < 1) {
    m_msg->aln.llen = 60;
  }

  if (m_msg->ldb_info.ldnaseq== SEQT_PROT) {
    m_msg->max_tot = MAXLIB_P;
  }

#if defined(FASTX) || defined(FASTY) || defined(TFAST)
  /* set up translation tables: faatran.c */
  aainit(ppst->tr_type,ppst->debug_lib);
#endif

/* a sanity check */
#if !defined(TFAST)
   if (m_msg->revcomp && m_msg->qdnaseq!=SEQT_DNA && m_msg->qdnaseq!=SEQT_RNA) {
     fprintf(stderr," cannot reverse complement protein\n");
     m_msg->revcomp = 0;
   }
#endif

   if (pgm_def_arr[pgm_id].ktup > 0) {

     if (ppst->param_u.fa.ktup < 0)
       ppst->param_u.fa.ktup = -ppst->param_u.fa.ktup;

     if (ppst->param_u.fa.ktup < 1 || ppst->param_u.fa.ktup > mktup) {
       fprintf(stderr," warning ktup = %d out of range [1..%d], reset to %d\n",
	       ppst->param_u.fa.ktup, mktup, ppst->param_u.fa.bktup);
       ppst->param_u.fa.ktup = ppst->param_u.fa.bktup;
     }

     if (ppst->sw_flag) {
       SAFE_STRNCPY(m_msg->f_id1,"sw",sizeof(m_msg->f_id1));
       SAFE_STRNCPY(m_msg->alabel, align_label[0], sizeof(m_msg->alabel));
     }
     else {
       SAFE_STRNCPY(m_msg->f_id1,"bs",sizeof(m_msg->f_id1));
       SAFE_STRNCPY(m_msg->alabel, align_label[1], sizeof(m_msg->alabel));
     }
   }

   if (pgm_id == TFA_PID) {
     m_msg->revcomp *= 3;
     if (m_msg->nframe == 3) m_msg->nframe += m_msg->revcomp;
   }
   else if (pgm_id == TFX_PID || pgm_id == TFY_PID) {
     if (m_msg->nframe == 1) m_msg->nframe += m_msg->revcomp;
   }

#if !defined(TFAST)
  /* for fasta/fastx searches, itt iterates the the query strand */
  m_msg->nitt1 = m_msg->qframe-1;
#else
  /* for tfasta/tfastxy searches, itt iterates the library frames */
  m_msg->nitt1 = m_msg->nframe-1;
#endif

  if (pgm_def_arr[pgm_id].ktup > 0) {	/* its FASTA, not SSEARCH */
    if (ppst->param_u.fa.ktup>=2 && !ppst->param_u.fa.optwid_set) {
      ppst->param_u.fa.optwid=16;
      switch (pgm_id) {
      case FA_PID:
      case FX_PID:
      case FY_PID:
	m_msg->thr_fact = 8;
	m_msg->thr_fact = 8;
	break;
      case TFA_PID:
      case TFX_PID:
      case TFY_PID:
	m_msg->thr_fact = 4;
	break;
      default:
	m_msg->thr_fact = 4;
      }
    }
    else { m_msg->thr_fact = 4;}
  }
  else {
#if !defined(SW_ALTIVEC) && !defined(SW_SSE2)
    m_msg->thr_fact = 1;	/* unvectorized SSEARCH  */
#else
    m_msg->thr_fact = 8;	/* vectorized SSEARCH */
#endif
  }

#ifdef PCOMPLIB
  m_msg->thr_fact = 1;	/* use much larger buffers */
#endif

#if defined(PRSS)
   if (m_msg->shuff_max < 10) m_msg->shuff_max = MAX_RSTATS;
   if (ppst->zsflag < 10) ppst->zsflag += 10;
   if (ppst->zs_win > 0) {
     m_msg->shuff_wid = ppst->zs_win;
   }
#endif

   if (pgm_def_arr[ppst->pgm_id].ktup > 0) {
     if (ppst->param_u.fa.iniflag) {
       ppst->score_ix = 1;
       SAFE_STRNCPY (m_msg->label, "initn init1", sizeof(m_msg->label));
     }
     else if (ppst->param_u.fa.optflag) {
       ppst->score_ix = 2;
       m_msg->stages = 1;
     }
   }

   if (!ppst->have_pam2) {
     alloc_pam (MAXSQ, MAXSQ, ppst);
     init_pam2(ppst);
   }
   init_pamx(ppst);

   if (ppst->pam_ms) {
     if (m_msg->qdnaseq == SEQT_PROT) {
       /* code to make 'L'/'I' identical scores */
       ix_l = pascii['L'];
       ix_i = pascii['I'];
       ppst->pam2[0][ix_l][ix_i] = ppst->pam2[0][ix_i][ix_l] =
	 ppst->pam2[0][ix_l][ix_l] = ppst->pam2[0][ix_i][ix_i] =
	 max(ppst->pam2[0][ix_l][ix_l],ppst->pam2[0][ix_i][ix_i]);
       for (i=1; i<=ppst->nsq; i++) {
	 ppst->pam2[0][i][ix_i] = ppst->pam2[0][i][ix_l] =
	   max(ppst->pam2[0][i][ix_l],ppst->pam2[0][i][ix_i]);
	 ppst->pam2[0][ix_i][i] = ppst->pam2[0][ix_l][i] =
	   max(ppst->pam2[0][ix_i][i],ppst->pam2[0][ix_l][i]);
       }

       /* code to make 'Q'/'K' identical scores */
       if (!shift_set) {
	 ix_l = pascii['Q'];
	 ix_i = pascii['K'];
	 ppst->pam2[0][ix_l][ix_i] = ppst->pam2[0][ix_i][ix_l] =
	   ppst->pam2[0][ix_l][ix_l] = ppst->pam2[0][ix_i][ix_i] =
	   (ppst->pam2[0][ix_l][ix_l]+ppst->pam2[0][ix_i][ix_i]+1)/2;
	 for (i=1; i<=ppst->nsq; i++) {
	   ppst->pam2[0][i][ix_i] = ppst->pam2[0][i][ix_l] =
	     (ppst->pam2[0][i][ix_l]+ppst->pam2[0][i][ix_i]+1)/2;
	   ppst->pam2[0][ix_i][i] = ppst->pam2[0][ix_l][i] =
	     (ppst->pam2[0][ix_i][i]+ppst->pam2[0][ix_l][i]+1)/2;
	 }
       }
     }
   }

   /*
   print_pam(ppst);
   */

   /* once we have a complete pam matrix, we can calculate Lambda and K 
      for "average" sequences */
   kar_p = NULL;
   init_karlin_a(ppst, aa0_f, &kar_p);
   do_karlin_a(ppst->pam2[0], ppst, aa0_f,
	       kar_p, &m_msg->Lambda, &m_msg->K, &m_msg->H);
   ppst->pLambda = m_msg->Lambda;
   ppst->pK = m_msg->K;
   ppst->pH = m_msg->H;
   ppst->LK_set = 1;
   free(kar_p);

#if defined(FASTF) || defined(FASTS) || defined(FASTM)
   if (ppst->ext_sq_set) {
     fprintf(stderr," -S not available on [t]fast[fs]\n");
     ppst->ext_sq_set = 0;
     ppst->nsq_e = ppst->nsq;

     /* reset sascii to ignore -S, map lc */
     init_ascii(0,lascii, ppst->nsq, 0);
     validate_novel_aa(lascii, ppst->nsq, 0);
   }
#endif
}

/* alloc_pam2p creates a profile structure */
int **
alloc_pam2p(int **pam2p, int len, int nsq) {
  int i, pam2p_len;
  int *pam2pp;

  if (pam2p == NULL) {
    if ((pam2p = (int **)calloc(len,sizeof(int *)))==NULL) {
      fprintf(stderr," Cannot allocate pam2p: %d\n",len);
      return NULL;
    }

    if((pam2p[0] = (int *)calloc((nsq+1)*len,sizeof(int)))==NULL) {
      fprintf(stderr, "Cannot allocate pam2p[0]: %d\n", (nsq+1)*len);
      free(pam2p);
      return NULL;
    }
  }
  else {
    pam2p_len = (nsq+1)*len*sizeof(int);
    pam2pp = pam2p[0];
    if ((pam2pp = (int *)realloc(pam2pp,pam2p_len))==NULL) {
      fprintf(stderr,
	      "Cannot reallocate pam2p[0]: %ld\n", (nsq+1)*len*sizeof(int));
      return NULL;
    }
    memset(pam2pp,0,pam2p_len);

    if ((pam2p = (int **)realloc(pam2p,len*sizeof(int *)))==NULL) {
      fprintf(stderr," Cannot reallocate pam2p: %d\n",len);
      return NULL;
    }
    pam2p[0] = pam2pp;
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

/* sortbest has now become comparison function specific so that we can use
   a different comparison for fasts/f 
*/
#if !defined(FASTS) && !defined (FASTF) && !defined(FASTM)
void
qshuffle() {}

#ifndef LALIGN	 /* LALIGN has last_calc() in last_thresh.c */
int
last_calc(
	  unsigned char *aa0, unsigned char *aa1, int maxn,
	  struct beststr **bestp_arr, int nbest,
	  struct mngmsg m_msg, struct pstruct *ppst
	  , void **f_str
	  , void *pstat_str)
{
  return nbest;
}
#endif

/* this function is almost never called, thus a slow shell sort */
void sortbest (bptr, nbest, irelv)
struct beststr **bptr;
int nbest, irelv;
{
    int gap, i, j;
    struct beststr *tmp;

    for (gap = nbest/2; gap > 0; gap /= 2)
	for (i = gap; i < nbest; i++)
	    for (j = i - gap; j >= 0; j-= gap) {
	      if (bptr[j]->rst.score[irelv] >= bptr[j + gap]->rst.score[irelv]) break;
	      tmp = bptr[j];
	      bptr[j] = bptr[j + gap];
	      bptr[j + gap] = tmp;
	    }
}

void show_aux(FILE *fp, struct beststr *bptr) {}
void header_aux(FILE *fp) {}

#else
/* this function is almost never called, thus a slow shell sort */
void sortbest (bptr, nbest, irelv)
struct beststr **bptr;
int nbest, irelv;
{
    int gap, i, j;
    struct beststr *tmp;

    for (gap = nbest/2; gap > 0; gap /= 2)
	for (i = gap; i < nbest; i++)
	    for (j = i - gap; j >= 0; j-= gap) {
	      if (bptr[j]->rst.escore < bptr[j + gap]->rst.escore) break;
	      tmp = bptr[j];
	      bptr[j] = bptr[j + gap];
	      bptr[j + gap] = tmp;
	    }
}

#if defined(FASTS) || defined(FASTM)

/* this shuffle is for FASTS  */
/* convert ',' -> '\0', shuffle each of the substrings */
void
qshuffle(unsigned char *aa0, int n0, int nm0, void *rand_state) {

  unsigned char **aa0start, *aap, tmp;
  int i,j,k, ns;

  if ((aa0start=(unsigned char **)calloc(nm0+1,
					 sizeof(unsigned char *)))==NULL) {
    fprintf(stderr,"cannot calloc for qshuffle %d\n",nm0);
    exit(1);
  }

  aa0start[0]=aa0;
  for (k=1,i=0; i<n0; i++) {
    if (aa0[i]==EOSEQ || aa0[i]==ESS) {
      aa0[i]='\0';
      aa0start[k++] = &aa0[i+1];
    }
  }  

  /* aa0start has the beginning of each substring */
  for (k=0; k<nm0; k++) {
    aap=aa0start[k];
    ns = strlen((const char *)aap);
    for (i=ns; i>1; i--) {
      j = my_nrand(i, rand_state);
      tmp = aap[j];
      aap[j] = aap[i-1];
      aap[i-1] = tmp;
    }
    aap[ns] = 0;
  }

  for (k=1; k<nm0; k++) {
/*  aap = aa0start[k];
    while (*aap) fputc(pst.sq[*aap++],stderr);
    fputc('\n',stderr);
*/
    aa0start[k][-1]=ESS;
  }

  free(aa0start);
}
#endif

#ifdef FASTF
void qshuffle(unsigned char *aa0, int n0, int nm0, void *rand_state) {

  int i, j, k, nmpos;
  unsigned char tmp;
  int nmoff;
  
  nmoff = (n0 - nm0 - 1)/nm0 + 1;

  for (i = nmoff-1 ; i > 0 ; i--) {

    /* j = nrand(i); if (i == j) continue;*/       /* shuffle columns */ 
    j = (nmoff -1 ) - i; 
    if (i <= j) break; /* reverse columns */

    /* swap all i'th column residues for all j'th column residues */
    for(nmpos = 0, k = 0 ; k < nm0 ; k++, nmpos += nmoff+1 ) {
      tmp = aa0[nmpos + i];
      aa0[nmpos + i] = aa0[nmpos + j];
      aa0[nmpos + j] = tmp;
    }
  }
}
#endif


/* show additional best_str values */
void show_aux(FILE *fp, struct beststr *bptr) {
  fprintf(fp," %2d %3d",bptr->rst.segnum,bptr->rst.seglen);
}

void header_aux(FILE *fp) {
  fprintf(fp, " sn  sl");
}
#endif

void
fill_pam(int **pam2p, int n0, int nsq, double **freq2d, double scale, int **no_remap) {
  int i, j, new_j;
  double freq;

  /* fprintf(stderr, "scale: %g\n", scale); */
  
  /* now fill in the pam matrix: */
  for (j = 1 ; j <=20 ; j++) {
    new_j = qascii[pssm_aa[j]];
    for (i = 0 ; i < n0 ; i++) {
      freq = scale * freq2d[i][j-1];
      if ( freq < 0.0) freq -= 0.5;
      else freq += 0.5;

      if (no_remap[i][j-1]) {
	pam2p[i][j] = (int)freq;
      }
      else {
	pam2p[i][new_j] = (int)(freq);
      }
    }
  }
}

double
get_lambda(int **pam2p, int n0, int nsq, unsigned char *query) {
  double lambda, H;
  double *pr, tot, sum;
  int i, ioff, j, min, max, q_i;

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
    for (j = 1; j < nsq ; j++) {
      if (min > pam2p[i][j])
	min = pam2p[i][j];
      if (max < pam2p[i][j])
	max = pam2p[i][j];
    }
  }

  /*  fprintf(stderr, "min: %d\tmax:%d\n", min, max); */
  
  if ((pr = (double *) calloc(max - min + 1, sizeof(double))) == NULL) {
    fprintf(stderr, "Couldn't allocate memory for score probabilities: %d\n", max - min + 1);
    exit(1);
  }

  tot = (double) rrtotal * (double) rrtotal * (double) n0;
  for (i = ioff ; i < n0 ; i++) {

    if (query[i] < 'A') {q_i = query[i];}
    else {q_i= aascii[query[i]];}

    for (j = 1; j < nsq ; j++) {
      pr[pam2p[i][j] - min] +=
	(double) ((double) rrcounts[q_i] * (double) rrcounts[j]) / tot;
    }
  }

  sum = 0.0;
  for(i = 0 ; i <= max-min ; i++) { 
    sum += pr[i];
    /*     fprintf(stderr, "%3d: %g %g\n", i+min, pr[i], sum); */
  }
  /*   fprintf(stderr, "sum: %g\n", sum); */

  for(i = 0 ; i <= max-min ; i++) { pr[i] /= sum; }

  if (!karlin(min, max, pr, &lambda, &H)) {
    fprintf(stderr, "Karlin lambda estimation failed\n");
  }

  /*   fprintf(stderr, "lambda: %g\n", lambda); */
  free(pr);

  return lambda;
}

/*
   *aa0 - query sequence
   n0   - length
   pamscale - scaling for pam matrix - provided by apam.c, either
              0.346574 = ln(2)/2 (P120, BL62) or
	      0.231049 = ln(2)/3 (P250, BL50) 
*/

void
scale_pssm(int **pssm2p, double **freq2d,
	   unsigned char *query, int n0,
	   int **pam2, double pamscale);

static unsigned char ustandard_aa[] ="\0ARNDCQEGHILKMFPSTWYV";

void
read_pssm(unsigned char *aa0, int n0, int nsq,
	  double pamscale, 
	  FILE *fp, int pgpf_type, struct pstruct *ppst) {
  int i, j, len, k;
  int qi, rj;	/* qi - index query; rj - index residues (1-20) */
  int **pam2p;
  int first, too_high;
  unsigned char *query, ctmp;
  char dline[512];
  double freq, **freq2d, lambda, new_lambda;
  double scale, scale_high, scale_low;

  pam2p = ppst->pam2p[0];

  if (pgpf_type == 0) {

    if (1 != fread(&len, sizeof(int), 1, fp)) {
      fprintf(stderr, "error reading from checkpoint file: %d\n", len);
      exit(1);
    }

    if (len != n0) {
      fprintf(stderr, "profile length (%d) and query length (%d) don't match!\n",
	      len,n0);
      exit(1);
    }

    /* read over query sequence stored in BLAST profile */
    if(NULL == (query = (unsigned char *) calloc(len+2, sizeof(char)))) {
      fprintf(stderr, "Couldn't allocate memory for query!\n");
      exit(1);
    }

    if (len != fread(query, sizeof(char), len, fp)) {
      fprintf(stderr, "Couldn't read query sequence from profile: %s\n", query);
      exit(1);
    }
  }
  else if (pgpf_type == 1) {

    if ((fgets(dline,sizeof(dline),fp) == NULL)  ||
	(1 != sscanf(dline, "%d",&len))) {
      fprintf(stderr, "error reading from checkpoint file: %d\n", len);
      exit(1);
    }

    if(len != n0) {
      fprintf(stderr, "profile length (%d) and query length (%d) don't match!\n",
	      len,n0);
      exit(1);
    }

    /* read over query sequence stored in BLAST profile */
    if(NULL == (query = (unsigned char *) calloc(len+2, sizeof(char)))) {
      fprintf(stderr, "Couldn't allocate memory for query!\n");
      exit(1);
    }

    if (fgets((char *)query,len+2,fp)==NULL) {
      fprintf(stderr, "Couldn't read query sequence from profile: %s\n", query);
      exit(1);
    }
  }  
  else {
    fprintf(stderr," Unrecognized PSSM file type: %d\n",pgpf_type);
    exit(1);
  }

  /* currently we don't do anything with query; ideally, we should
     check to see that it actually matches aa0 ... */

  /* quick 2d array alloc: */
  if((freq2d = (double **) calloc(n0, sizeof(double *))) == NULL) {
    fprintf(stderr, "Couldn't allocate memory for frequencies!\n");
    exit(1);
  }

  if((freq2d[0] = (double *) calloc(n0 * 20, sizeof(double))) == NULL) {
    fprintf(stderr, "Couldn't allocate memory for frequencies!\n");
    exit(1);
  }

  /* a little pointer arithmetic to fill out 2d array: */
  for (i = 1 ; i < n0 ; i++) {
    freq2d[i] = freq2d[i-1] + 20;
  }

  if (pgpf_type == 0) {
    for (qi = 0 ; qi < n0 ; qi++) {
      for (rj = 0 ; rj < 20 ; rj++) {
	if(1 != fread(&freq, sizeof(double), 1, fp)) {
	  fprintf(stderr, "Error while reading frequencies!\n");
	  exit(1);
	}
	freq2d[qi][rj] = freq;
      }
    }
  }
  else {
    for (qi = 0 ; qi < n0 ; qi++) {
      if ((fgets(dline,sizeof(dline),fp) ==NULL) ||
      (k = sscanf(dline,"%c %lg %lg %lg %lg %lg %lg %lg %lg %lg %lg %lg %lg %lg %lg %lg %lg %lg %lg %lg %lg\n",
		 &ctmp, &freq2d[qi][0], &freq2d[qi][1], &freq2d[qi][2], &freq2d[qi][3], &freq2d[qi][4], 
		 &freq2d[qi][5], &freq2d[qi][6], &freq2d[qi][7], &freq2d[qi][8], &freq2d[qi][9],
		 &freq2d[qi][10], &freq2d[qi][11], &freq2d[qi][12], &freq2d[qi][13], &freq2d[qi][14],
		      &freq2d[qi][15], &freq2d[qi][16], &freq2d[qi][17], &freq2d[qi][18], &freq2d[qi][19]))<1) {
	fprintf(stderr, "Error while reading frequencies: %d read!\n",k);
	exit(1);
      }
      for (rj=0; rj<20; rj++) { freq2d[qi][rj] /= 10.0; }	/* reverse scaling */
    }
  }

  scale_pssm(ppst->pam2p[0], freq2d, query, n0, ppst->pam2[0],pamscale);

  free(freq2d[0]);
  free(freq2d);

  free(query);
}

/* before fasta-36.3.6 (with reordered amino-acid mapping), scale_pssm()
   simply produced a log(q_ij/p_j) and put it into pam2p.

   But pssm's use pssm_aa encoding, while fasta-36.3.6 use NCBIstdaa
   encoding, so the pam2p must be re-mapped
*/

void
scale_pssm(int **pssm2p, double **freq2d, unsigned char *query, int n0, int **pam2, double pamscale) {
  int i, qi, rj;
  double freq, new_lambda, lambda;
  int first, too_high;
  double scale, scale_high, scale_low;
  int **no_remap;


  /* quick 2d array alloc: */
  if((no_remap = (int **) calloc(n0, sizeof(int *))) == NULL) {
    fprintf(stderr, "***error [%s:%d] Couldn't allocate memory for remap[%d]\n",__FILE__, __LINE__, n0);
    exit(1);
  }

  if((no_remap[0] = (int *) calloc(n0 * 20, sizeof(int))) == NULL) {
    fprintf(stderr, "***error [%s:%d] Couldn't allocate memory for remap[%d]\n",__FILE__, __LINE__, n0);
    exit(1);
  }

  for (qi=1; qi < n0; qi++) {
    no_remap[qi] = no_remap[qi-1]+20;
  }

  /* convert freq2d from frequences to log_scores; 
     fill zeros with BLOSUM62 values */

  for (rj = 0 ; rj < 20 ; rj++) {
    for (qi = 0 ; qi < n0 ; qi++) {
      if (freq2d[qi][rj] > 1e-20) {
	freq = log(freq2d[qi][rj] /((double) (rrcounts[rj+1])/(double) rrtotal));
	freq /= pamscale; /* this gets us close to originial pam scores */
	freq2d[qi][rj] = freq;
      }
      else {		
	/* when blastpgp decides to leave something out, it puts 0's in all the frequencies
	   in the binary checkpoint file.  In the ascii version, however, it uses BLOSUM62
	   values.  I will put in scoring matrix values as well */
	/* 11-Oct-2015 -- this does not work properly, because the
	   correct amino-acid ordering is not used -- pam2 uses
	   NCBIStdaa ordering, but the rest of the matrix uses pssm_aa
	   ordering, which is changed in fill_pam */

	no_remap[qi][rj] = 1;
	if (query[qi] < 'A') {
	  freq2d[qi][rj] = pam2[query[qi]][rj+1];
	}
	else {
	  freq2d[qi][rj] = pam2[aascii[query[qi]]][rj+1];
	}
      }
    }
  }

  /* now figure out the right scale */
  scale = 1.0;
  lambda = get_lambda(pam2, 20, 20, ustandard_aa);

#ifdef DEBUG
  /*
  fill_pam(pssm2p, n0, 20, freq2d, scale, no_remap);
  fprintf(stderr,"        ");
  for (rj = 1; rj <= 20; rj++) {
    fprintf(stderr,"   %c", NCBIstdaa[rj]);
  }
  fprintf(stderr,"\n");
  for (qi = 0 ; qi < n0 ; qi++) {
    fprintf(stderr, "%4d %c: ", qi+1, NCBIstdaa[query[qi]]);
    for (rj = 1 ; rj <= 20 ; rj++) {
      fprintf(stderr, "%4d", pssm2p[qi][rj]);
    }
    fprintf(stderr, "\n");
  }
  */
#endif

  /* should be near 1.0 because of our initial scaling by ppst->pamscale */
  /* fprintf(stderr, "real_lambda: %g\n", lambda); */

  /* get initial high/low scale values: */
  first = 1;
  while (1) {
    fill_pam(pssm2p, n0, 20, freq2d, scale, no_remap);
    new_lambda = get_lambda(pssm2p, n0, 20, query); 

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
    fill_pam(pssm2p, n0, 20, freq2d, scale, no_remap);
    new_lambda = get_lambda(pssm2p, n0, 20, query);
    
    if (new_lambda > lambda) scale_low = scale;
    else scale_high = scale;
  }

  scale = 0.5 * (scale_high + scale_low);
  fill_pam(pssm2p, n0, 20, freq2d, scale, no_remap);

  free(no_remap[0]);
  free(no_remap);

#ifdef DEBUG
  /*
    fprintf(stderr, "final scale: %g\n", scale);

    fprintf(stderr,"        ");
    for (rj = 1; rj <= 20; rj++) {
      fprintf(stderr,"   %c", NCBIstdaa[rj]);
    }
    fprintf(stderr,"\n");
    for (qi = 0 ; qi < n0 ; qi++) {
      fprintf(stderr, "%4d %c: ", qi+1, NCBIstdaa[query[qi]]);
      for (rj = 1 ; rj <= 20 ; rj++) {
	fprintf(stderr, "%4d", pssm2p[qi][rj]);
      }
      fprintf(stderr, "\n");
    }
  */
#endif

}

#if defined(CAN_PSSM)
int
parse_pssm_asn_fa(FILE *afd, int *n_rows, int *n_cols,
		  unsigned char **query,
		  double ***wfreqs, double ***freqs, int ***iscores,
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
  unsigned char *query, ctmp;
  char dline[512];
  char matrix[MAX_SSTR];
  double psi2_lambda;
  double freq, **wfreq2d=NULL, **freq2d=NULL, lambda, new_lambda;
  double scale, scale_high, scale_low;
  int **iscores2d=NULL;
  int gap_open, gap_extend;
  int n_rows, n_cols;


  pam2p = ppst->pam2p[0];

  /* get the information from the ASN.1 (binary) file */
  if (parse_pssm_asn_fa(fp, &n_rows, &n_cols, &query, &wfreq2d, &freq2d, &iscores2d,
			matrix, &gap_open, &gap_extend, &psi2_lambda)<=0) {
    return -1;
  }

  /* not using wfreq2d[][] right now, free it */
  if (wfreq2d != NULL) {
    if (wfreq2d[0] != NULL) {free(wfreq2d[0]);}
    free(wfreq2d);
  }

  /* do we have a query sequence */
  if (query == NULL) { query = aa0;}

  if (!gap_set) {
    if (gap_open) {
      if (gap_open > 0) {gap_open = -gap_open;}
      ppst->gdelval = gap_open;
    }
    else if (strncmp(matrix,"BLOSUM62",8)==0) {
      ppst->gdelval = -11;
    }
    gap_set = 1;
  }
  if (!del_set) {
    if (gap_extend) {
      if (gap_extend > 0) {gap_extend = -gap_extend;}
      ppst->ggapval = gap_extend;
    }
    else if (strncmp(matrix,"BLOSUM62",8)==0) {
      ppst->ggapval = -1;
    }
    del_set = 1;
  }

  if (strncmp(matrix, "BLOSUM62", 8)== 0 && !ppst->pam_set) {
    SAFE_STRNCPY(ppst->pamfile, "BL62", 120);
    SAFE_STRNCPY(ppst->pamfile_save, ppst->pamfile, 120);
    standard_pam(ppst->pamfile,ppst,del_set, gap_set);
    if (!ppst->have_pam2) {
     alloc_pam (MAXSQ, MAXSQ, ppst);
    }
    init_pam2(ppst);
    ppst->pam_set = 1;
  }

  if (n_cols < n0) { 
    fprintf(stderr, " query length: %d != n_cols: %d\n",n0, n_cols);
    exit(1);
  }

  /* try to just use the the iscore2d file */
  if (iscores2d != NULL) {
    for (qi = 0 ; qi < n0 ; qi++) {
      for (rj = 1 ; rj <= 24 ; rj++) {
	itmp = iscores2d[qi][rj];
	if (itmp < -256) itmp=0;
	pam2p[qi][rj] = itmp;
      }
    }
    /* all done, free it */
    free(iscores2d[0]);
    free(iscores2d);
  }
  else {
    scale_pssm(ppst->pam2p[0], freq2d, query, n0, ppst->pam2[0], pamscale);
  }

#if DEBUG
  if (ppst->debug_lib) {
    /*    fprintf(stderr, "final scale: %g\n", scale); */

    fprintf(stderr,"        ");
    for (rj = 1; rj <= 24; rj++) {
      fprintf(stderr,"  %c", NCBIstdaa[rj]);
    }
    fprintf(stderr,"\n");
    for (qi = 0 ; qi < n0 ; qi++) {
      fprintf(stderr, "%3d %c: ", qi+1, NCBIstdaa[aa0[qi]]);
      for (rj = 1 ; rj <= 24 ; rj++) {
	fprintf(stderr, "%3d", pam2p[qi][rj]);
      }
      fprintf(stderr, "\n");
    }
  }
#endif
  
  if (freq2d != NULL) {
    free(freq2d[0]);
    free(freq2d);
  }

  if (query != aa0) free(query);
  return 1;
}
#endif

/* last_params() sets up values in pstruct *ppst now that all
   parameters and data is available.

   It:
   (1) moves m_msg->n0 to ppst->n0 for statistics calculations
   (2) sets ppst->nsq_e
   (3) reads the PSSM file if one is being used
   (4) calculates m_msg->nm0 for FASTF/S/M
   (5) determines statistical strategy for FASTF/S, sets last_calc_flg
       and qshuffle
   (6) lowers ktup for short sequences
*/

void
last_params(unsigned char *aa0, int n0, 
	    struct mngmsg *m_msp,
	    struct pstruct *ppst
	    ) {
  int i, nsq;
  FILE *fp;
  int is_fastxy=0;
  int n0_eff;
  /* do_karlin_a() must be re-run everytime the scoring matrix changes */
  double *kar_p;
  double aa0_f[MAXSQ];

  if (n0 < 0) { return;}

  
  n0_eff = m_msp->n0;
  ppst->n0 = m_msp->n0;
#if !defined(TFAST) && (defined(FASTX) || defined(FASTY))
  n0_eff /= 3;
  is_fastxy = 1;
#endif

  /* reset the PAMFILE to the original value */
  if (strncmp(ppst->pamfile, ppst->pamfile_save,120)!=0) {
    SAFE_STRNCPY(ppst->pamfile, ppst->pamfile_save, 120);
    standard_pam(ppst->pamfile,ppst,del_set, gap_set);
    init_pam2(ppst);
    init_pamx(ppst);
  }

  /* **************************************************************** */
  /* adjust scoring matrix for short protein/translated protein queries */
  /* **************************************************************** */

  if (ppst->pam_variable) {
    if (min_pam_bits(n0_eff, DEF_MIN_BITS, ppst, del_set, gap_set)) {
      init_pam2(ppst);
      init_pamx(ppst);
      kar_p = NULL;
      init_karlin_a(ppst, aa0_f, &kar_p);
      do_karlin_a(ppst->pam2[0], ppst, aa0_f,
		  kar_p, &m_msp->Lambda, &m_msp->K, &m_msp->H);
      ppst->pLambda = m_msp->Lambda;
      ppst->pK = m_msp->K;
      ppst->pH = m_msp->H;
      ppst->LK_set = 1;
      free(kar_p);
    }
    else {
      fprintf(stderr,"+++ warning [%s:%d] - query too short [%d] for %d bit signal -- fasts36 may be more useful +++\n",
	      __FILE__, __LINE__, n0, DEF_MIN_BITS);
    }
  }

  if (ppst->ext_sq_set) { ppst->nsq_e = nsq = 2*ppst->nsq; }
  else {ppst->nsq_e = nsq = ppst->nsq;}

#if defined(CAN_PSSM)
  ppst->pam2p[0] = alloc_pam2p(ppst->pam2p[0],n0,MAXSQ);
  ppst->pam2p[1] = alloc_pam2p(ppst->pam2p[1],n0,MAXSQ);

  if (ppst->pam_pssm) {
    if ((ppst->pgpfile_type == 0) && (fp=fopen(ppst->pgpfile,"rb"))) {
      read_pssm(aa0, n0, ppst->nsq, ppst->pamscale, fp, 0, ppst);
      extend_pssm(aa0, n0, ppst);
    }
    else if ((ppst->pgpfile_type == 1) && (fp=fopen(ppst->pgpfile,"r"))) {
      read_pssm(aa0, n0, ppst->nsq, ppst->pamscale, fp, 1, ppst);
      extend_pssm(aa0, n0, ppst);
    }
    else if ((ppst->pgpfile_type == 2) && (fp=fopen(ppst->pgpfile,"rb"))) {
      if (read_asn_pssm(aa0, n0, ppst->nsq, ppst->pamscale, fp, ppst)>0) {
	extend_pssm(aa0, n0, ppst);
      }
      else {
	fprintf(stderr," Could not parse PSSM file: %s\n",ppst->pgpfile);
	ppst->pam_pssm = 0;
	return;
      }
    }
    else {
      fprintf(stderr," Could not open PSSM file: %s\n",ppst->pgpfile);
      ppst->pam_pssm = 0;
      return;
    }
  }
#endif

#if defined(FASTF) || defined(FASTS) || defined(FASTM)
  m_msp->nm0 = 1;
  for (i=0; i<n0; i++)
    if (aa0[i]==EOSEQ || aa0[i]==ESS) m_msp->nm0++;

/*
  for FASTS, we can do statistics in one of two different ways
  if there are <= 10 query fragments, then we calculate probabilistic
  scores for every library sequence.  If there are > 10 fragments, this
  takes much too long and too much memory, so we use the old fashioned
  raw score only z-score normalized method initially, and then calculate
  the probabilistic scores for the best hits.  To scale those scores, we
  also need a set of random probabilistic scores.  So we do the qshuffle
  to get them.

  For FASTF, precalculating probabilities is prohibitively expensive,
  so we never do it; FASTF always acts like FASTS with nfrags>10.

*/

#if defined(FASTS) || defined(FASTM)
  if (m_msp->nm0 > 10) m_msp->escore_flg = 0;
  else m_msp->escore_flg = 1;
#endif

  if (m_msp->escore_flg && (ppst->zsflag&1)) {
    m_msp->last_calc_flg = 0;
    m_msp->qshuffle = 0;
  }
  else {	/* need random query, second set of 2000 scores */
    m_msp->last_calc_flg = 1;
    m_msp->qshuffle = 1;
  }
#else
#ifndef LALIGN
  m_msp->last_calc_flg = 0;
#else
  m_msp->last_calc_flg = 1;	/* LALIGN needs last_calc for threshold */
#endif
  m_msp->qshuffle = 0;
  m_msp->escore_flg = 0;
  m_msp->nm0 = 1;
#endif

/* adjust the ktup if appropriate */  

  if (pgm_def_arr[ppst->pgm_id].ktup > 0) {
    if (!ktup_set ) {
      ppst->param_u.fa.ktup = pgm_def_arr[ppst->pgm_id].ktup;
      if (m_msp->qdnaseq == SEQT_PROT || is_fastxy) {
#if defined(FASTS) || defined(FASTM)
	if (n0_eff > 100 && ppst->param_u.fa.ktup > 2) ppst->param_u.fa.ktup = 2;
#endif
	if (n0_eff <= 40 && ppst->param_u.fa.ktup > 1) ppst->param_u.fa.ktup = 1;
      }
      else if (m_msp->qdnaseq == SEQT_DNA || m_msp->qdnaseq == SEQT_RNA) {
	if (n0_eff <= 20 && ppst->param_u.fa.ktup > 1) ppst->param_u.fa.ktup = 1;
#if defined(FASTS) || defined(FASTM)
	/* with the current (April 12 2005) dropfs2.c - ktup cannot be > 2 */
	else ppst->param_u.fa.ktup = 2;
#else
	else if (n0 < 50 && ppst->param_u.fa.ktup > 2) ppst->param_u.fa.ktup = 2;
	else if (n0 < 100 && ppst->param_u.fa.ktup > 3)  ppst->param_u.fa.ktup = 3;
#endif
      }
    }
    /* regardless of ktup state */
    if (ppst->param_u.fa.use_E_thresholds) {
      ppst->param_u.fa.use_E_thresholds = ppst->LK_set;
    }
    if (!E_cgap_set) {
      ppst->param_u.fa.E_join = ppst->param_u.fa.E_band_opt * 5;
    }
    else {
      if (ppst->param_u.fa.E_join > 1.0) {
	ppst->param_u.fa.E_join = ppst->param_u.fa.E_band_opt * ppst->param_u.fa.E_join;
      }
    }
  }
}

/* validate_params() exists because of bugs that keep appearing

   (1) pam2[0][x][0] or pam2[0][0][x] are not -BIGNUM
   (2) sascii[] (or qascii[], lascii[]) have values outside nsq_e.
 */

int
validate_params(const unsigned char *aa0, int n0, 
		const struct mngmsg *m_msg,
		const struct pstruct *ppst,
		const int *lascii, const int *pascii) {
  int good_params = 1;
  int i;

  /* check for -BIGNUM for boundaries of pam2[0][0:x][x:0] */

  for (i=0; i< ppst->nsq; i++) {
    if (ppst->pam2[0][0][i] > -1000) {
      fprintf(stderr," *** ERROR ***  pam2[0][0][%d/%c] == %d\n",
	      i,NCBIstdaa[i],ppst->pam2[0][0][i]);
      good_params = 0;
    }
    if (ppst->pam2[0][i][0] > -1000) {
      fprintf(stderr," *** ERROR ***  pam2[0][%d/%c][0] == %d\n",
	      i,NCBIstdaa[i],ppst->pam2[0][i][0]);
      good_params = 0;
    }
  }

  /* check for -BIGNUM for boundaries of pam2[1][0:x][x:0] */
  if (ppst->ext_sq_set) {
    for (i=0; i< ppst->nsqx; i++) {
      if (ppst->pam2[1][0][i] > -1000) {
	fprintf(stderr," *** ERROR ***  pam2[1][0][%d] == %d\n",
		i,ppst->pam2[1][0][i]);
      good_params = 0;
      }
      if (ppst->pam2[1][i][0] > -1000) {
	fprintf(stderr," *** ERROR ***  pam2[1][%d][0] == %d\n",
		i,ppst->pam2[1][i][0]);
      good_params = 0;
      }
    }
  }

  /* check for valid residues in query */
  for (i=0; i<n0; i++) {
    if (aa0[i] > ppst->nsq_e && aa0[i] != ESS) {
      fprintf(stderr," *** ERROR *** aa0[%d] = %c[%d > %d] out of range\n",
	      i, aa0[i], aa0[i], ppst->nsq_e);
      good_params = 0;
    }
  }

  for (i=0; i<128; i++) {
    if (lascii[i] < NA && lascii[i] > ppst->nsq_e) {
      fprintf(stderr," *** ERROR *** lascii [%c|%d] = %d > %d out of range\n",
	      i, i, lascii[i], ppst->nsq_e);
      good_params = 0;
    }

    /*  currently, pascii[] is not reset for upper-case only
    if (pascii[i] < NA && pascii[i] > ppst->nsq_e) {
      fprintf(stderr," *** WARNING *** pascii[%c|%d] = %d > %d out of range\n",
	      i, i, pascii[i], ppst->nsq_e);
    }
    */

  }

  return good_params;
}

/* given a good profile in ppst->pam2p[0], make an extended profile
   in ppst->pam2p[1]
*/
void
extend_pssm(unsigned char *aa0, int n0, struct pstruct *ppst) {

  int i, j, nsq;
  int sa_x, sa_t, sa_b, sa_z, sa_j;
  int **pam2p0, **pam2p1;

  nsq = ppst->nsq;

  pam2p0 = ppst->pam2p[0];
  pam2p1 = ppst->pam2p[1];

  sa_x = pascii['X'];
  sa_t = pascii['*'];
  if (sa_t >= ppst->nsq) {sa_t = sa_x;}
  sa_b = pascii['B'];
  sa_z = pascii['Z'];
  sa_j = pascii['J'];

  /* fill in boundaries, B, Z, *, X */
  for (i=0; i<n0; i++) {
    pam2p0[i][0] = -BIGNUM;
    pam2p0[i][sa_b] = (int)
      (((float)pam2p0[i][pascii['N']]+(float)pam2p0[i][pascii['D']]+0.5)/2.0);
    pam2p0[i][sa_z] = (int)
      (((float)pam2p0[i][pascii['Q']]+(float)pam2p0[i][pascii['E']]+0.5)/2.0);
    pam2p0[i][sa_j] = (int)
      (((float)pam2p0[i][pascii['I']]+(float)pam2p0[i][pascii['L']]+0.5)/2.0);
    pam2p0[i][sa_x] = ppst->pam_xm;
    pam2p0[i][sa_t] = ppst->pam_xm;
  }

  /* copy pam2p0 into pam2p1 */
  for (i=0; i<n0; i++) {
    pam2p1[i][0] = -BIGNUM;
    for (j=1; j<=ppst->nsq; j++) {
      pam2p1[i][j] = pam2p0[i][j];
    }
  }

  /* then fill in extended characters, if necessary */
  if (ppst->ext_sq_set) {
    for (i=0; i<n0; i++) {
      for (j=1; j<=ppst->nsq; j++) {
	pam2p0[i][nsq+j] = pam2p0[i][j];
	pam2p1[i][nsq+j] = ppst->pam_xm;
      }
    }
  }
}

void format_params(struct opt_def_str *opt_ptr, char *string) {
  
  if (opt_ptr->opt_char == 'r') {
    sprintf(string, " [+%d/%d]", opt_ptr->i_param1, opt_ptr->i_param2);
    return;
  }

  switch (opt_ptr->fmt_type) {

  case 1:
    sprintf(string, " [%d]", opt_ptr->i_param1); break;
  case 2:
    sprintf(string, " [%d,%d]", opt_ptr->i_param1, opt_ptr->i_param2); break;
  case 3:
    sprintf(string, " [%.4g]", opt_ptr->d_param1); break;
  case 4:
    sprintf(string, " [%.4g,%.4g]", opt_ptr->d_param1, opt_ptr->d_param2); break;
  case 5:
    sprintf(string, " [%s]", opt_ptr->s_param); break;
  case 0:
  default:
    string[0] = '\0'; break;
  }
}


#if defined(FASTX) || defined(FASTY)
static char *common_opts = "sfgjSEbdI";
#else
#if defined(LALIGN)
static char *common_opts = "sfgEZI";
#else
static char *common_opts = "sfgSbdI";
#endif
#endif

void
show_help(char *pgm_name, int pgm_id) {
  int i, j;
  int opt_line_cnt=0;
  char tmp_string[MAX_STR];
  struct opt_def_str *opt_ptr;

  printf("USAGE\n");
#ifndef LALIGN
  printf(" %s [-options] query_file library_file",pgm_name);
  if (pgm_def_arr[pgm_id].ktup > 0) {
    printf(" [ktup]\n");
  }
  else {printf("\n");}
#else
  printf(" %s [-options] seq_file1 seq_file2\n",pgm_name);
#endif
  printf(" %s -help for a complete option list\n",pgm_name);
  printf("\nDESCRIPTION\n");
  printf(" %s\n version: %s\n",pgm_def_arr[pgm_id].iprompt0, verstr);
  printf("\n");
  printf("COMMON OPTIONS (options must preceed query_file library_file)\n");

  for (i=0; i<strlen(common_opts); i++) {
    opt_ptr = g_options;
    for (j=0; opt_ptr[j].opt_char != '\0'; j++) {
      if (common_opts[i]==opt_ptr[j].opt_char) {
	format_params(&opt_ptr[j], tmp_string);
	printf(" -%c%c %s %s;",opt_ptr[j].opt_char, (opt_ptr[j].has_arg? ':' : ' '), 
	       tmp_string, opt_ptr[j].opt_descr_s);
	/* if ((++opt_line_cnt % 2) == 0) printf("\n"); */
	printf("\n");
	goto next_option;
      }
    }

    opt_ptr = f_options;
    for (j=0; opt_ptr[j].opt_char != '\0'; j++) {
      if (common_opts[i]==opt_ptr[j].opt_char) {
	format_params(&opt_ptr[j], tmp_string);
	printf(" -%c%c %s %s;",opt_ptr[j].opt_char, (opt_ptr[j].has_arg? ':' : ' '), 
	       tmp_string, opt_ptr[j].opt_descr_s);
	/* if ((++opt_line_cnt % 2)==0) printf("\n"); */
	printf("\n");
      }
    }
  next_option: continue;
  }
  if ((opt_line_cnt % 2) != 0) printf("\n");
  exit(0);
}

/* sorts a list of options, with upper and lower case characters
   sorted together */
void sort_opt_list(char *v, int n) {
  int gap, i, j, k;
  int incs[7] = { 336, 112, 48, 21, 7, 3, 1 };
  char tmp_c, tmp_u;
  int v_start;

  /* first shell sort the list using toupper() */
  for ( k = 0; k < 7; k++) {
    gap = incs[k];
    for (i = gap; i < n; i++) {
      tmp_c = v[i];
      tmp_u = toupper(tmp_c);
      j = i;
      while (j >= gap && toupper(v[j - gap]) > tmp_u) {
	v[j] = v[j - gap];
	j -= gap;
      }
      v[j] = tmp_c;
    }
  }
  /* then sort the toupper(==) pairs lower-case, upper-case */

  for (i=1; i<n; i++) {
    if (toupper(v[i])==toupper(v[i-1])) {
      if (v[i] > v[i-1]) { tmp_c = v[i]; v[i] = v[i-1]; v[i-1]=tmp_c;}
    }
  }
}

char *
sort_options (struct opt_def_str *g_options, struct opt_def_str *f_options) {
  struct opt_def_str *this_option;
  char *sorted_list, *sort_ptr;
  int i, opt_count;

  opt_count=0;
  this_option = g_options;
  while ((this_option++)->opt_char!='\0') { opt_count++;}
  this_option = f_options;
  while ((this_option++)->opt_char!='\0') { opt_count++;}

  if ((sorted_list = (char *)calloc(opt_count+1, sizeof(char)))==NULL) {
    fprintf(stderr," cannot allocate sorted_list[%d]\n",opt_count+1);
    exit(1);
  }

  sort_ptr = sorted_list;
  this_option = g_options;
  while (this_option->opt_char!='\0') {
    *sort_ptr++ = this_option->opt_char;
    this_option++;
  }

  this_option = f_options;
  while (this_option->opt_char!='\0') {
    *sort_ptr++ = this_option->opt_char;
    this_option++;
  }

  sort_opt_list(sorted_list, opt_count);

  return sorted_list;
}

void
show_all_help(char *pgm_name, int pgm_id) {
  int i, j;
  struct opt_def_str *opt_ptr;
  char tmp_string[MAX_STR];
  char *descr_ptr;
  char *sorted_list;
  
  sorted_list = sort_options(g_options, f_options);

  printf("USAGE\n");
  printf(" %s [-options] query_file library_file",pgm_name);
  if (pgm_def_arr[pgm_id].ktup > 0) {
    printf(" [ktup]\n");
  }
  else {printf("\n");}

  printf(" \"@\" query_file uses stdin; query_file:begin-end sets subset range\n");
  printf(" library file formats: 0:FASTA; 1:GenBankFF; 3:EMBL_FF; 7:FASTQ; 10:subset; 12:NCBI blastdbcmd;\n");
  printf("   alternate library formats: \"library_file 7\" for 7:FASTQ\n");

  printf("\nDESCRIPTION\n");
  printf(" %s\n version: %s\n",pgm_def_arr[pgm_id].iprompt0, verstr);
  printf("\n");
  printf("OPTIONS (options must preceed query_file library_file)\n");

  for (i=0; i<strlen(sorted_list); i++) {
    opt_ptr = g_options;
    for (j=0; opt_ptr[j].opt_char != '\0'; j++) {
      if (sorted_list[i]==opt_ptr[j].opt_char) {
	descr_ptr = (opt_ptr[j].opt_descr_l) ? opt_ptr[j].opt_descr_l : opt_ptr[j].opt_descr_s;
	format_params(&opt_ptr[j], tmp_string);
	printf(" -%c%c %s %s\n",opt_ptr[j].opt_char, (opt_ptr[j].has_arg? ':' : ' '),tmp_string, descr_ptr);
	goto next_option;
      }
    }

    opt_ptr = f_options;
    for (j=0; opt_ptr[j].opt_char != '\0'; j++) {
      if (sorted_list[i]==opt_ptr[j].opt_char) {
	descr_ptr = (opt_ptr[j].opt_descr_l) ? opt_ptr[j].opt_descr_l : opt_ptr[j].opt_descr_s;
	format_params(&opt_ptr[j], tmp_string);
	printf(" -%c%c %s %s\n",opt_ptr[j].opt_char, (opt_ptr[j].has_arg? ':' : ' '), tmp_string, descr_ptr);
      }
    }
  next_option: continue;
  }

  free(sorted_list);

  exit(0);
}
