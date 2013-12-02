/* Concurrent read version */

/* $Id: msg.h 625 2011-03-23 17:21:38Z wrp $ */
/* $Revision: 625 $  */

/* Cube definitions */

#ifdef PVM_SRC
#define FIRSTNODE	1
#define FIRSTWORK	1
#else
#define FIRSTNODE	1
#define FIRSTWORK	1
#endif

#define MAXNOD		128
#define ALLTYPES        -1
#define HOSTPID		0
#define MANAGEPID 	0
#define WORKPID 	0

#define MANAGER		0
#define ALLNODES        -1
#define ALLPIDS         -1

#define STARTTYPE0 	0	/* configuration buffer values */
#define STARTTYPE1	1	/* struct mngmsg m_msp */
#define STARTTYPE2	2	/* struct pstruct ppst */
#define STARTTYPE3	3	/* pam2[0,1] matrix */
#define STARTTYPE4	4	/* *pascii for fasty/tfasty */

#define QSEQTYPE0	5
#define QSEQTYPE1	6

#define MSEQTYPE0	10	/* cur_buf->hdr */
#define MSEQTYPE1	11	/* cur_buf->buf2_data * cur_buf->hdr.buf2_cnt */
#define MSEQTYPE2	12	/* bulk - seq_b, seq_record * hdr.buf2_cnt */
#define MSEQTYPE3	13	/* bulk - aa1b_start, aa1b_used+1 */
#define MSEQTYPE4	14	/* individ. - seq_record */
#define MSEQTYPE5	15	/* individ. - aa1b */
#define MSEQTYPE6	16

#define RES_TYPE0	20
#define RES_TYPE1	21
#define RES_TYPE2	22

#define ALN_TYPE0	30
#define ALN_TYPE1	31
#define ALN_TYPE2	32
#define ALN_TYPE3	33

#define FINISHED 	16384	/* this must be larger than BFR */

#define DO_SEARCH_FLG 0
#define DO_OPT_FLG 1
#define DO_ALIGN_FLG 2
#define DO_CALC_FLG 3
