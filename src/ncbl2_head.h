/* ncbl_head.h	header files for blast1.3 format */

/* $Id: ncbl2_head.h 625 2011-03-23 17:21:38Z wrp $ */
/* $Revision: 625 $  */

#define AMINO_ACID_SEQTYPE	1
#define AA_SEQTYPE	AMINO_ACID_SEQTYPE
#define AAFORMAT	AA_SEQTYPE

#define NUCLEIC_ACID_SEQTYPE	0
#define NT_SEQTYPE	NUCLEIC_ACID_SEQTYPE
#define NTFORMAT	NT_SEQTYPE

/* Filename extensions used by the two types of databases (a.a. and nt.) */
#define AA_LIST_EXT	"pal"
#define AA_HEADER_EXT	"phr"
#define AA_INDEX_EXT	"pin"
#define AA_SEARCHSEQ_EXT	"psq"

#define NT_LIST_EXT	"nal"
#define NT_HEADER_EXT	"nhr"
#define NT_INDEX_EXT	"nin"
#define NT_SEARCHSEQ_EXT	"nsq"

#define FORMATDBV3	3	/* formatdb version */
#define FORMATDBV4	4	/* formatdb version */
#define FORMATDBV5	5	/* formatdb version */

#define NULLB		'\0'	/* sentinel byte */

#ifndef CHAR_BIT
#define CHAR_BIT	8	/* these values should match blast */
#endif

#define NBPN		2
#define NSENTINELS	2
