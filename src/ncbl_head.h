/* ncbl_head.h	header files for blast1.3 format */

/* $Id: ncbl_head.h 625 2011-03-23 17:21:38Z wrp $ */
/* $Revision: 625 $  */

#define AMINO_ACID_SEQTYPE	1
#define AA_SEQTYPE	AMINO_ACID_SEQTYPE
#define NUCLEIC_ACID_SEQTYPE	2
#define NT_SEQTYPE	NUCLEIC_ACID_SEQTYPE

/* Filename extensions used by the two types of databases (a.a. and nt.) */
#define AA_HEADER_EXT	"ahd"
#define AA_TABLE_EXT	"atb"
#define AA_SEARCHSEQ_EXT	"bsq"
#define NT_HEADER_EXT	"nhd"
#define NT_TABLE_EXT	"ntb"
#define NT_SEARCHSEQ_EXT	"csq"

#define DB_TYPE_PRO	0x78857a4f	/* Magic # for a protein sequence database */
#define DB_TYPE_NUC	0x788325f8	/* Magic # for a nt. sequence database */

#define AAFORMAT	3	/* Latest a.a. database format ID number */
#define NTFORMAT	6	/* Latest nt. database format ID number */

#define NULLB		'\0'	/* sentinel byte */
#define NT_MAGIC_BYTE	0xfc	/* Magic byte at end of compressed nt db */

#ifndef CHAR_BIT
#define CHAR_BIT	8	/* these values should match blast */
#endif

#define NBPN		2
#define NSENTINELS	2
