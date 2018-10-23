/* map_db.c - read a FASTA or GCG format database and generate a list
   of indices for rapid memory mapping */

/* $Id: map_db.c 1239 2013-11-02 01:09:58Z wrp $ */

/* copyright (c) 1999, 2014 by William R. Pearson and The Rector &
   Visitors of the University of Virginia */

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

/* input is a lib_type 1,5, or 6 sequence database (lib_type specified after filename), 
   e.g. 'swissprot.lseg 1' */
/* map_db -n specifies a DNA database */

/* output is a BLAST2 formatdb type index file */

/* format of the index file:

1)  map_db version number ["MP"+2 bytes]
2)  number of sequences in database [4 bytes]
3)  total length of database        [8 bytes]  (MP1, 4 bytes for MP0)
4)  longest sequence in database    [8 bytes]  (MP1, 4 bytes for MP0)
5) list of offsets to definitions  [num_seq+1] int*8 (MP1, 4 bytes for MP0)
6) list of offsets to sequences    [num_seq+1] int*8 (MP1, 4 bytes for MP1)
7) list of flag characters for sequences [num_seq+1]bytes
    (used for GCG binary to encode 2bit or 4 bit representation)

    sequence files will be as defined by their format
*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <sys/types.h>
#include <sys/stat.h>

#ifndef FSEEK_T_DEF
#ifndef USE_FSEEKO
#define FSEEK_T_DEF
#define FSEEK fseek
#define FTELL ftell
typedef long fseek_t;
#else
#define FSEEK fseeko
#define FTELL ftello
typedef off_t fseek_t;
#endif
#endif

#define LASTLIB 6

int (*get_entry) ();

int a_get_ent(unsigned char *, int, fseek_t *, fseek_t *);
int gbf_get_ent(unsigned char *, int, fseek_t *, fseek_t *);

void src_int4_write(FILE *, int);
void src_int4_read(FILE *, int *);
void src_long4_write(FILE *, long);
void src_long4_read(FILE *, long *);
void src_long8_write(FILE *, int64_t);
void src_long8_read(FILE *, int64_t *);

void newname(char *nname, char *oname, char *suff, int maxn);
void init_ascii0(int *xascii, char *sq_map);

int (*get_ent_arr[LASTLIB+1])()={a_get_ent, gbf_get_ent, NULL, NULL, NULL,
				 NULL, NULL};

fseek_t openlib(char *, int);

#define NA 123
#define TERM 24
#define EL 125
#define ES 126
#define AAMASK 127

static int *sascii, aascii[128];
char *NCBIstdaa_ext = "-ABCDEFGHIKLMNPQRSTVWXYZU*OJ-abcdefghiklmnpqrstvwxyzu*oj";
char NCBIstdaa_ext_n = 56;

int
main(int argc, char **argv) {
  FILE *libi, *b_fd;  /* b_fd is used for both the binary data file and the index file */
  char lname[256];
  char bname[256];
  char iname[256];
  char format[4];
  char *bp;
  struct stat stat_buf;

  int i;
  int nlib;	/* number of entries */
  long max_len;	/* longest sequence */
  fseek_t tot_len;	/* total sequence length */

  int n1;
  
  fseek_t f_size;	/* file size from fstat() */
  int lib_size;	/* current space available - may be realloc'ed */
  int lib_inc;
  int lib_type; /* 1 for protein, 0 for DNA */
  int lib_aa;	/* dna=1; prot=0; */
  int build_binary; /* build binary sequence file */

  /* file offsets */
  fseek_t d_pos;	/* start of description */
  fseek_t s_pos;	/* start of sequence */
  fseek_t b_pos;	/* start of binary encoding */
  fseek_t *d_pos_arr;	/* array of description pointers */
  fseek_t *s_pos_arr;	/* array of ascii sequence pointers */
  fseek_t *b_pos_arr;	/* array of binary sequence pointers */
  unsigned char *zbuff;	/* tmp buffer for writes */
  unsigned char *sbuff;	/* sequence buffer */
  int sbuff_max, sbuff_dup;

  lib_type = 0;
  lib_size = 200000;
  lib_inc  = 100000;
  lib_aa = 1;
  sbuff = NULL;
  sbuff_max = 50000;
  sbuff_dup = 0;
  build_binary = 0;

  while (argc > 1 && *argv[1]=='-') {
    if (strcmp(argv[1],"-n")==0) lib_aa = 0;
    else if (strcmp(argv[1], "-b")==0) build_binary = 1;
    argv++;
    argc--;
  }

  /* open the database */
  if (argc > 1) strncpy(lname, argv[1],sizeof(lname));
  else {
    fprintf(stderr," Entry library name: ");
    fgets(lname,sizeof(lname),stdin);
    if ((bp=strchr(lname,'\n'))!=NULL) *bp='\0';
  }
    
  if ((bp=strchr(lname,' '))!=NULL) {
    lib_type = atoi(bp+1);
    *bp='\0';
  }
  else lib_type = 0;

  if (get_ent_arr[lib_type] == NULL) {
    fprintf(stderr," cannot index file %s type %d\n",lname,lib_type);
    exit(1);
  }
  
  if (lib_type == 6) lib_aa = 0;
  if (lib_type == 1) lib_aa = 0;
  
  if (lib_aa == 1) {
    init_ascii0(aascii, NCBIstdaa_ext);
  }
  else {
    if (build_binary) {
      fprintf(stderr,"*** WARNING ***  map_db -- binary files not available for DNA libraries\n");
      build_binary = 0;
    }
    init_ascii0(aascii, "\0ACGTURYMWSKDHVBNacgturymwskdhvbn");
    aascii['X'] = aascii['N'];
    aascii['x'] = aascii['n'];
  }
  sascii = &aascii[0];

  if ((f_size=openlib(lname,lib_type))==0) {
    fprintf(stderr," cannot open %s (type: %d)\n",lname,lib_type);
    exit(1);
  }

  if (build_binary) {
    newname(bname, lname, "bsq",sizeof(bname));
    if ((b_fd = fopen(bname, "w")) == NULL) {
      fprintf(stderr, "cannot open %s binary file for writing\n",bname);
      exit(1);
    }

    if ((zbuff=(unsigned char *)calloc(256,sizeof(char)))==NULL) {
      fprintf(stderr," cannot allocate zbuff[256]\n");
      exit(1);
    }

    if ((sbuff=(unsigned char *)calloc(sbuff_max,sizeof(char)))==NULL) {
      fprintf(stderr," cannot allocate sbuff[%d]\n",sbuff_max);
      exit(1);
    }

    /* write out the initial NULL */
    fwrite(zbuff,sizeof(char),1, b_fd);
  }
  b_pos = 1;	/* initialize whether used or not */

  /* allocate array of description pointers */
  if ((d_pos_arr=(fseek_t *)calloc(lib_size, sizeof(fseek_t)))==NULL) {
    fprintf(stderr," cannot allocate %d for desc. array\n",lib_size);
    exit(1);
  }
  /* allocate array of sequence pointers */
  if ((s_pos_arr=(fseek_t *)calloc(lib_size, sizeof(fseek_t)))==NULL) {
    fprintf(stderr," cannot allocate %d for seq. array\n",lib_size);
    exit(1);
  }
  /* allocate array of sequence pointers */
  if ((b_pos_arr=(fseek_t *)calloc(lib_size, sizeof(fseek_t)))==NULL) {
    fprintf(stderr," cannot allocate %d for binary array\n",lib_size);
    exit(1);
  }

  /* allocate array of sequence flags */

  nlib = 0; tot_len=0; max_len=-1;
  while ((n1=get_entry(sbuff, sbuff_max, &d_pos, &s_pos)) > 0) {
    if (build_binary) fwrite(sbuff, sizeof(char), n1+1, b_fd);

    d_pos_arr[nlib] = d_pos;
    s_pos_arr[nlib] = s_pos;
    b_pos_arr[nlib] = b_pos;

    b_pos += n1+1;

    nlib++;
    tot_len += n1;

    if (n1 > max_len) max_len = n1;
    if (nlib >= lib_size) { /* too many entries */

      lib_size += lib_inc;
      if ((d_pos_arr=(fseek_t *)realloc(d_pos_arr,lib_size*sizeof(fseek_t)))==NULL) {
	fprintf(stderr," cannot realloc allocate %d for desc.. array\n",
		lib_size);
	exit(1);
      }
      if ((s_pos_arr=(fseek_t *)realloc(s_pos_arr,lib_size*sizeof(fseek_t)))==NULL) {
	fprintf(stderr," cannot realloc allocate %d for seq. array\n",
		lib_size);
	exit(1);
      }
      if ((b_pos_arr=(fseek_t *)realloc(b_pos_arr,lib_size*sizeof(fseek_t)))==NULL) {
	fprintf(stderr," cannot realloc allocate %d for binary aseq. array\n",
		lib_size);
	exit(1);
      }
    }
  }

  if (build_binary) fclose(b_fd);

  if (stat(lname,&stat_buf)<0) {
    fprintf(stderr," cannot stat library: %s\n",lname);
    exit(1);
  }
  else {
    f_size = stat_buf.st_size;
  }

  d_pos_arr[nlib]= d_pos;	/* put in the end of the file */
  s_pos_arr[nlib]=0;
  b_pos_arr[nlib]= b_pos;

  /* all the information is in, write it out */
  
  newname(iname,lname,"xin",sizeof(iname));
  if ((libi=fopen(iname,"w"))==NULL) {
    fprintf(stderr," cannot open %s for writing\n",iname);
    exit(1);
  }

  /* write out format version, etc. to .xin file */
  format[0]='M';
  format[1]='P';
#ifdef BIG_LIB64
  format[2]= 1;		/* format 1,2 for 8-byte offsets */
#else
  format[2]='\0';	/* format '\0' for original 4-byte */
#endif

  format[3]=lib_type;
  fwrite(format,4,sizeof(char),libi);

  /* write out sequence type */
  src_int4_write(libi, lib_aa);
  /* write out file fstat as integrity check */
#ifdef BIG_LIB64
  src_long8_write(libi, f_size);
#else
  src_int4_write(libi, f_size);
#endif
  /* write out num_seq */
  src_int4_write(libi, nlib);

#ifdef BIG_LIB64
  /* write out tot_len, max_len */
  src_long8_write(libi, tot_len);
#else
  src_int4_write(libi, tot_len);
#endif
  src_int4_write(libi, max_len);

#ifdef BIG_LIB64
  for (i=0; i<=nlib; i++) src_long8_write(libi,d_pos_arr[i]);
  for (i=0; i<=nlib; i++) src_long8_write(libi,s_pos_arr[i]);
#else
  for (i=0; i<=nlib; i++) src_int4_write(libi,d_pos_arr[i]);
  for (i=0; i<=nlib; i++) src_int4_write(libi,s_pos_arr[i]);
#endif
  fclose(libi);

  if (build_binary) { /* do the same thing for the .xin_b file */
    if (stat(bname,&stat_buf)<0) {
      fprintf(stderr," cannot stat library: %s\n",bname);
      exit(1);
    }
    else {
      f_size = stat_buf.st_size;
    }

    newname(iname,lname,"xin_b",sizeof(iname));
    if ((b_fd=fopen(iname,"w"))==NULL) {
      fprintf(stderr," cannot open %s for writing\n",iname);
      exit(1);
    }

    /* write out format version, etc. to .xin file */
    format[0]='M';
    format[1]='P';
#ifdef BIG_LIB64
    format[2]= 2;		/* format 1,2 for 8-byte offsets */
#else
    format[2]='\0';	/* format '\0' for original 4-byte */
#endif

    format[3]=lib_type;
    fwrite(format,4,sizeof(char),libi);

    /* write out sequence type */
    src_int4_write(libi, lib_aa);
    /* write out file fstat as integrity check */
    src_long8_write(libi, f_size);
    /* write out num_seq */
    src_int4_write(libi, nlib);

    /* write out tot_len, max_len */
    src_long8_write(libi, tot_len);
    src_int4_write(libi, max_len);
    /* write out maximum length */
    src_int4_write(libi, sbuff_max);
    /* write out overlap */
    src_int4_write(libi, sbuff_dup);

    for (i=0; i<=nlib; i++) src_long8_write(libi,b_pos_arr[i]);
    fclose(libi);
  }

#ifdef BIG_LIB64
  fprintf(stderr," wrote %d sequences (tot=%lld, max=%ld) to %s\n",
	  nlib,tot_len,max_len,iname);
#else
  fprintf(stderr," wrote %d sequences (tot=%ld, max=%ld) to %s\n",
	  nlib,tot_len,max_len,iname);
#endif
  if (build_binary) {
    fprintf(stderr," binary file: %s.bsq written\n",lname);
  }

  exit(0);
}


FILE *libf=NULL;
fseek_t lpos;

#define MAXLINE 4096
char lline[MAXLINE+1];

fseek_t
openlib(char *lname, int lib_type)
{
  struct stat stat_buf;

  if (stat(lname,&stat_buf)<0) {
    fprintf(stderr," cannot stat library: %s\n",lname);
    return 0;
  }

  if ((libf=fopen(lname,"r"))==NULL) {
    fprintf(stderr," cannot open library: %s (type: %d)\n",
	    lname, lib_type);
    return 0;
  }
  
  get_entry = get_ent_arr[lib_type];

  lpos = FTELL(libf);
  if (fgets(lline,MAXLINE,libf)==NULL) return 0;
  return stat_buf.st_size;
}

int
a_get_ent(unsigned char *sbuff, int max_sbuff, 
	  fseek_t *d_pos, fseek_t *s_pos)
{
  char *cp;
  unsigned char *sptr, *sptr_max;
  int *ap, n1;

  sptr = sbuff;
  sptr_max = sbuff+max_sbuff;

  ap = sascii;

  while (lline[0]!='>' && lline[0]!=';') {
    lpos = FTELL(libf);
    if (fgets(lline,sizeof(lline),libf)==NULL) {
      *d_pos = lpos;
      return 0;
    }
  }

  *d_pos = lpos;

  /* make certain we have the end of the line */
  while (strchr((char *)lline,'\n')==NULL) {
    if (fgets(lline,sizeof(lline),libf)==NULL) break;
  }

  *s_pos = FTELL(libf);
  lline[0]='\0';
  n1 = 0;
  while (fgets(lline,sizeof(lline),libf)!=NULL) {
    if (lline[0]=='>') break;
    if (lline[0]==';') {
      if (strchr(lline,'\n')==NULL) {
	fprintf(stderr," excessive continuation\n%s",lline);
	return -1;
      }
    }

    if (sbuff) {
      for (cp=lline; *cp && sptr < sptr_max;) {
	if ((*sptr = ap[*cp++])<NA) {sptr++;}
      }

      if (sptr >= sptr_max) {
	fprintf(stderr," sequence too long: %ld\n",(long)(sptr-sbuff));
	exit(1);
      }
    }
    else {
      for (cp=lline; *cp; ) if (ap[*cp++]<NA) n1++;
    }

    lpos = FTELL(libf);
  }
  if (sbuff) {
    *sptr = '\0';
    return sptr - sbuff;
  }
  else {
    return n1;
  }
}

int
gbf_get_ent(unsigned char *sbuff, int max_sbuff, 
	    fseek_t *d_pos, fseek_t *s_pos)
{
  int n1;
  char *cp;
  register int *ap;

#if !defined(TFAST)
  ap = sascii;
#else
  ap = nascii;
#endif

  while (lline[0]!='L' || lline[1]!='O' || 
	 strncmp(lline,"LOCUS",5)) { /* find LOCUS */
    lpos = FTELL(libf);
    if (fgets(lline,MAXLINE,libf)==NULL) return (-1);
  }
  *d_pos=lpos;

  while (lline[0]!='O' || lline[1]!='R' ||
	 strncmp(lline,"ORIGIN",6)) { /* find ORIGIN */
    if (fgets(lline,MAXLINE,libf)==NULL) return (-1);
  }
  *s_pos = FTELL(libf);

  lline[0]='\0';
  n1=0;
  while (fgets(lline,MAXLINE,libf)!=NULL) {
    if (lline[0]=='/') break;
    for (cp=lline; *cp; ) if (ap[*cp++]<NA) n1++;
  }
  lpos = FTELL(libf);
  fgets(lline,MAXLINE,libf);

  return n1;
}

void src_int4_read(FILE *fd,  int *val)
{
#ifdef IS_BIG_ENDIAN
  fread((char *)val,(size_t)4,(size_t)1,fd);
#else
  unsigned char b[4];

  fread((char *)&b[0],(size_t)1,(size_t)4,fd);
  *val = 0;
  *val = (int)((int)((int)(b[0]<<8)+((int)b[1]<<8))+((int)b[2]<<8))
	  +(int)b[3];
#endif
}

void src_int4_write(FILE *fd,  int val)
{
#ifdef IS_BIG_ENDIAN
  fwrite(&val,(size_t)4,(size_t)1,fd);
#else
  unsigned char b[4];

  b[3] = val & 255;
  b[2] = (val=val>>8)&255;
  b[1] = (val=val>>8)&255;
  b[0] = (val=val>>8)&255;

  fwrite(b,(size_t)1,(size_t)4,fd);
#endif
}

void src_long8_write(FILE *fd,  int64_t val)
{
#ifdef IS_BIG_ENDIAN
  fwrite(&val,(size_t)8,(size_t)1,fd);
#else
  unsigned char b[8];

  b[7] = val & 255;
  b[6] = (val=val>>8)&255;
  b[5] = (val=val>>8)&255;
  b[4] = (val=val>>8)&255;
  b[3] = (val=val>>8)&255;
  b[2] = (val=val>>8)&255;
  b[1] = (val=val>>8)&255;
  b[0] = (val=val>>8)&255;

  fwrite(b,(size_t)1,(size_t)8,fd);
#endif
}

void
newname(char *nname, char *oname, char *suff, int maxn)
{
  strncpy(nname,oname,maxn-1);
  strncat(nname,".",1);
  strncat(nname,suff,maxn-strlen(nname));
}

/* init_ascii0 -- initializes an ascii mapping from a sequence
   ordering
*/
void 
init_ascii0(int *xascii, char *sq_map) {
  int i;
  int n_sq_map;

  n_sq_map = strlen(sq_map+1) + 1;

  /* first map everything as non-sequence */
  for (i=0; i<128; i++) {
    xascii[i] = NA;
  }

  /* then map the actual sequence letters */
  for (i = 1; i < n_sq_map; i++) {
    xascii[sq_map[i]] = i;
    if (n_sq_map <= 28) { /* only uppercase */
      xascii[sq_map[i]+32] = i;	/* map lowercase */
    }
  }

  /* then map the other stuff, EL  etc */
  xascii[0] = ES;
  xascii[10] = EL;
  xascii[13] = EL;
}

