/* list_db.c - report values from map_db.c  */

/* $Id: list_db.c $ */

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

/* format of the index file:

1)  map_db version number ["MP"+2 bytes]
2)  number of sequences in database [4 bytes]
3)  total length of database        [8 bytes]
4)  longest sequence in database    [8 bytes]
5) list of offsets to definitions  [num_seq+1] int*8
6) list of offsets to sequences    [num_seq+1] int*8
7) list of flag characters for sequences [num_seq+1] bytes
    (used for GCG binary to encode 2bit or 4 bit representation)

    sequence files will be as defined by their format
*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "uascii.h"
#include "mm_file.h"
/* #include "ncbl2_head.h" */

void src_int4_write(FILE *, int);
void src_int4_read(FILE *, int *);
void src_long4_read(FILE *, long *);
void src_long8_write(FILE *, int64_t);
void src_long8_read(FILE *, int64_t *);

void newname(char *nname, char *oname, char *suff, int maxn);

int
main(int argc, char **argv)
{
  FILE *libi;
  char lname[256];
  char iname[256];
  char format[4];
  char *bp;

  int i;
  int d_pos;	/* start of description */
  int s_pos;	/* start of sequence */
  int attr;	/* sequence attribute */
  int lib_aa;	/* 0 => DNA, 1 => protein */
  int nlib;	/* number of entries */
  fseek_t f_size;
  long lf_size;	/* (long) version of f_size */
  long max_len;	/* longest sequence */
  MM_OFF tot_len;	/* total sequence length */
  int n1;
  
  int lib_size;	/* current space available - may be realloc'ed */
  int lib_inc;
  int lib_type; /* 1 for protein, 0 for DNA */
  int lib_dna;	/* dna=1; prot=0; */
  fseek_t *d_pos_arr;	/* array of description pointers */
  fseek_t *s_pos_arr;	/* array of description pointers */
  long *tmp_pos_arr;
  char *attr_arr;	/* array of attribute chars */

  int mm64_flag;

  lib_type = 0;
  lib_dna = 0;

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

  newname(iname,lname,"xin",sizeof(iname));

  if ((libi=fopen(iname,"r"))==NULL) {
    fprintf(stderr," cannot open %s\n",iname);
    exit(1);
  }

  fread(format,1,sizeof(format),libi);
  printf("%c%c%d %d\n",format[0],format[1],format[2],format[3]);
  mm64_flag = (format[2]==1);

  src_int4_read(libi,&lib_aa);

  if (mm64_flag) src_long8_read(libi,&f_size);
  else {
    src_long4_read(libi,&lf_size);
    f_size = lf_size;
  }

  src_int4_read(libi,&nlib);

  if (mm64_flag) {
    src_long8_read(libi,&tot_len);
  }
  else {
    src_long4_read(libi,&lf_size);
    tot_len = lf_size;
  }
  src_long4_read(libi,&max_len);

  printf(" %d entries; tot: %lld; max: %ld\n",nlib,tot_len,max_len);

  /* allocate array of description pointers */
  if (!mm64_flag) {
    if ((tmp_pos_arr=(long *)calloc(nlib+1,sizeof(long)))==NULL) {
      fprintf(stderr," cannot allocate %d for tmp_pos array\n",
	      nlib+1);
    }
  }

  /* allocate array of description pointers */
  if ((d_pos_arr=(fseek_t *)calloc(nlib+1, sizeof(fseek_t)))==NULL) {
    fprintf(stderr," cannot allocate %d for desc. array\n",nlib+1);
    exit(1);
  }
  /* allocate array of sequence pointers */
  if ((s_pos_arr=(fseek_t *)calloc(nlib+1, sizeof(fseek_t)))==NULL) {
    fprintf(stderr," cannot allocate %d for seq. array\n",nlib+1);
    exit(1);
  }
  if ((attr_arr=(char *)calloc(nlib+1, sizeof(char)))==NULL) {
    fprintf(stderr," cannot allocate %d for attr. array\n",nlib+1);
    exit(1);
  }
  
  if (mm64_flag) {
    for (i=0; i<=nlib; i++) src_long8_read(libi,&d_pos_arr[i]);
    for (i=0; i<=nlib; i++) src_long8_read(libi,&s_pos_arr[i]);
  }
  else {
    for (i=0; i<=nlib; i++) src_long4_read(libi,&tmp_pos_arr[i]);
    for (i=0; i<=nlib; i++) d_pos_arr[i] = tmp_pos_arr[i];
    for (i=0; i<=nlib; i++) src_long4_read(libi,&tmp_pos_arr[i]);
    for (i=0; i<=nlib; i++) s_pos_arr[i] = tmp_pos_arr[i];
  }

  fread(attr_arr,nlib+1,sizeof(char),libi);
  fclose(libi);

  printf("header\tseq\n");

  for (i=0; i<nlib; i++) printf("%lld\t%lld\n",d_pos_arr[i],s_pos_arr[i]);

  exit(0);
}

void src_int4_read(FILE *fd,  int *val)
{
  int tval;
#ifdef IS_BIG_ENDIAN
  fread(&tval,(size_t)4,(size_t)1,fd);
  *val = tval;
#else
  unsigned char b[4];

  fread((char *)&b[0],(size_t)1,(size_t)4,fd);
  *val = 0;
  *val = (int)((int)((int)(b[0]<<8)+(int)b[1]<<8)+(int)b[2]<<8)
	  +(int)b[3];
#endif
}

void src_long4_read(FILE *fd,  long *val)
{
  int tval;
#ifdef IS_BIG_ENDIAN
  fread(&tval,(size_t)4,(size_t)1,fd);
  *val = tval;
#else
  unsigned char b[4];

  fread((char *)&b[0],(size_t)1,(size_t)4,fd);
  *val = 0;
  *val = (int)((int)((int)(b[0]<<8)+(int)b[1]<<8)+(int)b[2]<<8)
	  +(int)b[3];
#endif
}

void src_long8_read(FILE *fd,  int64_t *val)
{
#ifdef IS_BIG_ENDIAN
  fread((char *)val,(size_t)8,(size_t)1,fd);
#else
  unsigned char b[8];

  fread((char *)&b[0],(size_t)1,(size_t)8,fd);
  *val = 0;
  *val = (int)
    ((((((((int)b[0]<<8)+(int)b[1]<<8)+(int)b[2]<<8)+(int)b[3]<<8)+
		(int)b[4]<<8)+(int)b[5]<<8)+(int)b[6]<<8)+(int)b[7];
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

void
newname(char *nname, char *oname, char *suff, int maxn)
{
  strncpy(nname,oname,maxn-1);
  strncat(nname,".",1);
  strncat(nname,suff,maxn-strlen(nname));
}
