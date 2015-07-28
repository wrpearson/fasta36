/* $Id: mysql_lib.c 817 2011-08-02 03:54:02Z wrp $ */
/* $Revision: 817 $  */

/* mysql_lib.c copyright (c) 2000, 2014 by William R. Pearson and The
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

/* functions for opening, reading, seeking a mySQL database */

/* (close_tables added June, 2011)
  For the moment, this interface assumes that the file to be searched
  will be specified in a single, long, string with 4 required and 1
  optional parts:

  (1) a database open string. This string has four fields, separated by
      whitespace (' \t'):
        hostname:port dbname user password

   '--' dashes at the beginning of lines are ignored -
   thus the first line could be:
   -- hostname:port dbname user password

  (2) a database query string that will return an unique ID (not
      necessarily numberic, but it must be < 12 characters as libstr[12]
      is used) and a sequence string

  (2a) a series of mySQL commands that do not generate results
       starting with 'DO', followed by a select() statement.

  (3) a database select string that will return a description
      given a unique ID

  (4) a database select string that well return a sequence given a
      unique ID

  (5) [optional] an SQL statement to be run when closing the database
      (e.g. a DROP TABLE statement)

   Lines (3) and (4) are not required for pv34comp* libraries, but
   line (2) must generate a complete description as well as a sequence.

   18-July-2001
   Additional syntax has been added to support multiline SQL queries.

   If the host line begins with '+', then the SQL is openned on the same
   connection as the previous SQL file.

   If the host line contains '-' just before the terminal ';', then
   the file will not produce any output.

   This string can contain "\n". ";" are used to separate the four
   functions, which must be specified in the order shown above.
   The last (fourth) query must terminate with a ';' */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>

#include <mysql.h>
#define MYSQL_LIB 16

#include "defs.h"
#include "structs.h"
#include "mm_file.h"

#define XTERNAL
#include "uascii.h"
#define EOSEQ 0
/* #include "upam.h" */

char *alloc_file_name(char *f_name);
int mysql_getlib(unsigned char *, int, char *, int, fseek_t *, int *, struct lmf_str *, long *);
void mysql_ranlib(char *, int, fseek_t, char *, struct lmf_str *m_fd);

#define MYSQL_BUF 4096

struct lmf_str *
mysql_openlib(char *sname, int ldnaseq, int *sascii) {
  FILE *sql_file;
  char *tmp_str, *ttmp_str;
  int tmp_str_len;
  char *bp, *bps, *bdp, *tp, tchar;
  int i, qs_len, qqs_len;
  char *sql_db, *sql_host, *sql_dbname, *sql_user, *sql_pass;
  char *sql_do;
  int sql_do_cnt;
  int sql_port;
  struct lmf_str *m_fptr;

  /*  if (sql_reopen) return NULL; - should not be called for re-open */

  tmp_str_len = MYSQL_BUF;
  if ((tmp_str=(char *)calloc(tmp_str_len,sizeof(char)))==NULL) {
    fprintf(stderr,"cannot allocate %d for mySQL buffer\n",tmp_str_len);
    return NULL;
  }

  /* immediate mysql scripts start with '%' */
  if (sname[0] == '%') {
    strncpy(tmp_str,sname+1,tmp_str_len);
    tmp_str[sizeof(tmp_str)-1]='\0';
  }
  else {	/* read the script from a file */
    if ((sql_file=fopen(sname,"r"))==NULL) {
      fprintf(stderr," cannot open mySQL file: %s\n",sname);
      return NULL;
    }

    if ((qs_len=fread(tmp_str,sizeof(char),tmp_str_len-1,sql_file))<=0) {
      fprintf(stderr," cannot read mySQL file: %s\n",sname);
      return NULL;
    }
    else  {	/* read the entire file in MYSQL_BUF (4096) byte
		   chunks, reallocating as necessary */
      tmp_str[qs_len]='\0';
      qqs_len = qs_len;
      while (qqs_len >= tmp_str_len-1) {
	tmp_str_len += MYSQL_BUF;
	if ((tmp_str=(char *)realloc(tmp_str,tmp_str_len))==NULL) {
	  fprintf(stderr,
		  " cannot reallocate %d for mySQL buffer\n",tmp_str_len);
	  return NULL;
	}
	ttmp_str = &tmp_str[qqs_len];
	if ((qs_len=fread(ttmp_str,sizeof(char),MYSQL_BUF,sql_file))<0) {
	  fprintf(stderr," cannot read mySQL file: %s\n",sname);
	  return NULL;
	}
	ttmp_str[qs_len]='\0';
	qqs_len += qs_len;
      }
    }
    fclose(sql_file);
  }

  /* tmp_str has the entire contents of the file */
  bps = tmp_str;
  if ((bp=strchr(bps,';'))!=NULL) {
    /* get the connection info */
    *bp='\0';
    if ((sql_db=calloc(strlen(bps)+1,sizeof(char)))==NULL) {
      fprintf(stderr, " cannot allocate space for database name [%d], %s\n",
	      (int)strlen(bps),bps);
      return NULL;
    }
    /* have database name, parse the fields */
    else {  /* copy connection info into sql_db and parse */
      strcpy(sql_db,bps);	/* strcpy OK because allocated strlen(bps) */
      bps = bp+1;	/* points to next char after ';' */
      while (isspace(*bps)) bps++;
      *bp=';'; /* replace ; */
      bp = sql_db;
      while (*bp=='-') {*bp++ = ' ';}
      sql_host = strtok(bp," \t\n");
      sql_dbname = strtok(NULL," \t\n");
      sql_user = strtok(NULL," \t\n");
      sql_pass = strtok(NULL," \t\n");
      if ((tp=strchr(sql_host,':'))!=NULL) {
	*tp='\0';
	sql_port=atoi(tp+1);
      }
      else sql_port = 0;
    }
  }
  else {
    fprintf(stderr," cannot find database fields:\n%s\n",tmp_str);
    return NULL;
  }

  /* we have all the info we need to open a database, allocate lmf_str */
  if ((m_fptr = (struct lmf_str *)calloc(1,sizeof(struct lmf_str)))==NULL) {
    fprintf(stderr," cannot allocate lmf_str (%ld) for %s\n",
	    sizeof(struct lmf_str),sname);
    return NULL;
  }

  /* have our struct, initialize it */

  m_fptr->lb_name = alloc_file_name(sname);

  m_fptr->sascii = sascii;

  m_fptr->sql_db = sql_db;
  m_fptr->getlib = mysql_getlib;
  m_fptr->ranlib = mysql_ranlib;
  m_fptr->mm_flg = 0;
  m_fptr->sql_reopen = 0;
  m_fptr->lb_type = MYSQL_LIB;

  /* now open the database, if necessary */
  if ((m_fptr->mysql_conn=mysql_init(NULL))==NULL) {
    fprintf(stderr,"*** Error - mysql_init\n");
    goto error_r;
  }

  if (mysql_real_connect(m_fptr->mysql_conn,
			 sql_host,sql_user,sql_pass,
			 sql_dbname,
			 sql_port,
			 NULL,
			 0)==NULL)
    {
      fprintf(stderr,"*** Error %u - could  not open database:\n%s\n%s",
	      mysql_errno(m_fptr->mysql_conn),tmp_str,
	      mysql_error(m_fptr->mysql_conn));
      goto error_r;
    }
#ifdef DEBUG
  else {
    fprintf(stderr," Database %s opened on %s\n",sql_dbname,sql_host);
  }
#endif

  /* check for 'DO' command - copy to 'DO' string */
  while (*bps == '-') { *bps++=' ';}
  if (isspace(bps[-1]) && toupper(bps[0])=='D' &&
      toupper(bps[1])=='O' && isspace(bps[2])) {
    /* have some 'DO' commands */
    /* check where the end of the last DO statement is */

    sql_do_cnt = 1;	/* count up the number of 'DO' statements for later */
    bdp=bps+3;
    while ((bp=strchr(bdp,';'))!=NULL) {
      tp = bp+2; /* skip ;\n */
      while (isspace(*tp) || *tp == '-') {*tp++ = ' ';}
      if (toupper(*tp)=='D' && toupper(tp[1])=='O' && isspace(tp[2])) {
	sql_do_cnt++;		/* count the DO statements */
	bdp = tp+3;		/* move to the next DO statement */
      }
      else break;
    }
    if (bp != NULL) {	/* end of the last DO, begin of select */
      tchar = *(bp+1);
      *(bp+1)='\0';		/* terminate DO strings */
      if ((sql_do = calloc(strlen(bps)+1, sizeof(char)))==NULL) {
	fprintf(stderr," cannot allocate %d for sql_do\n",(int)strlen(bps));
	goto error_r;
      }
      else {
	strcpy(sql_do,bps);
	*(bp+1)=tchar;	/* replace missing ';' */
      }
      bps = bp+1;
      while (isspace(*bps)) bps++;
    }
    else {
      fprintf(stderr," terminal ';' not found: %s\n",bps);
      goto error_r;
    }
    /* all the DO commands are in m_fptr->sql_do in the form: 
     DO command1; DO command2; DO command3; */
    bdp = sql_do;
    while (sql_do_cnt-- && (bp=strchr(bdp,';'))!=NULL) {
      /* do the mysql statement on bdp+3 */
      /* check for error */
      *bp='\0';
      if (mysql_query(m_fptr->mysql_conn,bdp+3)) {
	fprintf(stderr,"*** Error %u - query failed:\n%s\n%s\n",
		mysql_errno(m_fptr->mysql_conn), bdp+3, mysql_error(m_fptr->mysql_conn));
	goto error_r;
      }
      *bp=';';
      bdp = bp+1;
      while (isspace(*bdp)) bdp++;
    }
  }

  /* copy 1st query field */
  if ((bp=strchr(bps,';'))!=NULL) {
    *bp='\0';
    if ((m_fptr->sql_query=calloc(strlen(bps)+1,sizeof(char)))==NULL) {
      fprintf(stderr, " cannot allocate space for query string [%d], %s\n",
	      (int)strlen(bps),bps);
      goto error_r;
    }
    /* have query, copy it */
    else {
      strcpy(m_fptr->sql_query,bps);
      *bp=';'; /* replace ; */
      bps = bp+1;
      while(isspace(*bps)) bps++;
    }
  }
  else {
    fprintf(stderr," cannot find database query field:\n%s\n",tmp_str);
    goto error_r;
  }

  /* copy get_desc field */
  if ((bp=strchr(bps,';'))!=NULL) {
    *bp='\0';
    if ((m_fptr->sql_getdesc=calloc(strlen(bps)+1,sizeof(char)))==NULL) {
      fprintf(stderr, " cannot allocate space for database name [%d], %s\n",
	      (int)strlen(bps),bps);
      goto error_r;
    }
    /* have get_desc, copy it */
    else {
      strcpy(m_fptr->sql_getdesc,bps);
      *bp=';'; /* replace ; */
      bps = bp+1;
      while(isspace(*bps)) bps++;
    }
  }
  else {
    fprintf(stderr," cannot find getdesc field:\n%s\n",tmp_str);
    goto error_r;
  }

  if ((bp=strchr(bps,';'))!=NULL) { *bp='\0';}

  if ((m_fptr->sql_getseq=calloc(strlen(bps)+1,sizeof(char)))==NULL) {
    fprintf(stderr, " cannot allocate space for database name [%d], %s\n",
	    (int)strlen(bps),bps);
    goto error_r;
  }

  if (strlen(bps) > 0) {
    strcpy(m_fptr->sql_getseq,bps);
    bps = bp+1;
  }
  else {
    fprintf(stderr," cannot find getseq field:\n%s\n",tmp_str);
    return 0;
  }
  if (bp!=NULL) *bp=';';

  /* check for close_table statement */
  if ((bp=strchr(bps,';'))!=NULL) {
    *bp='\0';
    if ((m_fptr->sql_close_tables=calloc(strlen(bps)+1,sizeof(char)))==NULL) {
      fprintf(stderr, " cannot allocate space for close_tables [%d], %s\n",
	      (int)strlen(bps),bps);
      goto error_r;
    }
    /* have get_desc, copy it */
    else {
      strcpy(m_fptr->sql_close_tables,bps);
      *bp=';'; /* replace ; */
      bps = bp+1;
      while(isspace(*bps)) bps++;
    }
  }

  /* now do the query */    

  if (mysql_query(m_fptr->mysql_conn,m_fptr->sql_query)) {
    fprintf(stderr,"*** Error %u - query failed:\n%s\n%s\n",
	    mysql_errno(m_fptr->mysql_conn), m_fptr->sql_query, mysql_error(m_fptr->mysql_conn));
    goto error_r;
  }

  if ((m_fptr->mysql_res = mysql_use_result(m_fptr->mysql_conn)) == NULL) {
    fprintf(stderr,"*** Error = use result failed\n%s\n",
	    mysql_error(m_fptr->mysql_conn));
    goto error_r;
  }
  return m_fptr;

 error_r:
  if (m_fptr->sql_close_tables) free(m_fptr->sql_close_tables);
  if (m_fptr->sql_getseq) free(m_fptr->sql_getseq);
  if (m_fptr->sql_getdesc) free(m_fptr->sql_getdesc);
  if (m_fptr->sql_query) free(m_fptr->sql_query);
  free(m_fptr);
  free(sql_db);
  return NULL;
}

struct lmf_str *
mysql_reopen(struct lmf_str *m_fptr) {
  m_fptr->sql_reopen = 1;
  return m_fptr;
}

void
mysql_closelib(struct lmf_str *m_fptr) {

  if (m_fptr == NULL) return;

  if (m_fptr->mysql_res != NULL)
    mysql_free_result(m_fptr->mysql_res);

  if (m_fptr->sql_close_tables) {
    if (mysql_query(m_fptr->mysql_conn,m_fptr->sql_close_tables)) {
      fprintf(stderr,"*** Error %u - close_tables failed:\n%s\n%s\n",
	      mysql_errno(m_fptr->mysql_conn), m_fptr->sql_close_tables,
	      mysql_error(m_fptr->mysql_conn));
    }
  }
  mysql_close(m_fptr->mysql_conn);
  m_fptr->sql_reopen=0;
}

/*
static char *sql_seq = NULL, *sql_seqp;
static int sql_seq_len;
static MYSQL_ROW sql_row;
*/

int
mysql_getlib( unsigned char *seq,
	      int maxs,
	      char *libstr,
	      int n_libstr,
	      fseek_t *libpos,
	      int *lcont,
	      struct lmf_str *lm_fd,
	      long *l_off)
{
  register unsigned char *cp, *seqp;
  register int *ap;
  unsigned char *seqm, *seqm1;
  char *bp;
  /*   int l_start, l_stop, len; */

  seqp = seq;
  seqm = &seq[maxs-9];
  seqm1 = seqm-1;

  ap = lm_fd->sascii;

  if (*lcont==0) {
    /* get a row, with UID, sequence */
    *l_off = 1;
    if ((lm_fd->mysql_row =mysql_fetch_row(lm_fd->mysql_res))!=NULL) {
      *libpos=(fseek_t)atol(lm_fd->mysql_row[0]);

      /* for @P:1-n removed */
      /*
      if ((bp=strchr(lm_fd->mysql_row[2],'@'))!=NULL &&
	  !strncmp(bp+1,"P:",2)) {
	sscanf(bp+3,"%d-%d",&l_start,&l_stop)
	l_start--;
	if (l_start < 0) l_start=0;
	if (l_stop > (len=strlen(lm_fd->mysql_row[1]))) l_stop= len-1;
	lm_fd->sql_seqp = lm_fd->mysql_row[1];
	lm_fd->sql_seqp[l_stop]='\0';
	lm_fd->sql_seqp += l_start;
      */

      if (lm_fd->mysql_row[2] == NULL) {
	fprintf(stderr," NULL comment at: [%s] %ld\n",
		lm_fd->mysql_row[0],*libpos);
      }
      else if ((bp=strchr(lm_fd->mysql_row[2],'@'))!=NULL &&
	  !strncmp(bp+1,"C:",2)) sscanf(bp+3,"%ld",l_off);
      else *l_off = 1;

      lm_fd->sql_seqp = lm_fd->mysql_row[1];

      /* because of changes in mysql_ranlib(), it is essential that
         libstr return the unique identifier; thus we must use
         sql_row[0], not sql_row[2]. Using libstr as the UID allows
         one to use any UID, not just numeric ones.  *libpos is not
         used for mysql libraries.
      */

      if (n_libstr <= MAX_UID) {
	/* the normal case returns only GID/sequence */
	strncpy(libstr,lm_fd->mysql_row[0],MAX_UID-1);
	libstr[MAX_UID-1]='\0';
      }
      else {
	/* here we do not use the UID in libstr, because we are not
           going back into the db */
	/* the PVM case also returns a long description */
	if (lm_fd->mysql_row[2]!=NULL) {
	  strncpy(libstr,lm_fd->mysql_row[2],n_libstr-1);
	}
	else {
	  strncpy(libstr,lm_fd->mysql_row[0],n_libstr-1);
	}
	libstr[n_libstr-1]='\0';
      }
    }
    else {
      mysql_free_result(lm_fd->mysql_res);
      lm_fd->mysql_res=NULL;
      *lcont = 0;
      *seqp = EOSEQ;
      return -1;
    }
  }

  for (cp=(unsigned char *)lm_fd->sql_seqp; seqp<seqm1 && *cp; ) {
    if ((*seqp++=ap[*cp++])<NA &&
	(*seqp++=ap[*cp++])<NA &&
	(*seqp++=ap[*cp++])<NA &&
	(*seqp++=ap[*cp++])<NA &&
	(*seqp++=ap[*cp++])<NA &&
	(*seqp++=ap[*cp++])<NA &&
	(*seqp++=ap[*cp++])<NA &&
	(*seqp++=ap[*cp++])<NA &&
	(*seqp++=ap[*cp++])<NA &&
	(*seqp++=ap[*cp++])<NA) continue;
    --seqp;
    if (*(cp-1)==0) break;
  }
  lm_fd->sql_seqp = (char *)cp;

  if (seqp>=seqm1) (*lcont)++;
  else {
    *lcont=0;
    if (lm_fd->sql_reopen) {
      mysql_free_result(lm_fd->mysql_res);
      lm_fd->mysql_res = NULL;
    }
  }

  *seqp = EOSEQ;
  /*   if ((int)(seqp-seq)==0) return 1; */
  return (int)(seqp-seq);
}

void
mysql_ranlib(char *str,
	     int cnt,
	     fseek_t libpos,
	     char *libstr,
	     struct lmf_str *lm_fd
	     )
{
  char tmp_query[1024], tmp_val[20];
  char *bp;

  str[0]='\0';

  /* put the UID into the query string - cannot use sprintf because of
     "%' etc */

  /*   sprintf(tmp_query,lm_fd->sql_getdesc,libpos); */

  if ((bp=strchr(lm_fd->sql_getdesc,'#'))==NULL) {
    fprintf(stderr, "no GID position in %s\n",lm_fd->sql_getdesc);
    goto next1;
  }
  else {
    *bp = '\0';
    strncpy(tmp_query,lm_fd->sql_getdesc,sizeof(tmp_query));
    tmp_query[sizeof(tmp_query)-1]='\0';
    /*    sprintf(tmp_val,"%ld",(long)libpos); */
    strncat(tmp_query,libstr,sizeof(tmp_query)-1);
    strncat(tmp_query,bp+1,sizeof(tmp_query)-1);
    *bp='#';
    lm_fd->lpos = libpos;
  }

  /*  fprintf(stderr," requesting: %s\n",tmp_query); */

  if (lm_fd->mysql_res !=NULL) {
    mysql_free_result(lm_fd->mysql_res);
    lm_fd->mysql_res = NULL;
  }

  if (mysql_query(lm_fd->mysql_conn,tmp_query)) {
    fprintf(stderr,"*** Error - query failed:\n%s\n%s\n",tmp_query,
	    mysql_error(lm_fd->mysql_conn));
    sprintf(str,"gi|%ld ***Error - query failed***",(long)libpos);
    goto next1;
  }

  if ((lm_fd->mysql_res = mysql_use_result(lm_fd->mysql_conn)) == NULL) {
/*     fprintf(stderr,"*** Error = use result failed\n%s\n", 
	   mysql_error(lm_fd->mysql_conn)); */
    sprintf(str,"gi|%ld ***use result failed***",(long)libpos);
    goto next0;
  }
  
  /* have the description */
  if ((lm_fd->mysql_row = mysql_fetch_row(lm_fd->mysql_res))==NULL) {
    /*    fprintf(stderr," cannot fetch description: %s\n",tmp_query); */
    sprintf(str,"gi|%ld ***cannot fetch description***",(long)libpos);
    goto next0;
  }
  
  if (lm_fd->mysql_row[1] != NULL) strncpy(str,lm_fd->mysql_row[1],cnt-1);
  else strncpy(str,lm_fd->mysql_row[0],cnt-1);
  str[cnt-1]='\0';
  while (strlen(str) < cnt-1 &&
	 (lm_fd->mysql_row = mysql_fetch_row(lm_fd->mysql_res))!=NULL) {
    strncat(str," ",cnt-2-strlen(str));
    if (lm_fd->mysql_row[1]!=NULL) 
      strncat(str,lm_fd->mysql_row[1],cnt-2-strlen(str));
    else break;
  }

  str[cnt-1]='\0';
  if ((bp = strchr(str,'\r'))!=NULL) *bp='\0';
  if ((bp = strchr(str,'\n'))!=NULL) *bp='\0';

 next0:
  mysql_free_result(lm_fd->mysql_res);
 next1: 
  lm_fd->mysql_res = NULL;

  /* get the sequence, set up for mysql_getseq() */
  /* put the UID into the query string */

  if ((bp=strchr(lm_fd->sql_getseq,'#'))==NULL) {
    fprintf(stderr, "no GID position in %s\n",lm_fd->sql_getseq);
    return;
  }
  else {
    *bp = '\0';
    strncpy(tmp_query,lm_fd->sql_getseq,sizeof(tmp_query));
    tmp_query[sizeof(tmp_query)-1]='\0';
    /*    sprintf(tmp_val,"%ld",(long)libpos); */
    strncat(tmp_query,libstr,sizeof(tmp_query));
    strncat(tmp_query,bp+1,sizeof(tmp_query));
    *bp='#';
  }

  if (mysql_query(lm_fd->mysql_conn,tmp_query)) {
    fprintf(stderr,"*** Error - query failed:\n%s\n%s\n",tmp_query,
	    mysql_error(lm_fd->mysql_conn));
  }

  if ((lm_fd->mysql_res = mysql_use_result(lm_fd->mysql_conn)) == NULL) {
    fprintf(stderr,"*** Error = use result failed\n%s\n",
	    mysql_error(lm_fd->mysql_conn));
  }
}
