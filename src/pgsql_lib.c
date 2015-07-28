/*  $Id: pgsql_lib.c 781 2011-06-20 10:31:40Z wrp $ */

/* pgsql_lib.c copyright (c) 2004, 2014 by William R. Pearson and
   The Rector & Visitors of the University of Virginia */
/*
     Licensed under the Apache License, Version 2.0 (the "License");
     you may not use this file except in compliance with the License.
     You may obtain a copy of the License at

       http://www.apache.org/licenses/LICENSE-2.0

     Unless required by applicable law or agreed to in writing,
     software distributed under this License is distributed on an "AS
     IS" BASIS, WITHOUT WRRANTIES OR CONDITIONS OF ANY KIND, either
     express or implied.  See the License for the specific language
     governing permissions and limitations under the License. 
*/

/* functions for opening, reading, seeking a pgsql database */

/*
  For the moment, this interface assumes that the file to be searched will
  be specified in a single, long, string with 4 parts:

  (1) a database open string. This string has four fields, separated by
      whitespace (' \t'):
        hostname:port dbname user password

   '--' dashes at the beginning of lines are ignored -
   thus the first line could be:
   -- hostname:port dbname user password

  (2) a database query string that will return an unique ID (not
      necessarily numberic, but it must be < 12 characters as libstr[12]
      is used) and a sequence string

  (2a) a series of pgsql commands that do not generate results
       starting with 'DO', followed by a select() statement.

  (3) a database select string that will return a description
      given a unique ID

  (4) a database select string that well return a sequence given a
      unique ID

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
   The last (fourth) query must terminate with a ';'

   19-July-2004

   This file is designed for PostgreSQL, which uses a different syntax
   for getting rows of data.  Specifically, a select statement must be
   associated with a "cursor", so that one can fetch a single row.

   This can be simply done with the statment:

   DECLARE next_seq CURSOR FOR "select statement ..."

   The need for a CURSOR complicates the getlib()/ranlib() design, which
   assumes that ranlib() can set something up that getlib() can read.
   This can be avoided by setting up an otherwise unnecessary cursor for
   the ranlib statement that gets a sequence.

*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>

#include <libpq-fe.h>
#define PGSQL_LIB 17

#include "defs.h"
#include "mm_file.h"

#define XTERNAL
#include "uascii.h"
#define EOSEQ 0
/* #include "upam.h" */

char *alloc_file_name(char *f_name);
int pgsql_getlib(unsigned char *, int, char *, int, fseek_t *, int *, struct lmf_str *, long *);
void pgsql_ranlib(char *, int, fseek_t, char *, struct lmf_str *m_fd);

#define PGSQL_BUF 4096

struct lmf_str *
pgsql_openlib(char *sname, int ldnaseq, int *sascii) {
  FILE *sql_file;
  PGconn *conn;
  PGresult *res;
  char *tmp_str, *ttmp_str;
  int tmp_str_len;
  char *bp, *bps, *bdp, *tp, tchar;
  int i, qs_len, qqs_len;
  char *sql_db, *sql_host, *sql_dbname, *sql_user, *sql_pass;
  char *sql_port;
  char *sql_do;
  int sql_do_cnt;
  struct lmf_str *m_fptr;

  /*  if (sql_reopen) return NULL; - should not be called for re-open */

  tmp_str_len = PGSQL_BUF;
  if ((tmp_str=(char *)calloc(tmp_str_len,sizeof(char)))==NULL) {
    fprintf(stderr,"cannot allocate %d for pgSQL buffer\n",tmp_str_len);
    return NULL;
  }

  if (sname[0] == '%') {
    strncpy(tmp_str,sname+1,tmp_str_len);
    tmp_str[sizeof(tmp_str)-1]='\0';
  }
  else {
    if ((sql_file=fopen(sname,"r"))==NULL) {
      fprintf(stderr," cannot open pgSQL file: %s\n",sname);
      return NULL;
    }

    if ((qs_len=fread(tmp_str,sizeof(char),tmp_str_len-1,sql_file))<=0) {
      fprintf(stderr," cannot read pgSQL file: %s\n",sname);
      return NULL;
    }
    else  {
      tmp_str[qs_len]='\0';
      qqs_len = qs_len;
      while (qqs_len >= tmp_str_len-1) {
	tmp_str_len += PGSQL_BUF;
	if ((tmp_str=(char *)realloc(tmp_str,tmp_str_len))==NULL) {
	  fprintf(stderr,
		  " cannot reallocate %d for pgSQL buffer\n",tmp_str_len);
	  return NULL;
	}
	ttmp_str = &tmp_str[qqs_len];
	if ((qs_len=fread(ttmp_str,sizeof(char),PGSQL_BUF,sql_file))<0) {
	  fprintf(stderr," cannot read pgSQL file: %s\n",sname);
	  return NULL;
	}
	ttmp_str[qs_len]='\0';
	qqs_len += qs_len;
      }
    }
    fclose(sql_file);
  }

  bps = tmp_str;
  if ((bp=strchr(bps,';'))!=NULL) {
    *bp='\0';
    if ((sql_db=calloc(strlen(bps)+1,sizeof(char)))==NULL) {
      fprintf(stderr, " cannot allocate space for database name [%d], %s\n",
	      (int)strlen(bps),bps);
      return NULL;
    }
    /* have database name, parse the fields */
    else {
      strcpy(sql_db,bps);	/* strcpy OK because allocated strlen(bps) */
      bps = bp+1;	/* points to next char after ';' */
      while (isspace(*bps)) bps++;
      *bp=';'; /* replace ; */
      bp = sql_db;
      while (*bp=='-') {*bp++ = ' ';}
      sql_host = strtok(bp," \t\n");
      if (sql_host[0]=='@') sql_host="";
      sql_dbname = strtok(NULL," \t\n");
      sql_user = strtok(NULL," \t\n");
      if (sql_user[0]=='@') sql_user="";
      sql_pass = strtok(NULL," \t\n");
      if (sql_pass[0]=='@') sql_pass="";
      if ((tp=strchr(sql_host,':'))!=NULL) {
	sql_port = tp+1;
	*tp='\0';
      }
      else sql_port = "";
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
  m_fptr->getlib = pgsql_getlib;
  m_fptr->ranlib = pgsql_ranlib;
  m_fptr->mm_flg = 0;
  m_fptr->sql_reopen = 0;
  m_fptr->lb_type = PGSQL_LIB;

  /* now open the database, if necessary */
  conn = PQsetdbLogin(sql_host,
		      sql_port,
		      NULL,
		      NULL,
		      sql_dbname,
		      sql_user,
		      sql_pass);

  if (PQstatus(conn) != CONNECTION_OK)     {
    fprintf(stderr, "Connection to database '%s' failed.\n", PQdb(conn));
    fprintf(stderr, "%s", PQerrorMessage(conn));
    PQfinish(conn);
    goto error_r;
  }
  else {
    m_fptr->pgsql_conn = conn;
    fprintf(stderr," Database %s opened on %s\n",sql_dbname,sql_host);
  }

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
      /* do the pgsql statement on bdp+3 */
      /* check for error */
      *bp='\0';
      res = PQexec(m_fptr->pgsql_conn,bdp+3);
      if (PQresultStatus(res) != PGRES_COMMAND_OK) {
	fprintf(stderr,"*** Error %s - query failed:\n%s\n",
		PQerrorMessage(m_fptr->pgsql_conn), bdp+3);
	PQclear(res);
	goto error_r;
      }
      PQclear(res);

      *bp=';';
      bdp = bp+1;
      while (isspace(*bdp)) bdp++;
    }
  }

  /* copy 1st query field */
  if ((bp=strchr(bps,';'))!=NULL) {
    *bp='\0';
    if ((m_fptr->sql_query=calloc(strlen(bps)+41,sizeof(char)))==NULL) {
      fprintf(stderr, " cannot allocate space for query string [%d], %s\n",
	      (int)strlen(bps),bps);
      goto error_r;
    }
    /* have query, copy it */
    else {
      strncpy(m_fptr->sql_query,"DECLARE next_seq CURSOR FOR ",40);
      strcat(m_fptr->sql_query,bps);
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
  }
  else {
    fprintf(stderr," cannot find getseq field:\n%s\n",tmp_str);
    return NULL;
  }
  if (bp!=NULL) *bp=';';

  /* now do the fetch */

  res = PQexec(m_fptr->pgsql_conn,"BEGIN;");
  if (PQresultStatus(res) != PGRES_COMMAND_OK) {
    fprintf(stderr,"*** Error %s - BEGIN failed:\n",
	    PQerrorMessage(conn));
    PQclear(res);
    goto error_r;
  }
  PQclear(res);

  res = PQexec(m_fptr->pgsql_conn, m_fptr->sql_query);
  if (PQresultStatus(res) != PGRES_COMMAND_OK) {
    fprintf(stderr,"*** Error %d:%s - query failed:\n%s\n",
	    PQresultStatus(res),PQerrorMessage(conn), m_fptr->sql_query);
    PQclear(res);
    goto error_r;
  }
  PQclear(res);
  m_fptr->pgsql_res=NULL;

  return m_fptr;

 error_r:
  free(m_fptr->sql_getseq);
  free(m_fptr->sql_getdesc);
  free(m_fptr->sql_query);
  free(m_fptr);
  free(sql_db);
  return NULL;
}

struct lmf_str *
pgsql_reopen(struct lmf_str *m_fptr) {
  m_fptr->sql_reopen = 1;
  return m_fptr;
}

void
pgsql_closelib(struct lmf_str *m_fptr) {

  if (m_fptr == NULL) return;
  if (m_fptr->pgsql_res != NULL) PQclear(m_fptr->pgsql_res);
  PQfinish(m_fptr->pgsql_conn);
  m_fptr->sql_reopen=0;
}

/*
static char *sql_seq = NULL, *sql_seqp;
static int sql_seq_len;
*/

int
pgsql_getlib( unsigned char *seq,
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
  PGresult *res;

  char *bp;
  /*   int l_start, l_stop, len; */

  seqp = seq;
  seqm = &seq[maxs-9];
  seqm1 = seqm-1;

  ap = lm_fd->sascii;

  if (*lcont==0) {
    /* get a row, with UID, sequence */
    *l_off = 1;

    /* check to see if we already have a valid result */
    if (lm_fd->pgsql_res==NULL) {
      res = PQexec(lm_fd->pgsql_conn,"FETCH next_seq");
      if (PQresultStatus(res) != PGRES_TUPLES_OK) {
	fprintf(stderr,"*** Error %s - getlib FETCH failed:\n%s\n",
		PQerrorMessage(lm_fd->pgsql_conn), lm_fd->sql_query);
	PQclear(res);
	lm_fd->pgsql_res = NULL;
	*lcont = 0;
	*seqp = EOSEQ;
	return -1;
      }
    }
    else {res = lm_fd->pgsql_res;}

    if (PQntuples(res)>0) {
      lm_fd->pgsql_res = res;
      *libpos=(fseek_t)atol(PQgetvalue(res,0,0));
	
      *l_off = 1;
      if (PQnfields(res) > 2 && (bp=strchr(PQgetvalue(res,0,2),'@'))!=NULL &&
	  !strncmp(bp+1,"C:",2)) sscanf(bp+3,"%ld",l_off);

      lm_fd->sql_seqp = PQgetvalue(res,0,1);
    
      /* because of changes in pgsql_ranlib(), it is essential that
         libstr return the unique identifier; thus we must use
         sql_row[0], not sql_row[2]. Using libstr as the UID allows
         one to use any UID, not just numeric ones.  *libpos is not
         used for pgsql libraries.
      */

      if (n_libstr <= MAX_UID) {
	/* the normal case returns only GID/sequence */
	strncpy(libstr,PQgetvalue(res,0,0),MAX_UID-1);
	libstr[MAX_UID-1]='\0';
      }
      else {
	/* here we do not use the UID in libstr, because we are not
           going back into the db */
	/* the PVM case also returns a long description */
	if (PQnfields(res)>2) {
	  strncpy(libstr,PQgetvalue(res,0,2),n_libstr-1);
	}
	else {
	  strncpy(libstr,PQgetvalue(res,0,0),n_libstr-1);
	}
	libstr[n_libstr-1]='\0';
      }
    }
    else {
      PQclear(lm_fd->pgsql_res);
      lm_fd->pgsql_res=NULL;
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
    PQclear(lm_fd->pgsql_res);
    lm_fd->pgsql_res = NULL;
  }

  *seqp = EOSEQ;
  /*   if ((int)(seqp-seq)==0) return 1; */
  return (int)(seqp-seq);
}

void
pgsql_ranlib(char *str,
	     int cnt,
	     fseek_t libpos,
	     char *libstr,
	     struct lmf_str *lm_fd
	     )
{
  char tmp_query[1024], tmp_val[20];
  PGresult *res;
  char *bp;

  str[0]='\0';

  /* put the UID into the query string - cannot use sprintf because of
     "%' etc */

  /*   sprintf(tmp_query,lm_fd->sql_getdesc,libpos); */

  if ((bp=strchr(lm_fd->sql_getdesc,'#'))==NULL) {
    fprintf(stderr, "no KEY position in %s\n",lm_fd->sql_getdesc);
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

  if (lm_fd->pgsql_res !=NULL) {
    PQclear(lm_fd->pgsql_res);
    lm_fd->pgsql_res = NULL;
  }

  res = PQexec(lm_fd->pgsql_conn,tmp_query);
  if (PQresultStatus(res) != PGRES_TUPLES_OK) {
    lm_fd->pgsql_res = NULL;

    sprintf(str,"gi|%ld ***Error - query failed***",(long)libpos);
    fprintf(stderr,"*** Error %s - ranlib DESC failed:\n%s\n",
	    PQerrorMessage(lm_fd->pgsql_conn), tmp_query);
    PQclear(res);
    goto next1;
  }

  if (PQntuples(res)<=0) {
/*     fprintf(stderr,"*** Error = use result failed\n%s\n", 
	   pgsql_error(lm_fd->pgsql_conn)); */
    sprintf(str,"gi|%ld ***use result failed***",(long)libpos);
    goto next0;
  }

  if (PQgetvalue(res,0,1)!= NULL) strncpy(str,PQgetvalue(res,0,1),cnt-1);
  else strncpy(str,PQgetvalue(res,0,0),cnt-1);
  str[cnt-1]='\0';
  /* change this later to support multiple row returns */
  /*
  while (strlen(str) < cnt-1 &&
	 (lm_fd->sql_row = pgsql_fetch_row(lm_fd->pgsql_res))!=NULL) {
    strncat(str," ",cnt-2-strlen(str));
    if (lm_fd->sql_row[1]!=NULL) 
      strncat(str,lm_fd->sql_row[1],cnt-2-strlen(str));
    else break;
  }
  */

  str[cnt-1]='\0';
  if ((bp = strchr(str,'\r'))!=NULL) *bp='\0';
  if ((bp = strchr(str,'\n'))!=NULL) *bp='\0';

 next0:
  PQclear(res);
 next1: 
  lm_fd->pgsql_res = NULL;

  /* get the sequence, set up for pgsql_getseq() */
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

  res = PQexec(lm_fd->pgsql_conn,tmp_query);
  if (PQresultStatus(res) != PGRES_TUPLES_OK) {
    PQclear(res);
    lm_fd->pgsql_res = NULL;
    fprintf(stderr,"*** Error - ranlib SEQ failed:\n%s\n%s\n",tmp_query,
	    PQerrorMessage(lm_fd->pgsql_conn));
    exit(1);
  }
  else {
    lm_fd->pgsql_res = res;
  }
}
