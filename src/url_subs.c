/* $Id: url_subs.c $ */

/* copyright (c) 1998, 1999, 2014 by William R. Pearson and the
   The Rector & Visitors of the University of Virginia */

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

/* 30 Dec 2004 - modify REF_URL to accomodate current Entrez */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "defs.h"
#include "structs.h"
#include "param.h"

#ifndef DEF_PROT_LIB
#define DEF_PROT_LIB "q"
#endif

extern int seq_pos(int pos, int rev, int off);

char *display_domains(char, struct annot_entry **s_annot_arr_p, int n_domains);
char *web_encode(const char *);

void encode_json_str(FILE *fp, const char *label, const char *value, int first) {
  if (!first) {fprintf(fp, ",\n");}
  fprintf(fp, " \"%s\": \"%s\"",label, value);
}

void encode_json_long(FILE *fp, const char *label, long value, int first) {
  if (!first) {fprintf(fp, ",\n");}
  fprintf(fp, " \"%s\": %ld",label, value);
}

void encode_json_dfmt(FILE *fp, const char *label, double value, char *fmt, int first) {
  fprintf(fp, fmt, label, value);
}

void encode_json_aln(FILE *fp, const struct a_struct *aln_p, long q_offset, long l_offset, int first) {
}

void encode_json_lines(FILE *fp, const char *label, const char *annot_s, int first) {
  char *obp, *bp;

  char *tmp_annot_s;
  int n_tmp_annot_s;
  
  n_tmp_annot_s = strlen(annot_s)+1;
  if ((tmp_annot_s = (char *)calloc(n_tmp_annot_s,sizeof(char)))==NULL) {
    fprintf(stderr,"*** error [%s:%d] *** cannot allocate tmp_annot_s[%d]\n",
	    __FILE__, __LINE__,n_tmp_annot_s);
    return;
  }

  SAFE_STRNCPY(tmp_annot_s, annot_s, n_tmp_annot_s);

  if (!first) {fprintf(fp, ",\n");}
  fprintf(fp, " \"%s\": [\n",label);

  obp = bp = tmp_annot_s;
  while ((bp = strchr(obp,'\n'))) {
    *bp='\0';
    if (obp != tmp_annot_s) fprintf(fp, ",\n");
    fprintf(fp," \"%s\"",obp);
    obp = bp+1;
  }
  fprintf(fp, "\n ]");
  free(tmp_annot_s);
}

void encode_json_domains(FILE *fp, const char *label, const struct annot_str *annot_p, int first) {
  int i;

  if (!first) {fprintf(fp, ",\n");}
  fprintf(fp, "\"%s\": [\n",label);
  for (i=0; i < annot_p->n_annot; i++) {
    if (annot_p->s_annot_arr_p[i]->label != '-') continue;
    if (i != 0) fprintf(fp, ",\n");
    fprintf(fp, "  { \"start\":%ld, \"stop\":%ld, \"description\":\"%s\" }",
	    annot_p->s_annot_arr_p[i]->pos+1,annot_p->s_annot_arr_p[i]->end+1,annot_p->s_annot_arr_p[i]->comment);
  }
  fprintf(fp,"\n  ]");
}

void do_url1(FILE *fp, const struct mngmsg *m_msp, const struct pstruct *ppst,
	     char *l_name, int n1,
	     const struct a_struct *aln_p, const char *annot_var_s,
	     const struct annot_str *q_annot_p,
	     const struct annot_str *l_annot_p )
{
  char my_q_name[200], my_l_name[200], json_l_name[200];
  char *db, *bp;
  char pgm[10], o_pgm[10], lib[MAX_LSTR];
  char *tmp_annot_s, *q_domain_s, *l_domain_s, *tmp_domain_s, *etmp_domain_s;
  int  n_tmp_annot_s, n_tmp_domain;
  long q_offset, l_offset;
  char *ref_url, *lbp=NULL;
  char *srch_url, *srch_url1, *dom_url;

  /* set the database */
  if (m_msp->ldb_info.ldnaseq==SEQT_DNA) db="nucleotide";
  else db="Protein";

  /* set the program type */
  if (strncmp(m_msp->f_id0,"rss",3)==0) {
    strncpy(pgm,"fa",sizeof(pgm));
  }
  else if (strncmp(m_msp->f_id0,"rfx",3)==0) {
    strncpy(pgm,"fx",sizeof(pgm));
  }
  else { strncpy(pgm,m_msp->f_id0,sizeof(pgm)); }

  SAFE_STRNCPY(o_pgm, pgm, sizeof(o_pgm));

  /* get a library name (probably does not work for %, + abbreviations */
  if (m_msp->lname[0]!='%') {
    SAFE_STRNCPY(lib,m_msp->lname,sizeof(lib));
  }
  else {
    SAFE_STRNCPY(lib,"%25",sizeof(lib));
    SAFE_STRNCAT(lib,&m_msp->lname[1],sizeof(lib));
  }
  lib[sizeof(lib)-1]='\0';

  if ((lbp = strchr(l_name,'|'))==NULL) {
    lbp = l_name;
  }
  else {
    lbp++;
  }

  SAFE_STRNCPY(my_q_name,m_msp->qtitle,sizeof(my_q_name));
  if ((bp=strchr(my_q_name,' '))!=NULL) *bp='\0';

  SAFE_STRNCPY(my_l_name,lbp,sizeof(my_l_name));

  if (pgm[0]=='t' || !strcmp(pgm,"fx") || !strcmp(pgm,"fy") ) {
    if ((lbp=strchr(my_l_name,':'))!=NULL) *lbp='\0';
    lbp = &my_l_name[strlen(my_l_name)-2];
    if ( *lbp == '_' ) *lbp = '\0';
  }

  /* change the program name for fastx, tfastx, tfasta */
  /* fastx returns proteins */
  if (strcmp(pgm,"fx")==0 || strcmp(pgm,"fy")==0) {SAFE_STRNCPY(pgm,"fa",sizeof(pgm));}
  else if (strcmp(pgm,"ff")==0) {SAFE_STRNCPY(pgm,"fa",sizeof(pgm));}
  else if (pgm[0]=='t') {
    SAFE_STRNCPY(pgm,"fx",sizeof(pgm));
    SAFE_STRNCPY(lib,DEF_PROT_LIB,sizeof(lib));
  }

  fflush(fp);

  q_offset = aln_p->q_offset;
  l_offset = aln_p->l_offset;

  /* set up ref_url, srch_url, srch_url1, dom_url */

  fflush(fp);

  ref_url = getenv("REF_URL");
  srch_url = getenv("SRCH_URL");
  srch_url1 = getenv("SRCH_URL1");
  dom_url = NULL;
  dom_url = getenv("DOMAIN_PLOT_URL");

  if (ref_url || srch_url || srch_url1 || dom_url) {
    fprintf(fp,"<!-- LINK_START %s -->",l_name);

  /* REF_URL should provide */
  /* "<A HREF=\"http://www.ncbi.nlm.nih.gov/entrez/query.fcgi?db=%s&fcmd=Search&doptcmd1=DocSum&term=%s\">Entrez lookup</A>&nbsp;&nbsp;" */
  if (ref_url != NULL) {fprintf(fp,ref_url,db,my_l_name);}

  /* SRCH_URL should provide */
  /* "<A HREF=\"http://localhost/fasta_www2/searchfa.cgi?query=%s&db=fasta_www.cgi&lib=%s&pgm=%s&start=%ld&stop=%ld&n1=%d&o_pgm=%s\">Re-search database</A>&nbsp;&nbsp;" */
  if (srch_url != NULL) {
    fprintf(fp,srch_url,my_l_name,db,lib,pgm,
	    l_offset+aln_p->amin1+1,l_offset+aln_p->amax1,n1,m_msp->f_id0);
  }

  /* SRCH_URL1 should provide: */
  /*  "<A HREF=\"http://localhost/fasta_www2/searchxf.cgi?query=%s&db=%s&lib=%s&pgm=%s&start=%ld&stop=%ld&n1=%d&o_pgm=%s\">General re-search</A>\n" */

  if (srch_url1 != NULL) {
    fprintf(fp,srch_url1,my_l_name,db,lib,pgm,
	    l_offset+aln_p->amin1+1,l_offset+aln_p->amax1,n1,m_msp->f_id0);
  }
  
  if (dom_url!=NULL) {
    if (annot_var_s && annot_var_s[0]) {
      tmp_annot_s = web_encode(annot_var_s);
    }
    else tmp_annot_s = "";

    q_domain_s = l_domain_s = NULL;

    if (q_annot_p && q_annot_p->n_domains > 0 && 
	(q_domain_s = display_domains('q',q_annot_p->s_annot_arr_p, q_annot_p->n_annot))!=NULL) {
    }
    if (l_annot_p && l_annot_p->n_domains > 0 && 
	(l_domain_s = display_domains('l',l_annot_p->s_annot_arr_p, l_annot_p->n_annot))!=NULL) {
    }

    /* combine domain strings */
    n_tmp_domain = 0;
    if (q_domain_s) n_tmp_domain += strlen(q_domain_s)+1;
    if (l_domain_s) n_tmp_domain += strlen(l_domain_s)+1;
    etmp_domain_s = "";
    if (n_tmp_domain > 0) {
      if ((tmp_domain_s=(char *)calloc(n_tmp_domain,sizeof(char)))==NULL) {
	fprintf(stderr,"*** error [%s:%d] *** cannot allocate tmp_domain_s[%d]\n",
		__FILE__, __LINE__,n_tmp_domain);
      }
      else {
	tmp_domain_s[0] = '\0';
	if (q_domain_s) SAFE_STRNCAT(tmp_domain_s, q_domain_s, n_tmp_domain);
	if (l_domain_s) SAFE_STRNCAT(tmp_domain_s, l_domain_s, n_tmp_domain);
	etmp_domain_s = web_encode(tmp_domain_s);
      }
    }

    /* appropriate format string: */
    /* 
       pgm=%s	    -- program abbrev that created alignment
       q_name=%s     -- query info
       q_cstart=%ld
       q_cstop=%ld
       q_astart=%ld
       q_astop=%ld
       l_name=%s     -- library info
       l_cstart=%ld
       l_cstop=%ld
       l_astart=%ld
       l_astop=%ld
       region=%s       -- aligned domain and variant information
       doms=%s

       DOMAIN_PLOT_URL = "pgm=%s;q_name=%s;q_cstart=%ld;q_cstop=%ld&q_astart=%ld&q_astop=%ld&l_name=%s&l_cstart=%ld&l_cstop=%ld&l_astart=%ld&l_astop=%ld&regions=%s&doms=%s"
    */

    /* think about the alternative of running a script
       rather than embedding it */

    fprintf(fp,dom_url,o_pgm,
	    my_q_name, q_offset+seq_pos(1,aln_p->qlrev,2),q_offset+seq_pos(m_msp->n0,aln_p->qlrev,2),
	    q_offset+seq_pos(aln_p->amin0+1,aln_p->qlrev,1), q_offset+seq_pos(aln_p->amax0, aln_p->qlrev,2),
	    my_l_name, l_offset+seq_pos(1,aln_p->llrev,2), l_offset+seq_pos(n1,aln_p->llrev,2),
	    l_offset+seq_pos(aln_p->amin1+1,aln_p->llrev,1),l_offset+seq_pos(aln_p->amax1,aln_p->llrev,2),
	    tmp_annot_s, etmp_domain_s);

    if (n_tmp_domain>0 && tmp_domain_s) {
      free(tmp_domain_s);
      free(etmp_domain_s);
    }
    if (l_annot_p && l_annot_p->n_domains && l_domain_s) {
      free(l_domain_s);
    }
    if (q_annot_p && q_annot_p->n_domains && q_domain_s) {
      free(q_domain_s);
    }
    if (annot_var_s && annot_var_s[0] && tmp_annot_s) free(tmp_annot_s);
  }

  fprintf(fp,"\n<!-- LINK_STOP -->");
  fflush(fp);
  }

  /*
    if ((srch_url2 = getenv("SRCH_URL2"))==NULL)
    fprintf(fp,"<A HREF=\"http://fasta.bioch.virginia.edu/fasta/cgi/lalignx.cgi?seq1=\"%s\"&in_seq1=\"FASTA\"&seq2=\"%s\"&in_seq2=\"Accession\"&ssr2=%ld:%ld\">lalign</A>\n<p>\n",my_l_name,db,lib,pgm,l_offset+aln_p->amin1+1,l_offset+aln_p->amax1,n1);
    else 
    fprintf(fp,srch_url1,my_l_name,db,lib,pgm,
    l_offset+aln_p->amin1+1,l_offset+aln_p->amax1,n1);
  */


  if (getenv("JSON_HTML")) {

    /* replace '|' with '_' */
    SAFE_STRNCPY(json_l_name, l_name, sizeof(json_l_name));
    for (bp=strchr(json_l_name,'|'); bp; bp=strchr(bp+1,'|')) { *bp = '_'; }

    /* replace '.' with '_' */
    for (bp=strchr(json_l_name,'.'); bp; bp=strchr(bp+1,'.')) { *bp = '_'; }

    fprintf(fp,"\n<script type=\"text/javascript\">\n//<![CDATA[\n var json_%s = {\n",json_l_name);
    encode_json_str(fp, "db", db, 1);
    encode_json_str(fp, "l_acc", l_name, 0);
    encode_json_str(fp, "acc", my_l_name, 0);
    encode_json_str(fp, "lib", lib, 0);
    encode_json_str(fp, "pgm", pgm, 0);
    encode_json_str(fp, "o_pgm", m_msp->f_id0, 0);
    encode_json_aln(fp, aln_p, q_offset, l_offset, 0);
    if (annot_var_s && annot_var_s[0]) { encode_json_lines(fp, "annot", annot_var_s, 0); }
    if (q_annot_p && q_annot_p->n_domains > 0) { encode_json_domains(fp, "q_domains", q_annot_p, 0); }
    if (l_annot_p && l_annot_p->n_domains > 0) { encode_json_domains(fp, "l_domains", l_annot_p, 0); }

    fprintf(fp, "\n}\n//]]>\n</script>");
    fflush(fp);
  }
}

char *display_domains(char target, struct annot_entry **annot_arr_p, int n_annots) {
  char *domain_s;
  char line[MAX_STR];
  int i, i_doms, n_domain_s = MAX_LSTR;

  /* since (currently) annot_var_s is MAX_LSOTR, do the same for domain_s */
  if ((domain_s = (char *)calloc(n_domain_s, sizeof(char)))==NULL) {
    fprintf(stderr,"*** error [%s:%d] *** cannot allocate domain_s[%d]\n",__FILE__, __LINE__,n_domain_s);
    return NULL;
  }

  for (i=0; i < n_annots; i++) {
    /* annot_arr_p[] has both domains and non domains, but n_domains only counts domains */
    if (annot_arr_p[i]->label != '-') continue;
    sprintf(line, "%cDomain:\t%ld-%ld\t%s\n",
	    target, annot_arr_p[i]->pos+1, annot_arr_p[i]->end+1, annot_arr_p[i]->comment);
    if (strlen(domain_s) + strlen(line)+1 > n_domain_s) {
      n_domain_s += n_domain_s/2;
      domain_s = realloc(domain_s, n_domain_s);
    }
    SAFE_STRNCAT(domain_s, line, n_domain_s);
  }

  domain_s = realloc(domain_s, (n_domain_s=strlen(domain_s))+1);
  domain_s[n_domain_s]='\0';

  return domain_s;
}

/* take an annotation string *annot_var_s and convert problematic characters to their web encoding */
/* ' ' (space) %20 */
/* '|' 	    %7C */
/* ';'	    %3B */
/* '='	    %3D */
/* '\n'	    %0A */

static char bad_chars[] = "\n =;|";

char *web_encode(const char *annot_var_s) {
  
  int n_tmp_annot_s;
  char *tmp_annot_s, *tmp_annot_d, *dp;
  const char *bp, *sp;
  int bad_cnt = 0;
  
  /* make string largest possible size */
  n_tmp_annot_s = strlen(annot_var_s)*3 + 1;
  if ((tmp_annot_s = (char *)calloc(n_tmp_annot_s,sizeof(char)))==NULL) {
    fprintf(stderr,"*** error [%s:%d] *** cannot allocate tmp_annot_s[%d]\n",__FILE__, __LINE__,n_tmp_annot_s);
    return NULL;
  }

  dp = tmp_annot_s;
  for (sp = annot_var_s; *sp ; sp++) {

    if ((*sp < '0') ||
	(*sp > 9 &&  *sp < 'A') ||
	(*sp > 'Z' &&  *sp < 'a') ||
	(*sp > 'z')) { sprintf(dp,"%%%02x",*sp); dp += 3;}
    else { *dp++ = *sp; }
  }

  n_tmp_annot_s = dp - tmp_annot_s;
  tmp_annot_s = realloc(tmp_annot_s, n_tmp_annot_s+1);
  tmp_annot_s[n_tmp_annot_s] = '\0';

  return tmp_annot_s;
}
