/* $Id: lib_sel.c 792 2011-06-26 18:20:30Z wrp $ */

/* copyright (c) 1996, 1997, 1998, 1999, 2014 by William R. Pearson and
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

/*	modified Dec 13, 1989 requires different FASTLIBS */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <ctype.h>
#include <string.h>

#include "defs.h"
#include "structs.h"

#ifdef NCBIBL13
#define LASTLIB NCBIBL13+1
#else
#define LASTLIB 11
#endif


struct lib_struct *get_lnames(char *tname, struct lib_struct *cur_lib_p);
struct lib_struct *add_file(char *name, char *env, struct lib_struct *cur_lib_p);
void lib_choice(char *lname, int nl, char *flstr, int ldnaseq);
void subs_env(char *dest, char *src, int dest_size);
char *ulindex(char *str, char *chr);

static char ldname[MAX_FN];
static char *libenv;

/* read in the library names.
   returns the beginning of the list of names, not the end

   if cur_lib_p is NULL, then allocates it and returns it.
   if cur_lib_p is not NULL, then links to cur_lib_p->next and returns cur_lib_p
*/
struct lib_struct *
get_lnames(char *iname, struct lib_struct *cur_lib_p)
{
  char *bp, tsave[MAX_STR], *tname;
  char lline[MAX_FN], *llp;
  struct lib_struct *new_lib_p;
  FILE *tptr;

  /* expand environment variables */
  
  tname = tsave;
  subs_env(tname, iname, sizeof(tsave));

  if (*tname != '@') {
    new_lib_p = add_file(tname,"\0",cur_lib_p);
    if (cur_lib_p == NULL) return new_lib_p;
    else return cur_lib_p;
  }
  else tname++;

  /* remove ' ' before deftype if present */
  if ((bp=strchr(tname,' '))!=NULL) *bp='\0';

  if ((tptr=fopen(tname,"r"))==NULL) {
    fprintf(stderr,"*** ERROR [%s:%d] could not open file of names: %s\n",__FILE__,__LINE__,tname);
    return NULL;
  }

  new_lib_p = cur_lib_p;
  while (fgets(lline,sizeof(lline),tptr)!=NULL) {
    if (lline[0]==';') continue;
    if ((bp=strchr(lline,'\n'))!=NULL) *bp='\0';
    subs_env(tsave, lline, sizeof(tsave));
    if (tsave[0]=='<') {
      strncpy(ldname,&tsave[1],sizeof(ldname));
      ldname[sizeof(ldname)-1]='\0';
      libenv=ldname;
    }
    else {
      new_lib_p = add_file(tsave,libenv,new_lib_p);
      if (cur_lib_p == NULL) cur_lib_p = new_lib_p;
    }
  }
  fclose(tptr);
  return cur_lib_p;
}

void
lib_choice(char *lname, int nl, char *flstr, int ldnaseq)
{
  FILE *fch;
  char line[MAX_STR], *bp;
  char *chstr[MAX_CH],*chfile[MAX_CH];
  char *chtmp, *charr;
  int i,j,k,chlen;

  charr = NULL;
  if (strlen(flstr)> (size_t)0) {
    chlen = MAX_CH*MAX_FN;
    if ((chtmp=charr=calloc((size_t)chlen,sizeof(char)))==NULL) {
      fprintf(stderr,"*** ERROR [%s:%d] cannot allocate choice file array\n",__FILE__,__LINE__);
      goto l1;
    }
    chlen--;
    if ((fch=fopen(flstr,"r"))==NULL) {
      fprintf(stderr,"*** ERROR [%s:%d] cannot open choice file: %s\n",__FILE__,__LINE__,flstr);
      goto l1;
    }
    fprintf(stderr,"\n Choose sequence library:\n\n");

    for (i=j=0; j<MAX_CH; i++) {
      if (fgets(line,sizeof(line),fch)==NULL) break;/* check for comment */
      if (line[0]==';') continue;
      if ((bp=strchr(line,'\n'))!=NULL) *bp='\0'; /* remove \n */
      if ((bp=strchr(line,'$'))==NULL) continue;  /* if no '$', continue */
      *bp++='\0';	      /* replace $ with \0, bp points to libtype */

      /* if libtypes don't match, continue */
      if ((*bp++ -'0')!=ldnaseq) continue;

      /* if the library file name is too long, quit */
      if ((k=strlen(line))>chlen) break;

      /* save the library file name */
      strncpy(chstr[j]=chtmp,line,chlen);
      chtmp += k+1; chlen -= k+1;

      if ((k=strlen(bp))>chlen) break;
      strncpy(chfile[j]=chtmp,bp,chlen);
      chtmp += k+1; chlen -= k+1;
      fprintf(stderr,"    %c: %s\n",*chfile[j++],line);
    }
  l2:  fprintf(stderr,"\n Enter library filename (e.g. %s), letter (e.g. P)\n",
	       (ldnaseq==0)? "prot.lib" : "dna.lib");
    fprintf(stderr," or a %% followed by a list of letters (e.g. %%PN): ");
    fflush(stderr);
    if (fgets(line,sizeof(line),stdin)==NULL) exit(0);
    if ((bp=strchr(line,'\n'))!=NULL) *bp='\0';
    if (strlen(line)==0) goto l2;
    strncpy(lname,line,nl);
  }
  else {
  l1: fprintf(stderr," library file name: ");
    fflush(stderr);
    if (fgets(line,sizeof(line),stdin)==NULL) exit(0);
    if ((bp=strchr(line,'\n'))!=NULL) *bp='\0';
    if (strlen(line)> (size_t)0) strncpy(lname,line,nl);
    else goto l1;
  }
  if (charr!=NULL) {
    fclose(fch);
    free(charr);
  }
}

/* lib_select parses the choices in char *lname and builds the list
   of library files.

   lib_select returns the head of the list of files

   lib_select recognizes:
   %prm  -- leading '%' indicates a list of one letter abbreviations
   +abrev1+abrev2+abrev3 -- leading '+' indicates a list of word abbreviations

   @library.fil -- leading '@' indicates a list of library files
   (possibly containing additional @files)

   Support for NCBI .pal/.nal files needs to be added
*/
struct lib_struct *
lib_select(char *lname, char *ltitle, const char *flstr, int ldnaseq)
{
  char line[MAX_FN*2], *bp, *bp1;
  char *llnames[MAX_LF]; /* pointers into new list of names */
  int new_abbr,ich, nch;	  /* use new multi-letter abbr */
  int ltmp;
  FILE *fch;
  struct lib_struct *cur_lib_p = NULL, *tmp_lib_p;

  new_abbr = 0;
  *ltitle = '\0';

  if (strlen(lname) > (size_t)1 && *lname != '%' && *lname != '+') {
    return get_lnames(lname, cur_lib_p); /* file name */ 
  }
  else {
    if (*flstr=='\0') {
      fprintf(stderr,"*** ERROR [%s:%d] abbrv. list request but FASTLIBS undefined, cannot use %s\n",__FILE__,__LINE__,lname);
      exit(1);
    }

    if (strchr(lname,'+')) {
      /* indicates list of database abbrevs (not files) */
      new_abbr=1;
      nch = 0;
      bp = lname+1; if (*bp == '+') bp++;
      for (bp1=bp; nch < MAX_LF && bp!=NULL && bp1!=NULL; bp=bp1+1) {
	if ((bp1=strchr(bp,'+'))!=NULL) *bp1='\0';
	llnames[nch++] = bp;
      }
    }
    else if (*lname=='%') {     /* list of single letter abbreviations */
      lname++;	/* bump over '%' to get letters */
    }

    /* else just use a single character abbreviation */

    if (strlen(flstr) > (size_t)0) {
      if ((fch=fopen(flstr,"r"))==NULL) {
	fprintf(stderr,"*** ERROR [%s:%d] cannot open choice file: %s\n",__FILE__,__LINE__,flstr);
	return NULL;
      }
    }

    /* read each line of FASTLIBS */
    while (fgets(line,sizeof(line),fch)!=NULL) { 
      if (line[0]==';') continue;	/* skip comments */
      if ((bp=strchr(line,'\n'))!=NULL) *bp='\0';	/* remove '\n' */
      if ((bp=strchr(line,'$'))==NULL) continue; /* no delim, continue */
      *bp++='\0';	/* point to library type */
      if ((*bp++ -'0')!=ldnaseq) continue; /* doesn't match, continue */

      /* if !new_abbr, match on one letter with ulindex() */
      if (!new_abbr) {
	if (*bp=='+') continue; /* not a +lib+ */
	else if (ulindex(lname,bp)!=NULL) { 
	  if (ltitle[0] == '\0') {
	    strncpy(ltitle,line,MAX_STR);
	  }
	  else {
	    ltmp = strlen(ltitle);
	    strncat(ltitle,",\n  ",MAX_STR-ltmp);
	    strncat(ltitle,line,MAX_STR-ltmp-4);
	  }
	  tmp_lib_p = get_lnames(bp+1, cur_lib_p);
	  if (tmp_lib_p) { cur_lib_p = tmp_lib_p;}
	}
      }
      else {
	if (*bp!='+') continue;
	else {
	  bp++;
	  if ((bp1 = strchr(bp,'+'))!=NULL) {
	    *bp1='\0';
	    for (ich = 0; ich<nch; ich++) {
	      if (strcmp(llnames[ich],bp)==0) {
		if (ltitle[0] == '\0') {
		  strncpy(ltitle,line,MAX_STR);
		}
		else {
		  ltmp = strlen(ltitle);
		  strncat(ltitle,",\n  ",MAX_STR-ltmp);
		  strncat(ltitle,line,MAX_STR-ltmp-4);
		}
		cur_lib_p = get_lnames(bp1+1, cur_lib_p);
		break;
	      }
	    }
	    *bp1='+';
	  }
	  else fprintf(stderr,"*** ERROR [%s:%d] %s missing final '+'\n",__FILE__,__LINE__,bp);
	}
      }
    }
    fclose(fch);
  }
  return cur_lib_p;
}

/* unlike lib_select() and get_lnames(), add_file() returns a new
   pointer, to which library files can be added
*/
struct lib_struct *
add_file(char *fname, char *env, struct lib_struct *cur_lib_p)
{
  char tname[MAX_STR], *bp, *bp1;
  char *lbptr;
  int len, lenv, l_size;
  struct lib_struct *this_lib_p;

  /*  check for default directory for files  */
  if (env != NULL && *env != '\0') lenv = strlen(env)+1;
  else lenv = 0;

  len=strlen(fname)+1+lenv;

  if (lenv > 1 && *fname != '#') {	/* add default directory to file name */
    strncpy(tname,env,sizeof(tname)-1);
#ifdef UNIX
    strcat(tname,"/");
#endif
    }
  else tname[0]='\0';

  /* get to the end of the current list */
  while (cur_lib_p && cur_lib_p->next) {cur_lib_p = cur_lib_p->next;}

  /* add fname to tname, allocate space, and move to space */
  strncat(tname,fname,sizeof(tname)-strlen(tname)-1);
  len=strlen(tname)+1;
  if ((lbptr=calloc(len,sizeof(char)))==NULL) {
    fprintf(stderr,"no more space for filenames: %s ignored\n",fname);
    return cur_lib_p;
  }
  else {
    strncpy(lbptr,tname,len);
    lbptr[len-1]='\0';
    /* have a file name to add, I need a lib_struct */
    if ((this_lib_p = (struct lib_struct *)calloc(1,sizeof(struct lib_struct)))==NULL) {
      fprintf (stderr,"*** Error -- Cannot allocate lib_struct for %s\n",tname);
      return NULL;
    }
    else {
      this_lib_p->file_name = lbptr;
      if (cur_lib_p != NULL) {cur_lib_p->next = this_lib_p;}
      return this_lib_p;
    }
  }
}

char *
ulindex(char *str, char *chr)
{
  char c;
 
  c = tolower((int)(*chr));

  while (*str != '\0' && tolower(*str) !=c ) str++;
  if (*str=='\0') return NULL;
  else return str;
}
