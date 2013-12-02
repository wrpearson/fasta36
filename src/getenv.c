#include <stdio.h>
#include <string.h>
#include <stdlib.h>

#define MAXENV 1024
char *envstr;

char *mgetenv(str)
char *str;
{
	static int EnvInit=0;

	char *eptr, *esptr, *bp;
	int i,esize;
	FILE *fenv;
 	
	if (EnvInit==0) {
		EnvInit=1;
		if ((fenv=fopen("environment","r"))!=NULL) {
			if ((envstr=malloc((size_t)(esize=MAXENV)))==NULL) {
				fclose(fenv); goto noenv;}
			esptr=envstr; esize -= 10;
			while (fgets(esptr,esize,fenv)!=NULL) {
				if ((bp=strchr(esptr,'\n'))!=NULL) *bp='\0';
				esize -= (i=strlen(esptr)+1);
				esptr += i;
				}
			fclose(fenv);
			esptr='\0';
			}
		else envstr=NULL;
	}
	
	if (envstr==NULL) return NULL;
	else {		
		for (eptr=envstr; *eptr; eptr += strlen(eptr)+1) {
			if (strncmp(str,eptr,(long)strlen(str))==0) {
				return strchr(eptr,'=')+1;
				}
			}
		return NULL;
		}
noenv:	envstr=NULL; return NULL;
	}

strnpcpy(to,from,max)
	char *to; Str255 from; size_t max;
{
	size_t i, n;
	
	n = (*from<max) ? *from : max;
	from++;

	for (i=0; i<n; i++) *to++ = *from++;
	if (n<max) *to='\0';
	}
