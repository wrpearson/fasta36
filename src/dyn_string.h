/* $Id: dyn_string.h 1197 2013-07-19 20:25:19Z wrp $ */
/* $Revision: 1197 $  */

/* structure, functions for dynamic strings */

struct dyn_string_str {
  char *string;
  int c_size;
  int mx_size;
  int inc;
};

/* initial allocation */
struct dyn_string_str *init_dyn_string(int size, int inc);

/* strcpy */
void dyn_strcpy(struct dyn_string_str *dyn_string, char *value);

/* strcat */
void dyn_strcat(struct dyn_string_str *dyn_string, char *value);

/* free */
void free_dyn_string(struct dyn_string_str *dyn_string);

/* reset */
void reset_dyn_string(struct dyn_string_str *dyn_string);

/* initialize to '\0' */
#define NULL_dyn_string(str) str->string[0]='\0'; str->c_size = 0;

