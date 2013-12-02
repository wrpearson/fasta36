
#define MAX_STR	512 /* standard label/message buffer */

#ifndef XTERNAL
long pminx, pmaxx, pminy, pmaxy;
int max_x=540, max_y=540;
double fxscal, fyscal, fxoff, fyoff;

int *linarr;
int nlinarr=5;

char lvstr[MAX_STR];

double elinval[4]={1e-4,1e-2,1.0,100.0};
double blinval[4]={40.0,30.0,20.0,10.0};
int ilinval[4]={200,100,50,25};
extern int have_bits, have_zdb;
#else
int have_bits=0, have_zdb=0;
long zdb_size=1;
#endif

void openplt(long, long, int, int, char *, char *);
void closeplt();
void drawdiag(long n0, long n1);
void closepl();
void move(int, int);
void cont(int, int);
void clsline();
void xaxis(long, int, char *);
void yaxis(long, int, char *);
void legend();
void linetype(int);
void opnline(int s, double bits);
void newline();
void clsline();
void move(int, int);
void draw(int, int);
void sxy_move(int, int);
void sxy_draw(int, int);
void draw_str(char *);
void draw_sstr(char *);

double bit_to_E(double bits);
