

/*****************************************************************************/
/***  (pmerge.c)                                                           ***/
/***                                                                       ***/
/*** merge fasta-formatted protein sequence databases                      ***/
/*** that have been differently masked                                     ***/
/*** into a single 'sequence-or'-ed output database                        ***/
/***                                                                       ***/
/*****************************************************************************/


#include "pgenwin.h"

/* forward declarations */
void printseq(struct Sequence *, FILE *);


int main(int argc, char *argv[]) {

   struct Database **dbs;
   struct Sequence *seq, *newseq;
   int i, ndbs;
   int j, len;

   if (argc<=1)
     {
      fprintf(stderr, "\n\
merge differently masked protein sequence databases\n\n\
 Usage: nmerge <db1> ... <dbn>\n\
 (output to stdout)\n\n\
");
      exit(1);
     }

   genwininit();
   ndbs = argc - 1;
   dbs = (struct Database **) malloc(ndbs*(sizeof(struct Database *)));
   for (i=0; i<ndbs; i++)
     {
      dbs[i] = opendbase(argv[i+1]);
      if (dbs[i]==NULL)
        {
         fprintf(stderr, "Error opening file %s\n", argv[i+1]);
         exit(1);
        }
     }

   for (newseq=firstseq(dbs[0]); newseq!=NULL; newseq=nextseq(dbs[0]))
     {
      len = newseq->length;
      for (i=1; i<ndbs; i++)
        {
         seq = nextseq(dbs[i]);
         if (seq==NULL)
           {
            fprintf(stderr, "premature end of file %s\n", argv[i+1]);
            fprintf(stderr, "%s\n\n", newseq->header);
            exit(1);
           }
         if (newseq->length!=seq->length)
           {
            fprintf(stderr, "sequence length discrepancy, file %s\n", argv[i+1]);
            fprintf(stderr, "%s\n%s\n\n", newseq->header, seq->header);
            exit(1);
           }
         for (j=0; j<len; j++)
           {
            if (toupper(newseq->seq[j])=='X') continue;
            if (toupper(seq->seq[j])=='X') newseq->seq[j] = 'x';
            else if (toupper(newseq->seq[j])!='X' &&
                     toupper(newseq->seq[j])!=toupper(seq->seq[j]))
              {
               fprintf(stderr, "sequence discrepancy, file %s, position %d\n", argv[i+1], j+1);
               fprintf(stderr, "%s\n%s\n\n", newseq->header, seq->header);
               exit(1);
              }
           }
         closeseq(seq);
        }
      printseq(newseq, stdout);
      closeseq(newseq);
     }

   for (i=0; i<ndbs; i++) closedbase(dbs[i]);
   exit(0);
  }

/*----------------------------------------------------------------(printseq)---*/

void printseq(struct Sequence *seq, FILE *fpout) {

   int i, len;
   int line, linelen;

   linelen = 60;

   fprintf(fpout, "%s\n", seq->header);
   for (i=0, line=0; i<seq->length; i++, line++)
     {
      if (line==linelen)
        {
         fputc('\n', fpout);
         fflush(fpout);
         line = 0;
        }
      fputc(seq->seq[i], fpout);
     }
   fputc('\n', fpout);

   return;
  }

/*-----------------------------------------------------------------------------*/
