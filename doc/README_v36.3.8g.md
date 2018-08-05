

## The FASTA package - protein and DNA sequence similarity searching and alignment programs

Changes in **fasta-36.3.8g** released 5-Aug-2018

1. (Apr 2018) incorporation of "-t t" temrination codes ("*") in -m 8CB, -m 8CC, and -m9C so that aligned termination codons are indicated as "**" (-m8CB) or
"*1" (-m8CC, -m9C).

2. (Mar 2018) Updates to scripts/annot_blast_btop2.pl to provide
subalignment scoring for blastp searches (BLOSUM62 only).  (see
doc/readme.v36)

3. (Feb. 2018) a new extended option, -XB, which causes percent
identity, percent similarity, and alignment length to be calculated
using the BLAST model, which does not count gaps in the alignment
length.

see readme.v36 for other bug fixes.

Changes in **fasta-36.3.8g** released 31-Dec-2017

1. (December, 2017) -- Make statistical thresholds more robust for
   small E()-values with normally distributed scores (ggsearch36,
   glsearch36).

2. (September, 2017) Treat all lower-case queries as uppercase with -S option.

3. (May, 2017) Improvements/fixes to sub-alignment scoring strategies.

4. Improvements/fixes to psisearch2 scripts.

For more detailed information, see `doc/readme.v36`.

