
## The FASTA package - protein and DNA sequence similarity searching and alignment programs

Changes in **fasta-36.3.8g** released 23-Oct-2018

  1. (Oct. 2018) Improvements to scripts in the `psisearch2/` directory:

    1. `psisearch2/m89_btop_msa2.pl`

      1. the `--clustal` option produces a "CLUSTALW (1.8)", which is required for some downstream programs

      2. the `--trunc_acc` option removes the database and accession from identifiers of the form:
         `sp|P09488|GSTM1_HUMAN` to produce `GSTM1_HUMAN`.

      3. the `--min_align` option specifies the fraction of the query sequence that must be aligned 
         (q_end-q_start+1)/q_length)

      Together, these changes make it possible for the output of `m89_btop_msa2.pl` to be used by
the EMBOSS program `fprotdist`.

    2. A more general implementation of `psisearch2_msa_iter.sh`, which does `psisearch2` one iteration at a time, and a new equivalent `psisearch2_msa_iter_bl.sh`, which uses `psiblast` to do the search.

  2. (Oct. 2018) A small restructuring of the `make/Makefiles` to remove the `-lz` dependence for non-debugging scripts (and add it back when -DDEBUG is used).

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

