
## The FASTA package - protein and DNA sequence similarity searching and alignment programs

This directory contains the source code for the FASTA package of
programs (W. R. Pearson and D. J. Lipman (1988), "Improved Tools
for Biological Sequence Analysis", *PNAS 85:2444-2448*).  The current verion of the program is `fasta-36.3.8i`.

If you are reading this at
[fasta.bioch.virginia.edu/wrpearson/fasta/fasta36](https://fasta.bioch.virginia.edu/wrpearson/fasta/fasta36),
links are available to executable binaries for Linux, MacOS, and
Windows.  The source code is also available from
[github.com/wrpearson/fasta36](https://github.com/wrpearson/fasta36).

The FASTA package offers many of the same programs as `BLAST`, but
takes a different approach to statistical estimates, and provides
additional optimal programs for local (`ssearch36`) and global
(`ggsearch36`, `glsearch36`) alignment, and for non-overlapping
internal local alignments (`lalign36`).

The programs available include:

<table>
<tr>
<tr><td colspan=3><hr/></td><tr>
<th>FASTA </th><th> BLAST </th><th> description </th></tr>
<tr><td colspan=3><hr/></td><tr>
<tr>
<td> fasta36 </td><td> blastp/blastn </td><td> Protein and DNA local similarity search </td></tr>
<tr>
<td> ssearch36 </td><td> </td><td> optimal Smith-Waterman search -- vectorized on Intel and Arm architectures </td></tr>
<tr>
<td> ggsearch36 </td><td> </td><td> optimal global Needleman-Wunsche search -- vectorized on Intel and Arm architectures </td></tr>
<tr>
<td> glsearch36 </td><td> </td><td> optimal global(query)/local (library) search -- vectorized on Intel and Arm architectures </td></tr>
<tr>
<td> fastx36 / fasty36 </td><td> blastx </td><td> DNA query search against protein sequence database. (fasty36 uses a slower, more sophisticated frame shift aligner) </td></tr>
<tr>
<td> tfastx36 / tfasty36</td><td> tblastn </td><td> protein query search against DNA database</td></tr>
<tr><td colspan=3><hr/></td><tr>
<tr>
<td> fastf36 / tfastf36 </td><td> </td><td> compares an ordered peptide mixture against a protein (fastf36) or DNA (tfastf36) database </td></tr>
<tr>
<td> fastm36 / tfastm36 </td><td> </td><td> compares a set of ordered peptide against a protein (fastf36) or DNA (tfastf36) database or oligonucleotides against a DNA database</td></tr>
<tr>
<td> fasts36 / tfasts36 </td><td> </td><td> compares an unordered set of peptides against a protein (fasts36) or DNA (tfasts36) database </td></tr>
<tr><td colspan=3><hr/></td><tr>
<tr>
<td> lalign36 </td><td> </td><td> look for non-overlapping internal alignments, similar to a "dot-plot," but with statistical signficance </td></tr>
<tr><td colspan=3><hr/></td><tr>
</table>

Changes in **fasta-36.3.8i** May, 2023

1. restore the default `-s BL62` gap penalties to -8, -1 (they were -11, -1, matching `-s BP62`

2. restore functionality of `-A` option, which forces Smith-Waterman final display alignments with DNA (normally banded Smith-Waterman is used)

3. add `--id` option to `scripts/get_protein.py` to add a custom identifier

Changes in **fasta-36.3.8i** Nov, 2022

1. bug fix to remove duplicate variant annotations

2. update to scripts/get_protein.py and annotation scripts.

3. modify code to reduce mktemp compilation warning messages

4. changes to annotation scripts for Pfam shutdown; new ann_pfam_www.py, ann_pfam_sql.py

Changes in **fasta-36.3.8i** Sept, 2021

1. Enable translation table -t 9 for Echinoderms.  This bug has existed
   since alternate translation tables were first made available.

Changes in **fasta-36.3.8i** May, 2021

1. Add an option, -Xg, that preserves the gi|12345 string the score
   summary and alignment output.

Changes in **fasta-36.3.8i** Nov, 2020

1. fasta-36.3.8i (November, 2020) incorporates the SIMDe
   (SIMD-everywhere,
   https://github.com/simd-everywhere/simde/blob/master/simde/x86/sse2.h)
   macro definitions that allow the smith\_waterman\_sse2.c,
   global\_sse2.c, and glocal\_sse2.c code to be compiled on non-Intel
   architectures (currently tested on ARM/NEON).  Many thanks to
   Michael R. Crusoe (https://orcid.org/0000-0002-2961-9670) for the
   SIMDE code converstion, and to Evan Nemerson for creating SIMDe.

2. The code to read FASTA format sequence files now ignores lines with
   '#' at the beginning, for compatibility with PSI Extended FASTA
   Format (PEFF) files (http://www.psidev.info/peff).

Changes in **fasta-36.3.8h** May, 2020

1. fasta-36.3.8h (May 2020) fixes a bug that appeared when
multiple query sequences were searched against a large library
that would not fit in memory. In that case, the number of
library sequences and residues increased by the library size
with each new search.

2. More consistent formats for **ERROR** and **Warning** messages.

3. Corrections to code to address compiler warnings with gcc8/9.

4. addition of 's' option to show similarity in -m8CBls (or -m8CBs, -m8CBsl) and 'd' option to show raw (unaligned) domain information.

Changes in **fasta-36.3.8h** February, 2020

1. The license for Michael Farrar's Smith-Waterman sse2 code and global/glocal sse2 code is now open source (BSD), see COPYRIGHT.sse2 for details.

Changes in **fasta-36.3.8h** August, 2019

1. Modifications to support makeblastdb format v5 databases. Currently, only simple database reads have been tested.

Changes in **fasta-36.3.8h** March, 2019

1. Translation table 1 (`-t 1`) now translates 'TGA'->'U' (selenocysteine).

2. New script for extracting DNA sequences from genomes (`scripts/get_genome_seq.py`).  Currently works with human  (hg38), mouse (mm10), and rat (rn6).

Changes in **fasta-36.3.8h** January, 2019

1. Bug fixes: `fastx`/`tfastx` searches done with the `-t t` option  (which adds a `*` to protein sequences so that termination codons can  be matched), did not work properly with the `VT` series of matrices,  particularly `VT10`.  This has been fixed.

2. New features: Both query and library/subject sequences can be generated by specifying a program script, either by putting a `!` at the start of the query/subject file name, or by specifying library type `9`. Thus, `fasta36 \\!../scripts/get_protein.py+P09488+P30711 /seqlib/swissprot.fa` or `fasta36 "../scripts/get_protein.py+P09488+P30711 9" /seqlib/swissprot.fa` will compare two query sequences, `P09488` and `P30711`, to SwissProt, by downloading them from Uniprot using the `get_protein.py` script (which can download sequences using either Uniprot or RefSeq protein accessions). Often, the leading `!` must be escaped from shell interpretation with `\\!`.

New scripts that return FASTA sequences using accessions or genome coordinates are available in `scripts/`. `get_protein.py`, `get_uniprot.py`, `get_up_prot_iso_sql.py` and `get_refseq.py`. `get_refseq.py` can download either protein or mRNA RefSeq entries. `get_up_prot_iso_sql.py` retrieves a protein and its isoforms from a MySQL database.

`get_genome_seq.py` extracts genome sequences using coordinates from local reference genomes (`hg38` and `mm10` included by default).

Changes in **fasta-36.3.8h** December, 2018

The `scripts/ann_exons_up_www.pl` and `ann_exons_up_sql.pl` now include the option `--gen_coord` which provides the associated genome coordinate (including chromosome) as a feature, indicated by `'<'` (start of exon) and `'>'` (end of exon).

Changes in **fasta-36.3.8h** released November, 2018

**fasta-36.3.8h** provides new scripts and modifications to the   `fasta` programs that normalize the process of merging sub-alignment   scores and region information into both FASTA and BLAST results.  To   move BLASTP towards FASTA with respect to alignment annotation and   sub-alignment scoring:

1. The `blastp_annot_cmd.sh` runs a blast search, finds and scores   domain information for the alignments, and merges this information   back into the blast output `.html` file.  This script uses: 

   1. `annot_blast_btab2.pl --query query.file --ann_script annot_script.pl --q_ann_script annot_script.pl blast.btab_file > blast.btab_file_ann` (a blast tabular file with one or two new fields, an annotation field and (optionally with --dom_info) a raw domain content field.
   2. `merge_blast_btab.pl --btab blast.btab_file_ann blast.html > blast_ann.html`  (merge the annotations and domain content information in the `blast.btab_file_ann` file together with the standard blast output file to produce annotated alignments.
   3. In addition, `rename_exons.py` is available to rename exons (later other domains) in the subject sequences to match the exon labeling in the aligned query sequence.
   4. `relabel_domains.py` can be used to adjust color sets for homologous domains.

2.  There is also an equivalent `fasta_annot_cmd.sh` script that provides similar funtionality for the FASTA programs.  This script does not need to use `annot_blast_btab2.pl` to produce domain subalignment scores (that functionality is provided in FASTA), but it also can use `merge_fasta_btab.pl` and `rename_exons.py` to modify the names of the aligned exons/domains in the subject sequences.

3. To support the independence of the `blastp`/`fasta` output from html annotation, the FASTA package includes some new options:

   1. The `-m 8CBL` option includes query sequence length and subject sequence length in the blast tabular output.  In addition, if domain annotations are available, the raw domain coordinates are provided in an additional field after the annotation/subalignment scoring field.  `-m 8CBl` provides the sequence lengths, but does not add the raw domain coordinates.

   2. The `-Xa` option prevents annotation information from being included in the html output -- it is only available in the `-m 8CB`  (or `-m 8CBL/l`) output

   3. To reduce problems with spaces in script arguements, annotation scripts with spaces separating arguments can use '+' instead of ' '.

   4. The `fasta_annot_cmd.sh` script produces both a conventional alignment on `stdout` and a `-m 8CBL` alignment, which is sent to a separate file, which is separated from the `-m F8CBL` option with a `=`, thus `-m F8CBL=tmp_output.blast_tab`.

Changes in **fasta-36.3.8g** released 23-Oct-2018

1. (Oct. 2018) Improvements to scripts in the `psisearch2/` directory:

   1. `psisearch2/m89_btop_msa2.pl`
      1. the `--clustal` option produces a "CLUSTALW (1.8)", which is required for some downstream programs
      2. the `--trunc_acc` option removes the database and accession from identifiers of the form: `sp|P09488|GSTM1_HUMAN` to produce `GSTM1_HUMAN`.
      3. the `--min_align` option specifies the fraction of the query sequence that must be aligned `(q_end-q_start+1)/q_length)`
   Together, these changes make it possible for the output of `m89_btop_msa2.pl` to be used by the EMBOSS program `fprotdist`.

   2. A more general implementation of `psisearch2_msa_iter.sh`, which does `psisearch2` one iteration at a time, and a new equivalent `psisearch2_msa_iter_bl.sh`, which uses `psiblast` to do the search.

* (Oct. 2018) A small restructuring of the `make/Makefiles` to remove the `-lz` dependence for non-debugging scripts (and add it back when -DDEBUG is used).

Changes in **fasta-36.3.8g** released 5-Aug-2018

1. (Apr 2018) incorporation of `-t t1` termination codes ("*") in `-m 8CB`, `-m 8CC`, and `-m9C` so that aligned termination codons are indicated as `**` (`-m8CB`) or `*1` (`-m8CC`, `-m9C`).

2. (Mar 2018) Updates to scripts/annot_blast_btop2.pl to provide subalignment scoring for blastp searches (BLOSUM62 only).  (see doc/readme.v36)

3. (Feb. 2018) a new extended option, `-XB`, which causes percent identity, percent similarity, and alignment length to be calculated using the BLAST model, which does not count gaps in the alignment length.

see readme.v36 for other bug fixes.

Changes in **fasta-36.3.8g** released 31-Dec-2017

1. (December, 2017) -- Make statistical thresholds more robust for small E()-values with normally distributed scores (`ggsearch36`,`glsearch36`).

2. (September, 2017) Treat lower-case queries with no upper-case residues as uppercase with `-S` option.

3. (May, 2017) Improvements/fixes to sub-alignment scoring strategies.

4. Improvements/fixes to psisearch2 scripts.

For more detailed information, see `doc/readme.v36`.

