
## The FASTA package - protein and DNA sequence similarity searching and alignment programs

The **FASTA** (pronounced FAST-Aye, not FAST-Ah) programs are a
comprehensive set of similarity searching and alignment programs for
searching protein and DNA sequence databases.  Like the **BLAST** programs `blastp` and `blastn`, the `fasta` program itself uses a rapid heuristic strategy for finding similar regions in protein and DNA sequences.  But in
addition to heuristic similarity searching, the FASTA package provides
programs for rigorous local (`ssearch`) and global (`ggsearch`)
similarity searching, as well as a program for finding non-overlapping
sequence similarities (`lalign`).  Like BLAST, the FASTA package also
includes programs for aligning translated DNA sequences against
proteins (`fastx`, `fasty` are equivalent to `blastx`, `tfastx`,
`tfasty` are similar to `tblastn`).

####September, 2016

The current FASTA version is fasta-36.3.8e.

The fasta-36.3.6e version includes a new directory, `psisearch2`, with
scripts to run iterative PSSM (PSI-BLAST or SSEARCH36) searches using
an improved strategy for reducing PSSM contamination due to alignment
over-extension.

As of November, 2014, the FASTA program code is available under the
Apache 2.0 open source license.

Up-to-date release notes are available in the file `doc/readme.v36`.

Documentation on the FASTA programs is available in the files:

dir/file | description
----------|------------
`doc/fasta36.1` | (unix man page)
`doc/changes_v36.html` | (short descriptions of enhancements to FASTA programs)
`doc/readme.v36` | (text descriptions of bug fixes and version history)
`doc/fasta_guide.tex` | (Latex file which describes fasta36, and provides an introduction to the FASTA programs, their use and installation.)
`doc/fasta_guide.pdf1` | (printable/viewable description of fasta-36)

`fasta_guide.pdf` provides background information on installing the
fasta programs (in particular, the `FASTLIBS` file), that new users of
the fasta3 package may find useful.

Parts of the FASTA package are distributed across several sub-directories

dir | description
----|------------
`bin/` | (pre-compiled binaries for some architectures)
`conf/` | example `FASTLIBS` files (files for finding libraries)
`data/` | scoring matrices
`doc/` | documentation files
`make/` | make files
`misc/` | perl scripts to reformat -m 9 output, convert -R search.res files for 'R', and embed domains in shuffled sequences
`psisearch2/` | perl/python scripts implementing the new `psisearch2_msa` iterative PSSM search
`scripts/` |  perl scripts for -V (annotate alignments) and -E (expand library) options
`seq/` | test sequences
`src/` | source code
`sql/` | sql files and scripts for using the sql database access
`test/` | test scripts

For some binary distributions, only the `doc/`, `data/`, `seq/`, and `bin/`,
directories are provided.

To make the standard FASTA programs:
```
   cd src
   make -f ../make/Makefile.linux_sse2 all
```
where `../make/Makefile.linux_sse2` is the appropriate Makefile for your system. 

The executable programs will then be found in `../bin`
(e.g. `../bin/fasta36`, etc.)

For a simple test of a program, try (from the src directory)
```
   ../bin/fasta36 -q ../seq/mgstm1.aa ../seq/prot_test.lseg
```

