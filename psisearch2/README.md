
## PSISEARCH2 - iterative PSSM-based similarity searching using PSIBLAST or SSEARCH36

#### September, 2016

`psisearch2_msa.pl` and `psisearch2_msa.py` (both scripts have
identical arguments and functionality) perform iterative searches
using PsiBLAST or ssearch36, but with additional options that
dramatically improve search selectivity.  In tests with challenging
queries, `psisearch2_msa.pl/py` searches often reduce the number of
false-positives more than ten fold, and sometimes 100-fold or more.

For a simple test of `psisearch2,` try (from the `psisearch2/` directory):

```
  ./psisearch2_msa.pl --query ../seq/mgstm1.aa --db ../seq/prot_test.lseg
```

This command should produce the output:
```
#./psisearch2_msa.pl --query ../seq/mgstm1.aa --db ../seq/prot_test.lseg
./psisearch2_msa.pl ssearch ../seq/mgstm1.aa ../seq/prot_test.lseg converged (2 iterations)
```
as well as four files:
```
mgstm1.aa.it1
mgstm1.aa.it1.bnd_out
mgstm1.aa.it2
mgstm1.aa.it2.bnd_out
```

Real iterative searches must be run against comprehensive sequence
databases, like SwissProt, RefSeq proteins, or Uniprot, e.g.:

```
   ./psisearch2_msa.pl --query ../seq/mgstm1.aa --db /slib2/swissprot.lseg
```

## More selective searches

By default, `psisearch2_msa.pl` simple runs a search program
(`ssearch36` by default, use `--pgm psisblast` to run `psiblast`),
scans the output to produce a multiple sequence alignment, which is
then used to build a `PSSM` for the next iteration.  Running
`psisearch2_msa.pl` for five iterations should produce results very
similar to running `psiblast` for five iterations.

`psisearch2_msa.pl` can perform much more selective searches using the
`--query_seed` option, which is equivalent to the `--int_seed=query`
and `--end_seed=query` options.  The `--query_seed` option causes the
`m89_btop_msa2.pl` program to insert query residues into the gapped
positions of subject sequences in the sequence library used to produce
the `PSSM`.  This `PSSM` is slightly less sensitive than the normal
model, but it is much less likely to produce alignment-overextension,
so it is much less likely for alignments to extend into neighboring,
non-homologous regions and contaminate the `PSSM` model.

In addition to `--query_seed`, two other options: `--align`, and
`--domain` can also be used to reduce alignment overextension.
`--align` causes `psisearch2_msa.pl` to include only portion of a
subject sequence that aligned the first time it shares significant
similarity to the query PSSM; additional residues from the sequence
that are aligned in later iterations are not included. This option can
be used with both `--pgm ssearch` and `--pgm psiblast`.

The `--domain` option uses a more sophisticated strategy for including
additional residues in a PSSM.  They are included only if the
similarity score across the domain has a probability less than 0.001
(q-value 30). This option is only available for `--pgm ssearch`, and
requires that a second option, `--annot_db`, be specified.  Typically
`--annot_db=pfam`.

## Customizing psisearch2

`psisearch2_msa.pl` is a script that uses other programs for
similarity searching and constructing the PSSM (and for annotating
domains in alignments if the `--domain` option is used).  The location
of these programs is defined in the `psisearch2_msa.pl` and
`psisearch2_msa.py` scripts using the `$pgm_bin` and `$pgm_data`
variable (perl, `pgm_bin` and `pgm_data` for python).  You will
probably need to modify those variables for your installation.  In
particular, the NCBI `datatool` program is required for producing the
asn binary files required by `ssearch36`.

## `psisearch2_msa.pl/py` output

In this current version, the `psisearch2_msa.pl` and
`psisearch2_msa.py` programs *ONLY* produce tab-delimited BTOP output.
The programs do not produce the traditional `psiblast` alignment,
which includes a list of hits with E()-values, and a set of
alignments.  More alignment output flexibility will be available soon.

