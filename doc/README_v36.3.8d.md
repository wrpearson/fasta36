

## The FASTA package - protein and DNA sequence similarity searching and alignment programs

Changes in **fasta-36.3.8d** released 13-April-2016:

1. Various bug fixes to `pssm_asn_subs.c` that avoid coredumps when
   reading NCBI PSSM ASN.1 binary files.  `pssm_asn_subs.c` can now read
   UUPACAA sequences.

2. default gap penalties for VT40 (from -14/-2 to -13/-1), VT80 (from
   -14/-2 to -11/-1), and VT120 (from -10/-1 to 11/-1) have changed
   slightly.

3. Introduction of `scripts/m9B_btop_msa.pl` and
   `scripts/m8_btop_msa.pl`, which uses the BTOP (`-m 9B` or `-m 8CB`)
   encoded alignment strings to produce a query driving multiple
   sequence alignment (MSA) in ClustalW format.  This MSA can be used
   as input to `psiblast` to produce an ASN.1 PSSM.

4. The `scripts/annot_blast_btop2.pl` script replaces
   `scripts/annot_blast_btop.pl` and allows annotation of both the query
   and subject sequences.

5. Various domain annotation scripts have been renamed for clarity.
   For example, `ann_feats_up_sql.pl` uses an SQL implementation of
   Uniprot features tables to annotate domains.  Likewise,
   `ann_pfam_www.pl` gets domain information from the Pfam web site,
   while `ann_pfam27.pl` gets the information from the downloaded
   Pfam27 mySQL tables, and `ann_pfam28.pl` uses the Pfam28 mySQL
   tables.

6. percent identity in sub-alignment scores is calculated like a BLAST
   percent identity -- gaps are not included in the denominator.

For more detailed information, see `doc/readme.v36`.

