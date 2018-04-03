#!/bin/sh

################
# example that runs psisearch2_msa.pl iteratively through 5 iterations.
# Equivalent to:
# psisearch2_msa.pl --query CL0238_emb.fa --num_iter 5 --db /slib2/fa_dbs/rpd3_pfam28_lib.lseg
#

# iteration 1:

psisearch2_msa.pl --query CL0238_emb.fa --num_iter 1 --db /slib2/fa_dbs/rpd3_pfam28_lib.lseg --int_mask query --end_mask query --out_suffix qm9hi --m_format m9

# iteration 2 - 5
for it in 2 3 4 5; do
  prev=$(($it-1))
  psisearch2_msa.pl --query CL0238_emb.fa --num_iter 1 --db /slib2/fa_dbs/rpd3_pfam28_lib.lseg --int_mask query --end_mask query --out_suffix qm9hi --this_iter $it --prev_m89res CL0238_emb.fa.it${prev}.qmi --m_format m9
done
