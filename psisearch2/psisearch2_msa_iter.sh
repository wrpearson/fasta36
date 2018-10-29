#!/bin/sh

################
# example that runs psisearch2_msa.pl iteratively through 5 iterations.
# Equivalent to:
# psisearch2_msa.pl --query CL0238_emb.fa --num_iter 5 --db /slib2/fa_dbs/rpd3_pfam28_lib.lseg
#


PS_BIN=~/Devel/fa36_v3.8/psisearch2
Q_DIR="../seq"
FA_DB=/slib2/fa_dbs/qfo78.lseg
BL_DB=/slib2/bl_dbs/qfo78
DB=$FA_DB

OUT_SUFF='qm8CB'

M_FORMAT='m8CB'
ITERS='2 3 4 5'

for q_file_p in $*; do

    q_file=${q_file_p##*/}
    echo $q_file

    # iteration 1:

    $PS_BIN/psisearch2_msa.pl --query $Q_DIR/$q_file  --num_iter 1 --db $DB --int_mask query --end_mask query --out_suffix $OUT_SUFF --m_format $M_FORMAT

    # iteration 2 - 5
    for it in $ITERS; do
	prev=$(($it-1))
	$PS_BIN/psisearch2_msa.pl --query $Q_DIR/$q_file --num_iter 1 --db $DB --int_mask query --end_mask query --out_suffix $OUT_SUFF --this_iter $it --prev_m89res $q_file.it${prev}.$OUT_SUFF --m_format $M_FORMAT
    done

done
