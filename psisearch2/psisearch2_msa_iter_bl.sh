#!/bin/sh

################
# example that runs psisearch2_msa.pl iteratively through 5 iterations using psiblast instead of ssearch
# Equivalent to:
# psisearch2_msa.pl --pgm psiblast --query query.aa --num_iter 5 --db /slib2/bl_dbs/qfo78
#

PS_BIN=~/Devel/fa36_v3.8/psisearch2
q_file=$1
m_format='m8CB'
SRC_QDIR=../hum_1dom200_queries

iters='2 3 4 5'
# iters=''

for q_file_p in $*; do

    q_file=${q_file_p##*/}
    echo $q_file

    # iteration 1:
#    echo "$PS_BIN/psisearch2_msa.pl --pgm psiblast --query $SRC_QDIR/$q_file --num_iter 1 --db /slib2/bl_dbs/qfo78 --int_mask query --end_mask query --out_suffix q_pblt --m_format $m_format --save_list asnbin"
    $PS_BIN/psisearch2_msa.pl --pgm psiblast --query $SRC_QDIR/$q_file --num_iter 1 --db /slib2/bl_dbs/qfo78 --int_mask query --end_mask query --out_suffix q_pblt --m_format $m_format --save_list asntxt

    # iteration 2 - 5
    for it in $iters; do
	prev=$(($it-1))
	$PS_BIN/psisearch2_msa.pl --pgm psiblast --query $SRC_QDIR/$q_file --num_iter 1 --db /slib2/bl_dbs/qfo78 --int_mask query --end_mask query --out_suffix q_pblt --this_iter $it --prev_m89res $q_file.it${prev}.q_pblt --m_format $m_format --save_list asntxt
    done
done
