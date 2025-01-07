#!/bin/sh
echo ""
echo "STARTING FASTA36" `date` "on" `hostname`
echo `uname -a`
echo ""
if [ ! -d results ]; then
  mkdir results
fi

export FA_DB=$SLIB2/fa_dbs/qfo20.lseg

echo "starting fasta36 - protein" `date`
../bin/fasta36 -q -mBB -S -sBP62 ../seq/mgstm1.aa $FA_DB > results/test_m1.ok2_BB
echo "done"
echo "starting fastxy36" `date`
../bin/fastx36 -mBB -S -q -sBP62 ../seq/mgtt2_x.seq $FA_DB > results/test_t2.xk2_BB
../bin/fasty36 -mBB -S -q -sBP62 ../seq/mgtt2_x.seq $FA_DB > results/test_t2.xk2_BB
echo "done"
echo "starting ssearch36" `date`
../bin/ssearch36 -mBB -S -q -sBP62 ../seq/mgstm1.aa  $FA_DB > results/test_m1.ss_BB
echo "done"
echo "starting fasta36 - DNA" `date`
../bin/fasta36 -q -mBB ../seq/mgstm1.nt %RMB 4 > results/test_m1.ok4_BB
../bin/fasta36 -q -mBB ../seq/mgstm1.rev %RMB 4 > results/test_m1.ok4r_BB
echo "done"
#echo "starting tfasta36" `date`
#tfasta36 -q ../seq/mgstm1.aa %RMB > results/test_m1.tk2
#echo "done"
echo "starting tfastxy36" `date`
../bin/tfastx36 -q -mBB -i -3 -N 5000 -sBP62 ../seq/mgstm1.aa %p > results/test_m1.tx2_BB
../bin/tfasty36 -q -mBB -i -3 -N 5000 -sBP62 ../seq/mgstm1.aa %p > results/test_m1.ty2_BB
echo "done"
echo "FINISHED" `date`
