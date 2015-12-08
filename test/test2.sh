#!/bin/sh
echo ""
echo "starting FASTA36" `date` "on" `hostname`
echo `uname -a`
echo ""
echo "starting fasta36 - protein" `date`
if [ ! -d results ]; then
 mkdir results
fi
../bin/fasta36 -S -z 21 -s BP62 -c O ../seq/mgstm1.aa q > results/test2o_m1.ok2_bp62
../bin/fasta36 -S -z 21 -s BP62 ../seq/mgstm1.aa q > results/test2_m1.ok2_bp62
../bin/fasta36 -S -z 21 -c O ../seq/mgstm1.aa q > results/test2o_m1.ok2_z21
../bin/fasta36 -S -z 21 ../seq/mgstm1.aa q > results/test2_m1.ok2_z21
../bin/fasta36 -S -m BB  ../seq/mgstm1.aa q > results/test2_m1.ok2mB
../bin/fasta36 -S -m 8CC  ../seq/mgstm1.aa q > results/test2_m1.ok2m8CC
../bin/fasta36 -S -m 8CB ../seq/mgstm1.aa q > results/test2_m1.ok2m8CB
echo "done"
echo "starting fastxy36" `date`
../bin/fastx36 -m 9c -S -c O -q ../seq/mgtt2_x.seq q > results/test2o_t2.xk2m9c
../bin/fastx36 -m 8C -S -q ../seq/mgtt2_x.seq q > results/test2_t2.xk2m8C
../bin/fastx36 -m 8CB -S -q ../seq/mgtt2_x.seq q > results/test2_t2.xk2m8CB
../bin/fastx36 -m BB -S -q -H ../seq/mgtt2_x.seq q > results/test2_t2.xk2mB
../bin/fasty36 -S -c O -q ../seq/mgtt2_x.seq q > results/test2o_t2.yk2
../bin/fasty36 -S -c O -q ../seq/mgtt2_x.seq q > results/test2_t2.yk2
../bin/fasty36 -S -m8 -q ../seq/mgtt2_x.seq q > results/test2_t2.yk2m8
../bin/fastx36 -m 9c -c O  -S -q -z 22 ../seq/mgstm1.esq q > results/test2o_m1.xk2m9cz22
../bin/fastx36 -m 8C -S -q  ../seq/mgstm1.esq q > results/test2_m1.xk2m8Cz22
../bin/fastx36 -m 8CB -S -q ../seq/mgstm1.esq q > results/test2_m1.xk2m8CBz22
../bin/fasty36 -S -c O -q -z 21 ../seq/mgstm1.esq q > results/test2o_m1.yk2z21
../bin/fasty36 -S -q -z 21 ../seq/mgstm1.esq q > results/test2_m1.yk2z21
echo "done"
echo "starting ssearch36" `date`
../bin/ssearch36 -m 9c -S -z 22 -q ../seq/mgstm1.aa  q > results/test2_m1.ssm9cz22
../bin/ssearch36 -m 9C -S -z 21 -q ../seq/mgstm1.aa  q > results/test2_m1.ssm9Cz21
../bin/ssearch36 -m 8C -S -q ../seq/mgstm1.aa  q > results/test2_m1.ssm8C
../bin/ssearch36 -m 8CB -S -q ../seq/mgstm1.aa  q > results/test2_m1.ssm8CB
echo "done"
echo "starting fasta36 - DNA" `date`
../bin/fasta36 -S -q -c O -r "+2/-4" ../seq/mgstm1.nt %RMB 4 > results/test2o_m1.ok4z1r24
../bin/fasta36 -S -q -r "+2/-4" ../seq/mgstm1.nt %RMB 4 > results/test2_m1.ok4z1r24
../bin/fasta36 -S -q -c O ../seq/mgstm1.rev %RMB 4 > results/test2o_m1.ok4r
../bin/fasta36 -S -q  ../seq/mgstm1.rev %RMB 4 > results/test2_m1.ok4r
../bin/fasta36 -m BB -q  ../seq/mgstm1.rev %RMB > results/test2_m1.ok6mB
../bin/fasta36 -m 8C -q  ../seq/mgstm1.rev %RMB > results/test2_m1.ok6m8C
../bin/fasta36 -m 8CB -q  ../seq/mgstm1.rev %RMB > results/test2_m1.ok6m8CB
echo "done"
echo "starting tfastxy36" `date`
../bin/tfastx36 -c O -m 9c -q -i -3 ../seq/mgstm1.aa %p > results/test2o_m1.tx2
../bin/tfastx36 -m 8C -q -i -3 ../seq/mgstm1.aa %p > results/test2_m1.tx2_m8C
../bin/tfastx36 -m 8CB -q -i -3 ../seq/mgstm1.aa %p > results/test2_m1.tx2_m8CB
../bin/tfasty36 -c O  -q -i -3 -N 5000 ../seq/mgstm1.aa %p > results/test2o_m1.ty2
../bin/tfasty36 -q -i -3 -N 5000 ../seq/mgstm1.aa %p > results/test2_m1.ty2
echo "done" `date`
