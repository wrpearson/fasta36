#!/bin/csh -f
echo "starting fasta36 - protein" `date`
foreach z ( 1 2 3 11 12 21 22)
../bin/fasta36 -q  -z $z -d 0 ../seq/mgstm1.aa a > results/test_m1_a.ok2_${z}
../bin/fasta36 -q  -z $z -d 0 ../seq/oohu.aa a > results/test_m1_b.ok2_${z}
../bin/fasta36 -q -S -z $z -d 0 ../seq/prio_atepa.aa a > results/test_m1_c.ok2S_${z}
../bin/fasta36 -q -S -z $z -d 0 ../seq/h10_human.aa a > results/test_m1_d.ok2S_${z}
../bin/fasta36 -c -1 -q  -z $z -d 0 ../seq/mgstm1.aa a > results/test_m1_a.c1_ok2_${z}
../bin/fasta36 -c -1 -q  -z $z -d 0 ../seq/oohu.aa a > results/test_m1_b.c1_ok2_${z}
../bin/fasta36 -c -1 -q -S -z $z -d 0 ../seq/prio_atepa.aa a > results/test_m1_c.c1_ok2S_${z}
../bin/fasta36 -c -1 -q -S -z $z -d 0 ../seq/h10_human.aa a > results/test_m1_d.c1_ok2S_${z}
end
echo "done"
echo "starting ssearch36" `date`
foreach z ( 1 2 3 11 21 22)
../bin/ssearch36 -q  -z $z -d 0 ../seq/mgstm1.aa a > results/test_m1_a.ssS_${z}
../bin/ssearch36 -q  -z $z -d 0 ../seq/oohu.aa a > results/test_m1_b.ssS_${z}
../bin/ssearch36 -q -sBL62 -d 0 -S -f -11 -z $z ../seq/prio_atepa.aa a > results/test_m1_c.ssSbl62_${z}
../bin/ssearch36 -q -sBL62 -d 0 -S -f -11 -z $z ../seq/h10_human.aa a > results/test_m1_d.ssSbl62_${z}
end
echo "done"
echo "done" `date`
