#!/bin/sh
echo ""
echo "STARTING FASTA36" `date` "on" `hostname`
echo `uname -a`
echo ""
if [ ! -d results ]; then
  mkdir results
fi
echo "starting fasta36 - protein" `date`
../bin/fasta36 -XM2G -q -m 6 -Z 100000 ../seq/mgstm1.aa:1-100 q > results/test2G_m1.ok2.html
../bin/fasta36 -q -XM2G -S -z 11 -O results/test2G_m1.ok2_p25 -s P250 ../seq/mgstm1.aa:100-218 q
echo "done"
echo "starting fastxy36" `date`
../bin/fastx36 -q -XM2G -m 9c -S ../seq/mgtt2_x.seq q 1 > results/test2G_t2.xk1
../bin/fasty36 -q -XM2G -S ../seq/mgtt2_x.seq q > results/test2G_t2.yk2
../bin/fastx36 -q -XM2G -m 9c -S -z 2 ../seq/mgstm1.esq a > results/test2G_m1.xk2z2
../bin/fasty36 -q -XM2G -S -z 2 ../seq/mgstm1.esq a > results/test2G_m1.yk2z2
echo "done"
echo "starting fastxy36 rev" `date`
../bin/fastx36 -q -XM2G -m 9c -m 5 ../seq/mgstm1.rev q > results/test2G_m1.xk2r
../bin/fasty36 -q -XM2G -m 5 -M 200-300 -z 2 ../seq/mgstm1.rev q > results/test2G_m1.yk2rz2
../bin/fasty36 -q -XM2G -m 5 -z 11 ../seq/mgstm1.rev q > results/test2G_m1.yk2rz11
echo "done"
echo "starting ssearch36" `date`
../bin/ssearch36 -q -XM2G -m 9c -S -z 3 ../seq/mgstm1.aa  q > results/test2G_m1.ssz3
../bin/ssearch36 -q -XM2G -M 200-300 -z 2 -Z 100000 -s P250 ../seq/mgstm1.aa q > results/test2G_m1.ss_p25
echo "done"
if [ -e ../bin/ssearch36s ]; then
    echo "starting ssearch36s" `date`
    ../bin/ssearch36s -q -XM2G -m 9c -S -z 3 ../seq/mgstm1.aa  q > results/test2G_m1.sssz3
    ../bin/ssearch36s -q -XM2G -M 200-300 -z 2 -Z 100000 -s P250 ../seq/mgstm1.aa q > results/test2G_m1.sss_p25
    echo "done"
fi
echo "starting prss36(ssearch/fastx)" `date`
../bin/ssearch36 -q -XM2G -k 1000 -a ../seq/mgstm1.aa ../seq/xurt8c.aa  > results/test2G_m1.rss
../bin/fastx36 -q -XM2G -k 1000 ../seq/mgstm1.esq ../seq/xurt8c.aa > results/test2G_m1.rfx
echo "done"
echo "starting ggsearch36/glsearch36" `date`
../bin/ggsearch36 -q -XM2G -m 9i -w 80 ../seq/hahu.aa q > results/test2G_h1.gg
../bin/glsearch36 -q -XM2G -m 9i -w 80 ../seq/hahu.aa q > results/test2G_h1.gl
../bin/ggsearch36 -q -XM2G ../seq/gtt1_drome.aa q > results/test2G_t1.gg
../bin/glsearch36 -q -XM2G ../seq/gtt1_drome.aa q > results/test2G_t1.gl
echo "done"
echo "starting fasta36 - DNA" `date`
../bin/fasta36 -S -q -XM2G ../seq/mgstm1.nt %RMB 4 > results/test2G_m1.ok4
../bin/fasta36 -S -q -XM2G ../seq/mgstm1.rev %RMB 4 > results/test2G_m1.ok4r
echo "done"
#echo "starting tfasta36" `date`
#tfasta36 -q -XM2G ../seq/mgstm1.aa %RMB > results/test2G_m1.tk2
#echo "done"
echo "starting tfastxy36" `date`
../bin/tfastx36 -m 9c -q -XM2G -i -3 -m 6 ../seq/mgstm1.aa %p > results/test2G_m1.tx2.html
../bin/tfasty36 -q -XM2G -i -3 -N 5000 ../seq/mgstm1.aa %p > results/test2G_m1.ty2
echo "done"
echo "starting fastf36" `date`
../bin/fastf36 -q -XM2G ../seq/m1r.aa q > results/test2G_mf.ff
../bin/fastf36 -q -XM2G ../seq/m1r.aa q > results/test2G_mf.ff_s
echo "done"
echo "starting tfastf36" `date`
../bin/tfastf36 -q -XM2G ../seq/m1r.aa %r > results/test2G_mf.tfr
echo "done"
echo "starting fasts36" `date`
../bin/fasts36 -q -XM2G -V '*?@' ../seq/ngts.aa q > results/test2G_m1.fs1
../bin/fasts36 -q -XM2G ../seq/ngt.aa q > results/test2G_m1.fs
../bin/fasts36 -q -XM2G -n ../seq/mgstm1.nts m > results/test2G_m1.nfs
echo "starting fastm36" `date`
../bin/fastm36 -q -XM2G ../seq/ngts.aa q > results/test2G_m1.fm
../bin/fastm36 -q -XM2G -n ../seq/mgstm1.nts m > results/test2G_m1.nfm
echo "done"
echo "starting tfasts36" `date`
../bin/tfasts36 -q -XM2G ../seq/n0.aa %r > results/test2G_m1.ts_r
echo "starting lalign36" `date`
../bin/lalign36 -q -XM2G -k 1000 ../seq/mchu.aa ../seq/mchu.aa > results/test2G_mc.lal
../bin/lalign36 -q -XM2G -z 3 ../seq/mchu.aa ../seq/mchu.aa > results/test2G_mc.lal_z3
../bin/lalign36 -q -XM2G -s BL62 -f -11 -g -1  ../seq/mchu.aa ../seq/mchu.aa > results/test2G_mc.lal_bl62
../bin/lalign36 -q -XM2G -k 1000 ../seq/mwkw.aa ../seq/mwkw.aa > results/test2G_mw.lal
../bin/lalign36 -z 3 -q -XM2G ../seq/mwkw.aa ../seq/mwkw.aa > results/test2G_mw.lal_z3
../bin/lalign36 -q -XM2G -s BL62 -f -11 -g -1 ../seq/mwkw.aa ../seq/mwkw.aa > results/test2G_mw.lal_bl62
echo "FINISHED" `date`
