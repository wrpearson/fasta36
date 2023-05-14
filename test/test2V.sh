#!/bin/sh
echo ""
echo "starting FASTA36" `date` "on" `hostname`
echo `uname -a`
echo ""
echo "starting fasta36 - protein" `date`

FA_DB=$SLIB2/fa_dbs/qfo20.lseg

if [ ! -d results ]; then
 mkdir results
fi
../bin/fasta36 -V q\!../scripts/ann_feats_up_www2.pl -V \!../scripts/ann_feats_up_www2.pl -S -z 21 -s BP62 ../seq/gstm1_human.vaa $FA_DB > results/test2V_m1.ok2_bp62
../bin/fasta36 -V q\!../scripts/ann_feats_up_www2.pl -V \!../scripts/ann_feats_up_www2.pl -S -z 21 ../seq/gstm1_human.vaa $FA_DB > results/test2V_m1.ok2_z21
../bin/fasta36 -V q\!../scripts/ann_feats_up_www2.pl -V \!../scripts/ann_feats_up_www2.pl -S -m BB  ../seq/gstm1_human.vaa $FA_DB > results/test2V_m1.ok2mB
echo "done"
echo "starting fastxy36" `date`
../bin/fastx36 -V \!../scripts/ann_feats_up_www2.pl -m 9c -S -q ../seq/mgtt2_x.seq $FA_DB > results/test2V_t2.xk2m9c
../bin/fastx36 -V \!../scripts/ann_feats_up_www2.pl -m BB -S -q ../seq/mgtt2_x.seq $FA_DB > results/test2V_t2.xk2mB
../bin/fastx36 -V \!../scripts/ann_feats_up_www2.pl -m 9c -S -q -z 22 ../seq/gstm1b_human.nt $FA_DB > results/test2V_m1.xk2m9cz22
../bin/fasty36 -V \!../scripts/ann_feats_up_www2.pl -S -q -z 21 ../seq/gstm1b_human.nt $FA_DB > results/test2V_m1.yk2z21
echo "done"
echo "starting ssearch36" `date`
../bin/ssearch36 -V q\!../scripts/ann_pfam_www.pl -V \!../scripts/ann_pfam_www.pl -m 9c -S -z 22 -q ../seq/gstm1_human.vaa $FA_DB > results/test2V_m1.ssm9cz22
../bin/ssearch36 -V q\!../scripts/ann_pfam_www.pl -V \!../scripts/ann_pfam_www.pl -m 9C -S -z 21 -q ../seq/gstm1_human.vaa $FA_DB > results/test2V_m1.ssm9Cz21
../bin/ssearch36 -V q\!../scripts/ann_pfam_www.pl -V \!../scripts/ann_pfam_www.pl -m 8CC -S  -q ../seq/gstm1_human.vaa $FA_DB > results/test2V_m1.ssm8CC
echo "done" `date`
echo "starting ssearch36" `date`
../bin/ggsearch36 -V q\!../scripts/ann_feats_up_www2.pl -V \!../scripts/ann_feats_up_www2.pl -m 9c -S -q ../seq/gstm1_human.vaa $FA_DB > results/test2V_m1.ggm9c
../bin/ggsearch36 -V q\!../scripts/ann_feats_up_www2.pl -V \!../scripts/ann_feats_up_www2.pl -m 9C -S -z 21 -q ../seq/gstm1_human.vaa $FA_DB > results/test2V_m1.ggm9Cz21
echo "done" `date`
