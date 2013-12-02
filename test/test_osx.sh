#!/bin/csh -f
echo ""
echo "starting fasta36_t - protein" `date` "on" `hostname`
echo `uname -a`
echo ""
fasta36_t -q -m 6 -Z 100000 mgstm1.aa:1-100 q > test_m1.ok2_t.html
fasta36_t -S -q -z 11 -O test_m1.ok2_t_p25 -s P250 mgstm1.aa:100-218 q
echo "done"
echo "starting fastxy36_t" `date`
fastx36_t -m 9 -S -q mgtt2_x.seq q 1 > test_t2.xk1_t
fasty36_t -S -q mgtt2_x.seq q > test_t2.yk2_t
fastx36_t -m 9 -S -q -z 2 mgstm1.esq a > test_m1.xk2_tz2
fasty36_t -S -q -z 2 mgstm1.esq a > test_m1.yk2_tz2
echo "done"
echo "starting fastxy36_t rev" `date`
fastx36_t -m 9 -q -m 5 mgstm1.rev q > test_m1.xk2r_t
fasty36_t -q -m 5 -M 200-300 -z 2 mgstm1.rev q > test_m1.yk2r_tz2
fasty36_t -q -m 5 -z 11 mgstm1.rev q > test_m1.yk2rz11_t
echo "done"
echo "starting ssearch36_t" `date`
ssearch36_t -m 9 -S -z 3 -q mgstm1.aa  q > test_m1.ss_tz3
ssearch36_t -q -M 200-300 -z 2 -Z 100000 -s P250 mgstm1.aa q > test_m1.ss_t_p25
echo "done"
echo "starting ssearch36_t" `date`
ssearch36s_t -m 9 -S -z 3 -q mgstm1.aa  q > test_m1.sss_tz3
ssearch36s_t -q -M 200-300 -z 2 -Z 100000 -s P250 mgstm1.aa q > test_m1.sss_t_p25
echo "done"
echo "starting prss36(ssearch/fastx)" `date`
ssearch36_t -q -k 1000 -A mgstm1.aa xurt8c.aa  > test_m1.rss
fastx36_t -q -k 1000 -A mgstm1.esq xurt8c.aa > test_m1.rfx
echo "done"
echo "starting ggsearch36/glsearch36" `date`
ggsearch36_t -q -m 9i -w 80 hahu.aa q > test_h1.gg_t
glsearch36_t -q -m 9i -w 80 hahu.aa q > test_h1.gl_t
ggsearch36_t -q gtt1_drome.aa q > test_t1.gg_t
glsearch36_t -q gtt1_drome.aa q > test_t1.gl_t
echo "done"
echo "starting fasta36_t - DNA" `date`
fasta36_t -q -z 2 mgstm1.seq %M 4 > test_m1.ok4_tz2
fasta36_t -q mgstm1.rev %M 4 > test_m1.ok4r_t
echo "done"
echo "starting tfastxy36_t" `date`
tfastx36_t -m 9 -q -i -3 -m 6 mgstm1.aa %m > test_m1.tx2_t.html
tfasty36_t -q -3 -N 5000 mgstm1.aa %m > test_m1.ty2_t
echo "done"
echo "starting fastf36_t" `date`
fastf36_t -q m1r.aa q > test_mf.ff_s
echo "done"
echo "starting tfastf36_t" `date`
tfastf36_t -q m1r.aa %m > test_mf.tf_r
echo "done"
echo "starting fasts36_t" `date`
fasts36_t -q n0.aa q > test_m1.fs_s
echo "done"
echo "starting tfasts36_t" `date`
tfasts36_t -q n0.aa %m > test_m1.ts_r
echo "done"
echo "done with threaded section: " `date`
echo "starting lalign36" `date`
lalign36 -k 1000 -q mchu.aa mchu.aa > test_mc.lal
lalign36 -z 3 -q mchu.aa mchu.aa > test_mc.lal_z3
lalign36 -s BL62 -f -11 -g -1  -q mchu.aa mchu.aa > test_mc.lal_bl62
lalign36 -k 1000 -q mwkw.aa mwkw.aa > test_mw.lal
lalign36 -z 3 -q mwkw.aa mwkw.aa > test_mw.lal_z3
lalign36 -s BL62 -f -11 -g -1  -q mwkw.aa mwkw.aa > test_mw.lal_bl62
echo "done"
echo "starting fasta36 - protein" `date`
fasta36 -q -z 2 mgstm1.aa q > test_m1.ok2z2
fasta36 -q -s P250 mgstm1.aa q > test_m1.ok2_p25 
echo "done"
echo "starting fastx36" `date`
fastx36 -m 9 -q mgstm1.esq q > test_m1.ok2x 
echo "done"
echo "starting fasty36" `date`
fasty36 -q mgstm1.esq q > test_m1.ok2y 
echo "done"
echo "starting fasta36 - DNA " `date`
fasta36 -m 9 -q mgstm1.seq %m 4 > test_m1.ok4 
echo "done"
echo "starting ssearch36" `date`
ssearch36 -S -q -z 2 mgstm1.aa q > test_m1.ss_z2
ssearch36 -q -s P250 mgstm1.aa q > test_m1.ss_p25 
echo "done"
echo "starting ssearch36" `date`
ssearchs36 -S -q -z 2 mgstm1.aa q > test_m1.sss_z2
ssearchs36 -q -s P250 mgstm1.aa q > test_m1.sss_p25 
echo "done"
echo "starting tfastxy36" `date`
tfastx36 -q mgstm1.aa %m > test_m1.tx2 
tfasty36 -m 9 -q mgstm1.aa %m > test_m1.ty2 
echo "done"
echo "starting fasts36" `date`
fasts36 -q -V '@?*' ngts.aa q > test_m1.fs1
fasts36 -q ngt.aa q > test_m1.fs
echo "done" `date`
