rem ""
rem "starting fasta36_t - protein on win32"
rem ""
..\bin\fasta36_t -q -m 6 -Z 100000 ..\seq\mgstm1.aa:1-100 q > results\test_m1.ok2_t.html
..\bin\fasta36_t -S -q -z 11 -O results\test_m1.ok2_t_p25 -s P250 ..\seq\mgstm1.aa:100-218 q
rem "done"
rem "starting fastxy36_t"
..\bin\fastx36_t -m 9c -S -q ..\seq\mgtt2_x.seq q 1 > results\test_t2.xk1_t
..\bin\fasty36_t -S -q ..\seq\mgtt2_x.seq q > results\test_t2.yk2_t
..\bin\fastx36_t -m 9c -S -q -z 2 ..\seq\mgstm1.esq a > results\test_m1.xk2_tz2
..\bin\fasty36_t -S -q -z 2 ..\seq\mgstm1.esq a > results\test_m1.yk2_tz2
rem "done"
rem "starting fastxy36_t rev"
..\bin\fastx36_t -m 9c -q -m 5 ..\seq\mgstm1.rev q > results\test_m1.xk2r_t
..\bin\fasty36_t -q -m 5 -M 200-300 -z 2 ..\seq\mgstm1.rev q > results\test_m1.yk2r_tz2
..\bin\fasty36_t -q -m 5 -z 11 ..\seq\mgstm1.rev q > results\test_m1.yk2rz11_t
rem "done"
rem "starting ssearch36_t"
..\bin\ssearch36_t -m 9c -S -z 3 -q ..\seq\mgstm1.aa  q > results\test_m1.ss_tz3
..\bin\ssearch36_t -q -M 200-300 -z 2 -Z 100000 -s P250 ..\seq\mgstm1.aa q > results\test_m1.ss_t_p25
rem "done"
rem "starting prss/prfx36"
..\bin\ssearch36_t -q -k 1000 -A ..\seq\mgstm1.aa ..\seq\xurt8c.aa  > results\test_m1.rss
..\bin\fastx36_t -q -k 1000 -A ..\seq\mgstm1.esq ..\seq\xurt8c.aa > results\test_m1.rfx
rem "done"
rem "starting fasta36_t - DNA"
..\bin\fasta36_t -S -q -z 2 ..\seq\mgstm1.seq %M 4 > results\test_m1.ok4_tz2
..\bin\fasta36_t -S -q ..\seq\mgstm1.rev %M 4 > results\test_m1.ok4r_t
rem "done"
rem "starting tfastxy36_t"
..\bin\tfastx36_t -m 9c -q -i -3 -m 6 ..\seq\mgstm1.aa %m > results\test_m1.tx2_t.html
..\bin\tfasty36_t -q -i -3 -N 5000 ..\seq\mgstm1.aa %m > results\test_m1.ty2_t
rem "done"
rem "starting fastf36_t"
..\bin\fastf36_t -q ..\seq\m1r.aa q > results\test_mf.ff_t
..\bin\fastf36 -q ..\seq\m1r.aa q > results\test_mf.ff_s
rem "done"
rem "starting tfastf36_t"
..\bin\tfastf36_t -q ..\seq\m1r.aa %m > results\test_mf.tf_tr
rem "done"
rem "starting fasts36_t"
..\bin\fasts36_t -q -V '*?@' ..\seq\ngts.aa q > results\test_m1.fs1_t
..\bin\fasts36_t -q ..\seq\ngt.aa q > results\test_m1.fs_t
..\bin\fasts36_t -q -n ..\seq\mgstm1.nts m > results\test_m1.nfs_t
rem "done"
rem "starting tfasts36_t"
..\bin\tfasts36_t -q ..\seq\n0.aa %m > results\test_m1.ts_r
rem "done"
rem "starting fasta36 - protein"
..\bin\fasta36 -q -z 2 ..\seq\mgstm1.aa q 1 > results\test_m1.ok1z2
..\bin\fasta36 -q -s P250 ..\seq\mgstm1.aa q > results\test_m1.ok2_p25 
rem "done"
rem "starting fastx3"
..\bin\fastx36 -m 9c -q ..\seq\mgstm1.esq q > results\test_m1.ok2x 
rem "done"
rem "starting fasty3"
..\bin\fasty36 -q ..\seq\mgstm1.esq q > results\test_m1.ok2y 
rem "done"
rem "starting fasta36 - DNA "
..\bin\fasta36 -m 9c -q ..\seq\mgstm1.seq M 4 > results\test_m1.ok4 
rem "done"
rem "starting ssearch3"
..\bin\ssearch36 -S -q -z 2 ..\seq\mgstm1.aa q > results\test_m1.ss_z2
..\bin\ssearch36 -q -s P250 ..\seq\mgstm1.aa q > results\test_m1.ss_p25 
rem "done"
rem "starting tfastxy3"
..\bin\tfastx36 -q ..\seq\mgstm1.aa M > results\test_m1.tx2 
..\bin\tfasty36 -m 9c -q ..\seq\mgstm1.aa M > results\test_m1.ty2 
rem "done"
rem "starting fasts36"
..\bin\fasts36 -q -V '@?*' ..\seq\ngts.aa q > results\test_m1.fs1
..\bin\fasts36 -q ..\seq\ngt.aa q > results\test_m1.fs
rem "done"
