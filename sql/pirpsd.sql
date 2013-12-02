xdb.wrplab PIRPSD seq_demo demo_pass;
SELECT PIRID, SEQUENCES, PIRID
 FROM c_psdsequence;
SELECT PIRID, concat(PIRID," ",TITLE) FROM c_psdmain
 WHERE PIRID='#';
SELECT PIRID, SEQUENCES, PIRID
 FROM c_psdsequence
 WHERE PIRID='#';
