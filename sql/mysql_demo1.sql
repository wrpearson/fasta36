seqdb_host seqdb_demo seqdb_user password;
SELECT acc, protein.seq, sp_name
 FROM annot INNER JOIN protein USING(prot_id) WHERE annot.db='sp' LIMIT 50000;
SELECT acc, concat('sp|',acc,'|',sp_name,' ',descr) FROM annot WHERE acc='#' AND db='sp';
SELECT acc,protein.seq FROM protein INNER JOIN annot USING(prot_id)
 WHERE annot.acc='#' AND db='sp';
