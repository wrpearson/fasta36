
DROP DATABASE seq_demo;
CREATE DATABASE seq_demo;

USE seq_demo;

CREATE TABLE prot (
id INT UNSIGNED NOT NULL AUTO_INCREMENT PRIMARY KEY,
seq TEXT NOT NULL,
bin BLOB NOT NULL,
len INT UNSIGNED NOT NULL
);

CREATE TABLE annot (
prot_id INT UNSIGNED NOT NULL,
gi INT UNSIGNED NOT NULL PRIMARY KEY,
db ENUM("gb","emb","dbj","prf","ref","pdb","pir","sp") NOT NULL,
descr TEXT NOT NULL,

INDEX (prot_id),
INDEX (db)
);

CREATE TABLE sp (
    gi INT UNSIGNED NOT NULL,
    acc   VARCHAR(10),
    name  VARCHAR(10),

    PRIMARY KEY (gi)
);
