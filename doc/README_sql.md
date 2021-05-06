## Sequence libraries in relational databases

This has been tested in Debian 10 (Buster).

### PostgreSQL

To compile the FASTA programs:
* install the package `libpq-dev` (`sudo apt-get install libpq-dev`)
* `cd src` and `make -f ../make/Makefile.linux_pgsql all`

Library file example:
```
<postgresql hostname>:<port> <db name> <db user> <db password>
DO SELECT 1;
SELECT id,seq FROM sequences;
SELECT description FROM sequences WHERE id=#;
SELECT seq FROM sequences WHERE id=#;
```

Program exection:
```
../bin/fasta36 -q ../seq/mgstm1.aa "/path/to/library_file 17"
```