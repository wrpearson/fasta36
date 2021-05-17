## Sequence libraries in relational databases

This has been tested in Debian 10 (Buster).

### Library file example
```
<postgresql/mariadb hostname>:<port> <db name> <db user> <db password>
DO SELECT 1;
SELECT id,seq FROM sequences;
SELECT description FROM sequences WHERE id=#;
SELECT seq FROM sequences WHERE id=#;
```

### PostgreSQL

To compile the FASTA programs:
* install the package `libpq-dev` (`sudo apt-get install libpq-dev`)
* `cd src` and `make -f ../make/Makefile.linux_pgsql all`

Program execution example:
```
../bin/fasta36 -q ../seq/mgstm1.aa "/path/to/library_file 17"
```

### MariaDB (fork of MySQL)

To compile the FASTA programs:
* install the package `libmariadb-dev` (`sudo apt-get install libmariadb-dev`)
* `cd src` and `make -f ../make/Makefile.linux_mariadb all`

Program execution example:
```
../bin/fasta36 -q ../seq/mgstm1.aa "/path/to/library_file 16"
```