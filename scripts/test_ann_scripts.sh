#!/bin/sh

## a script to test the annotation scripts ann_*.pl
## acc_examples contains:
# P09488  -- new NCBI format
# sp|P09488 --more traditional format
# up|P09488|GSTM1_HUMAN
# SP:GSTM1_HUMAN P09488  ---ebi searches with accession
# SP:GSTM1_HUMAN -- ebi searches without accession
##

if [ ! $1=='' ]; then
   script_file=$1
else
    script_file=ann_script_list
fi

if [ ! $1=='' ]; then
   ex_file=$1
else
    ex_file=acc_examples
fi

for script in `cat $script_file `; do
  for acc_type in `cat $ex_file`; do
      echo $script ${acc_type}
      $script ${acc_type}
  done
  echo '***DONE***'  $script `date`
done
