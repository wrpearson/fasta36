#!/bin/sh

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
done
