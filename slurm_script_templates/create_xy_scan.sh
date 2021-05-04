#!/bin/bash

## to creat input and running files for multiple tip-descents        ##
# first try whatever point (in this nomenclature it will be 0 0 ))   #
# to prepare files for the first y line then have - true - in the if #
if true ; then
  Bf="1 2 3 4 5"
  Af=0
  for b in $Bf; do
    echo $b
    y="`echo $b`";
    x="`echo $Af`";
    sed "s/BBB/$y/g"  pdtt.template > pdtt.tmp ;
    sed "s/AAA/$x/g"  pdtt.tmp > pdtt_$b.$Af.py ;
    rm pdtt.tmp
    sed "s/BBB/$y/g"  pdtt_run.template > pdtt_run.tmp ;
    sed "s/AAA/$x/g"  pdtt_run.tmp > rw_$b.$Af.slrm ;
    rm pdtt_run.tmp
  done
fi
# to prepare the rest of the scan (a > 0, b  >= 0) , the have - true - in the if ##
if true ; then
  Bi="0 1 2 3 4 5"
  Ai="1 2 3 4 5"
  for b in $Bi; do
    echo _b_ $b
    for a in $Ai; do 
      echo _a_ $a
      y="`echo $b`";
      x="`echo $a`";
      sed "s/BBB/$y/g"  pdtt.template > pdtt.tmp ;
      sed "s/AAA/$x/g"  pdtt.tmp > pdtt_$b.$a.py ;
      rm pdtt.tmp
      sed "s/BBB/$y/g"  pdtt_run.template > pdtt_run.tmp ;
      sed "s/AAA/$x/g"  pdtt_run.tmp > rw_$b.$a.slrm ;
      rm pdtt_run.tmp
   done
  done
fi

