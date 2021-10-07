#!/bin/bash

#######################################################################
## to submit computations for multiple tip-descents                  ##
# first try whatever point (in this nomenclature it will be 0 0 ))   #
# to prepare files for the first y line then have - true - in the if #
######################################################################
# This example is for (y) 3 x (x) 15 scan points for hBN AFM calc    #
#********************************************************************##
# first line over y ##
if true ; then
  Bf="1 2"
  Af=0
  for b in $Bf; do
    echo $b
    y="`echo $b`";
    x="`echo $Af`";
    sbatch rw_$b.$Af.slrm ;
  done
fi
# to prepare the rest of the scan (a > 0, b  >= 0) , the have - true - in the if ##
if true ; then
  Bi="0 1 2"
  Ai="1 2 3 4 5 6 7 8 9 10 11 12 13 14"
  for b in $Bi; do
    echo _b_ $b
    for a in $Ai; do 
      echo _a_ $a
      y="`echo $b`";
      x="`echo $a`";
      sbatch rw_$b.$a.slrm ;
   done
  done
fi
