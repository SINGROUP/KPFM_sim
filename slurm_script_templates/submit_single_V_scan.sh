#!/bin/bash

## to creat input and running files for 2 tip-descents with different voltages ##
# first try whatever point (in this nomenclature it will be 0 0 ))   #
# to prepare files for the first y line then have - true - in the if #
if true ; then
  b=0
  Af=0
  Vf="1000.0 2.0 1.7 1.5 1.3 1 0.7 0.5 0.31 0.3 0.2 0.1 0.05 0.01 -0.01 -0.05 -0.1 -0.2 -0.3 -0.31 -0.5 -0.7 -1 -1.3 -1.5 -1.7 -2.0"
  for Vi in $Vf; do
    sbatch rw_$b.$Af.$Vi.slrm ;
  done
fi
if false ; then
  b=0
  Af=1
  Vf="2.0 1.7 1.5 1.3 1 0.7 0.5 0.3 -0.3 -0.5 -0.7 -1.0 -1.3 -1.5 -1.7 -2.0"
  for Vi in $Vf; do
    echo   rw_$b.$Af.$Vi.slrm ;
    sbatch rw_$b.$Af.$Vi.slrm ;
  done
fifi
