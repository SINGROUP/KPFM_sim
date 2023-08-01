#!/bin/bash

## to creat input and running files for multiple tip-descents        ##
# first try whatever point (in this nomenclature it will be 0 0 ))   #
# to prepare files for the first y line then have - true - in the if #
#
#mpicE_run.template  mpicE.template
#
if true ; then
  for i in {0..20}; do
    echo $i
    ix="`echo $i`";
    sed "s/AAAA/$ix/g"  mpicE_run.template > mpicE_$i.sh ;
    sed "s/AAAA/$ix/g"  mpicE.template > mpicE_$i.py ;
  done
fi
