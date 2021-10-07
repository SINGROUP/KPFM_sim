#!/bin/bash

#######################################################################
## to control calculations for  multiple tip-descents                ##
# first try whatever point (in this nomenclature it will be 0 0 ))   #
# to prepare files for the first y line then have - true - in the if #
######################################################################
# This example is for (y) 3 x (x) 15 scan points for hBN AFM calc    #
#********************************************************************##

name="krejcion" # user name

squeue -u $name

echo "***** All running jobs above, now all the direct results vs. running jobs are bellow *****"

for i in {0..2};
do
    for j in {0..14}
    do
       str="afm_res_$i"
       str+="_$j.db"
       # echo $i $j $str
       file=`grep "result_file" sbatch-*.out | grep $str | tail -1 | cut -d':' -f 1`
       done=`grep "Done, done" $file`
       error=`grep "=== ERROR ==" $file`
       echo $i $j $file $done $error
       if test -z "$done"
       then
           if test -z "$error"
           then
               tmp=`echo $file | cut -d'.' -f 1 | cut -d'-' -f 2`
               squeue -u $name | grep $tmp
           fi
       fi
    done
done
