#!/bin/bash

## to creat input and running files for 2 tip-descents and multiple Voltages ##
# first try whatever point (in this nomenclature it will be 0 0 ))   #
# to prepare files for the first y line then have - true - in the if #
grdb=glob_res/kpfm.db
#global results db file#
cpdb=false 
# copy database #
dbf="kpfm_res_"
if true ; then
  b=0
  Af=0
  #Vf="2.0 1.7 1.5 1.3 1.0 0.7 0.5 0.3 -0.3 -0.5 -0.7 -1.0 -1.3 -1.5 -1.7 -2.0"
  Vf="1000.0 2.0 1.7 1.5 1.3 1 0.7 0.5 0.31 0.3 0.2 0.1 0.05 0.01 -0.01 -0.05 -0.1 -0.2 -0.3 -0.31 -0.5 -0.7 -1 -1.3 -1.5 -1.7 -2.0"
  for Vi in $Vf; do
    echo $b $Af $Vi
    y="`echo $b`";
    x="`echo $Af`";
    V="`echo $Vi`"
    #sed "s/BBB/$y/g"  pdtt_V.template > pdtt.tmp ;
    #sed "s/AAA/$x/g"  pdtt.tmp > pdtt_$b.$Af.tmp ;
    #sed "s/CCC/$V/g"  pdtt_$b.$Af.tmp > pdtt_$b.$Af.$Vi.py;
    #rm pdtt.tmp pdtt_$b.$Af.tmp
    sed "s/BBB/$y/g"  pdtt_V_run.template > pdtt_run.tmp ;
    sed "s/AAA/$x/g"  pdtt_run.tmp > rw_$b.$Af.tmp ;
    sed "s/CCC/$V/g"  rw_$b.$Af.tmp > rw_$b.$Af.$Vi.slrm;
    rm pdtt_run.tmp rw_$b.$Af.tmp
    if $cpdb; then
        tdbf=$dbf
        tdbf+=$b
        tdbf+=_$Af
        tdbf+=_$Vi
        tdbf+=.db
        echo $tdbf
        cp $grdb $tdbf
    fi
  done
fi
if false ; then
  b=0
  Af=1
  Vf="2.0 1.7 1.5 1.3 1.0 0.7 0.5 0.3 -0.3 -0.5 -0.7 -1.0 -1.3 -1.5 -1.7 -2.0"
  #Vf="2.0 1.7 1.5 1.3 0.7 0.5 0.3 -0.3 -0.5 -0.7 -1.3 -1.5 -1.7 -2.0"
  for Vi in $Vf; do
    echo $b $Af $Vi
    y="`echo $b`";
    x="`echo $Af`";
    V="`echo $Vi`"
    sed "s/BBB/$y/g"  pdtt_V.template > pdtt.tmp ;
    sed "s/AAA/$x/g"  pdtt.tmp > pdtt_$b.$Af.tmp ;
    sed "s/CCC/$V/g"  pdtt_$b.$Af.tmp > pdtt_$b.$Af.$Vi.py;
    rm pdtt.tmp pdtt_$b.$Af.tmp
    sed "s/BBB/$y/g"  pdtt_V_run.template > pdtt_run.tmp ;
    sed "s/AAA/$x/g"  pdtt_run.tmp > rw_$b.$Af.tmp ;
    sed "s/CCC/$V/g"  rw_$b.$Af.tmp > rw_$b.$Af.$Vi.slrm;
    rm pdtt_run.tmp rw_$b.$Af.tmp
    if $cpdb; then
        tdbf=$dbf
        tdbf+=$b
        tdbf+=_$Af
        tdbf+=_$Vi
        tdbf+=.db
        echo $tdbf
        cp $grdb $tdbf
    fi
  done
fi
