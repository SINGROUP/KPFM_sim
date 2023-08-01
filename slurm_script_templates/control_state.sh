#!/bin/bash

name="krejcion" # user name

squeue -u $name

### AFM adjustable part here: ### 

check_AFM=false   # true/false - whether to check AFM part only
ienda=10          # i will check range(0,ienda)
jenda=10          # j will check range(0,jenda)
str_ba="afm_res_" #beginning of the results db file#

### KPFM adjustable part here: ###

check_KPFM=true   # true/false - whether to check AFM part only
Vf="1000.0 2.0 1.7 1.5 1.3 1 0.7 0.5 0.31 0.3 0.2 0.1 0.05 0.01 -0.01 -0.05 -0.1 -0.2 -0.3 -0.31 -0.5 -0.7 -1 -1.3 -1.5 -1.7 -2.0"
iendk=0          # i will check range(0,iendk)
jendk=0          # j will check range(0,jendk)
str_bk="kpfm_res_" #beginning of the results db file#

### AFM Running part here: ###

if $check_AFM;
then
for i in $(seq 0 $ienda);
do
    for j in $(seq 0 $jenda)
    do
        str=$str_ba
        str+="$i"
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
fi

### KPFM running part here: ###

if $check_KPFM;
then
for i in $(seq 0 $iendk);
do
    for j in $(seq 0 $jendk);
    do 
        for Vi in $Vf
        do
            str=$str_bk
            str+="$i"
            str+="_$j"
            str+="_$Vi.db"
            # echo $i $j $Vi $str  
            file=`grep "result_file" sbatch-*.out | grep $str | tail -1 | cut -d':' -f 1`
            done=`grep "Done, done" $file`
            error=`grep "=== ERROR ==" $file`
            echo $i $j $Vi $file $done $error
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
done
fi
### ****** ****** ###

echo "That is all!"
