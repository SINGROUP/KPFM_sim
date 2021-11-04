#!/bin/bash
module load gcc/9.3.0  openmpi/4.0.3 cp2k/7.1-omp 
export OMP_NUM_THREADS=1
ulimit -s unlimited
export ppath=/scratch/project_2003835/programs_KPFM
export PYTHONPATH=$PYTHONPATH:$ppath/CP2k_mtools3:$ppath/DFT_gridIO3:$ppath/KPFM_FEM3:$ppath/KPFM_sim3
KPFM_GLOBAL_SCRIPTS=$ppath/KPFM_sim3/scripts

# **************************************************************************************************** #
#                                                                                                      #
# script for putting all the results from all the (AFM) workers to be saved in the global results db   #
#                                                                                                      #
# **************************************************************************************************** #

echo "start combining results from different workers"

GR='glob_res/afm_copy.db' # adjust the global results name #
bu=$GR
bu+='.bak'

echo $GR $bu

cp $GR $bu


for i in {0..3}; # adjust amount of y points #
do
    for j in {0..14}; # adjust amout of x points #
    do
         ar=afm_res_$i
         ar+=_$j.db  # adjust the input db file name(s) #
         echo $ar
         python3 $KPFM_GLOBAL_SCRIPTS/copy_scan_points_and_wfn_files.py -i $ar -o $GR # -w worker_0_0/ wfn_data # -k
    done
done

echo "Done, done"

