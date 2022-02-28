#!/bin/bash
module load gcc/9.3.0  openmpi/4.0.3 cp2k/7.1-omp 
export OMP_NUM_THREADS=1
ulimit -s unlimited
export ppath=/scratch/project_2003835/programs_KPFM
export PYTHONPATH=$PYTHONPATH:$ppath/CP2k_mtools3:$ppath/DFT_gridIO3:$ppath/KPFM_FEM3:$ppath/KPFM_sim3
KPFM_GLOBAL_SCRIPTS=$ppath/KPFM_sim3/scripts

# **************************************************************************************************** #
#                                                                                                      #
# script for putting all the results from all the KPFM workers to be saved in the global results db    #
#                                                                                                      #
# **************************************************************************************************** #

echo "start combining results from different workers"

GR='glob_res/afm_copy.db' # adjust the global results name #
bu=$GR
bu+='.bak'

echo $GR $bu

cp $GR $bu

# adjust the voltages, where it was calculated  #
vol="0.31 0.2 -0.2 -0.31 -2.0 -1.7 -1.5 -1.3 -1 -0.7 -0.5 -0.3 -0.1 -0.05 -0.01 0.01 0.05 0.1 0.3 0.5 0.7 1 1.3 1.5 1.7 2.0"
vol+=" 1000.0"


for i in {0..0}; # adjust amount of y points #
do
    for j in {0..0}; # adjust amout of x points #
    do
        for vi in $vol;
        do
            ar=kpfm_res_$i
            ar+=_$j
            ar+=_$vi.db  # adjust the input db file name(s) #
            echo $ar
            wr=worker_$i
            wr+=_$j
            wr+=_$vi
            python3 $KPFM_GLOBAL_SCRIPTS/copy_scan_points_and_wfn_files.py -i $ar -o $GR # -w $vi wfn_data # -k ## -k -- if kpts are necessary ##
        done
    done
done

echo "Done, done"

