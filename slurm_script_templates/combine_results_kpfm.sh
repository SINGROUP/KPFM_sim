#!/bin/bash
module load gcc/9.3.0  openmpi/4.0.3 cp2k/7.1-omp 
export OMP_NUM_THREADS=1
ulimit -s unlimited
export ppath=/scratch/project_2003835/programs_KPFM
export PYTHONPATH=$PYTHONPATH:$ppath/CP2k_mtools3:$ppath/DFT_gridIO3:$ppath/KPFM_FEM3:$ppath/KPFM_sim3
KPFM_GLOBAL_SCRIPTS=$ppath/KPFM_sim3/scripts

echo "start combining results from different workers"

#files='kpfm_res_0_0_-0.3.db  kpfm_res_0_0_-0.5.db  kpfm_res_0_0_-0.7.db  kpfm_res_0_0_-1.3.db  kpfm_res_0_0_-1.5.db  kpfm_res_0_0_-1.7.db  kpfm_res_0_0_-1.db  kpfm_res_0_0_-2.0.db
#kpfm_res_0_0_0.3.db   kpfm_res_0_0_0.5.db   kpfm_res_0_0_0.7.db   kpfm_res_0_0_1.3.db   kpfm_res_0_0_1.5.db   kpfm_res_0_0_1.7.db   kpfm_res_0_0_1.db   kpfm_res_0_0_2.0.db'

vol="0.31 0.2 -0.2 -0.31"  #'-2.0 -1.7 -1.5 -1.3 -1.0 -0.7 -0.5 -0.3 0.3 0.5 0.7 1 1.3 1.5 1.7 2.0'

GR='glob_res/afm.db'
#GR='e_fields.db'
bu=$GR
bu+='.bak'

echo $GR $bu
cp $GR $bu

for vi in $vol;
do
    #for j in {0..10};
    #do
         ar=kpfm_res_0_0_$vi
         ar+=.db
         wr=worker_0_0_$vi
         echo $ar $wr
         python3 $KPFM_GLOBAL_SCRIPTS/copy_scan_points_and_wfn_files.py -k -i $ar -o $GR -w $wr wfn
    #done
done


# put results of worker_1 and worker_2 together
# Usage: python combine_workers_to_global.py (optionally -k or --kpts if k-points are necessary) <worker1_res_db_file> <worker1_path> <worker2_res_db_file> <worker2_path> <to: global_results_db_file> <global_results_wfn_path>  
# python3 combine_workers_to_global.py --kpts afm_res.db worker_1 afm_res2.db worker_2 glob_res/afm.db wfn_data

## Usage: python (optionally -k or --kpts if k-points are necessary) copy_scan_points_and_wfn_files.py from_db_file <worker_path> to_db_file <results_wfn_path> 
## put results from worker_y-Cl to the global results
#python3 copy_scan_points_and_wfn_files.py --kpts afm_res_y-Cl.db worker_y-Cl/ glob_res/afm.db wfn_data

echo "Done, done"

