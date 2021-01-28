#!/bin/bash
module purge
module load cp2k/6.1-openmpi-scalapack-popt anaconda/2020-04-tf2
export OMP_NUM_THREADS=1
ulimit -s unlimited

#source /home/krejcio1/.bashrc
eval "$(conda shell.bash hook)"
conda activate KPFM1
progdir=/home/krejcio1/wrkdir/Programs/ana-20.04tf2_py3.7.7

export PYTHONPATH=$PYTHONPATH:$progdir/KPFM_sim3:$progdir/CP2k_mtools3:$progdir/DFT_gridIO3

echo "start combining results from different workers"

# put results of worker_1 and worker_2 together
# Usage: python combine_workers_to_global.py (optionally -k or --kpts if k-points are necessary) <worker1_res_db_file> <worker1_path> <worker2_res_db_file> <worker2_path> <to: global_results_db_file> <global_results_wfn_path>  
python3 combine_workers_to_global.py afm_res.db worker_1 afm_res2.db worker_2 glob_res/afm.db wfn_data

# Usage: python (optionally -k or --kpts if k-points are necessary) copy_scan_points_and_wfn_files.py from_db_file <worker_path> to_db_file <results_wfn_path> 
# put results from worker_y-Cl to the global results
python3 copy_scan_points_and_wfn_files.py afm_res_y-Cl.db worker_y-Cl/ glob_res/afm.db wfn_data

echo "Done, done"

