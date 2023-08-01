#!/bin/bash
#SBATCH --account=juritala
#SBATCH -p medium           # test - for testing 1h ; medium - up to 20 nodes/36 hours ; large - 20-200 nodes/36 hours #
#SBATCH --time=12:00:00     # dd-hh:mm:ss
#SBATCH -J npy2cube         # name 
#SBATCH -o sbatch-%j.out    # where the outputs & errors are written
#SBATCH -e sbatch-%j.err    # where the errpr is written
#SBATCH --nodes=1           # n processes to run (N x 64 = n); max 192 ; max debug 48
#SBATCH --ntasks-per-node=1

# check modules.txt for updating following lines: # module purge
module purge
module load python-env/3.8.6 openmpi/4.0.3 gcc/9.3.0
#module load gcc/9.3.0  openmpi/4.0.3 cp2k/7.1-omp   # The full cp2k is not needed any more #
export OMP_NUM_THREADS=1
ulimit -s unlimited
export ppath=./programs_KPFM
export PYTHONPATH=$PYTHONPATH:$ppath/CP2k_mtools3:$ppath/DFT_gridIO3:$ppath/KPFM_FEM3:$ppath/KPFM_sim3:$ppath/ase
KPFM_GLOBAL_SCRIPTS=$ppath/KPFM_sim3/scripts
ORIG_DIR=`pwd`

project_path=./
slurm_id=$SLURM_JOB_ID

# ******************** safety check *************************
echo " *** safety check ***"
echo "PYTHONPATH:" $PYTHONPATH
echo "KPFM_GLOBAL_SCRIPTS:" $KPFM_GLOBAL_SCRIPTS
echo "ORIG_DIR:" $ORIG_DIR
echo "glob_file:" $glob_file
echo "worker_name:" $worker_name
echo "task_db_file" $task_db_file
echo "slurm_id (2x)" $slurm_id $SLURM_JOB_ID
echo "result_file" $result_file
echo " *** safety check ended ***"
# ***************** safety check ended **********************

python3 npy_to_cube.py

echo "Done, done"
