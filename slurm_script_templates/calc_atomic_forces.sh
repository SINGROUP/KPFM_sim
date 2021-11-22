#!/bin/bash
#SBATCH --account=project_2003835
#SBATCH -p medium # test - for testing 1h ; medium - up to 20 nodes/36 hours ; large - 20-200 nodes/36 hours #
#SBATCH --time=36:00:00      # dd-hh:mm:ss
#SBATCH -J forces        # name 
#SBATCH -o sbatch-%j.out # where the outputs & errors are written
#SBATCH -N 5                  # N nodes to run (N x 64 = n); max 192 ; max debug 48
#SBATCH --ntasks-per-node=128  # to run properly!!!!


# check modules.txt for updating following lines: # module purge
module load gcc/9.3.0  openmpi/4.0.3 cp2k/7.1-omp 
export OMP_NUM_THREADS=1
ulimit -s unlimited
export ppath=/scratch/project_2003835/programs_KPFM
export PYTHONPATH=$PYTHONPATH:$ppath/CP2k_mtools3:$ppath/DFT_gridIO3:$ppath/KPFM_FEM3:$ppath/KPFM_sim3
KPFM_GLOBAL_SCRIPTS=$ppath/KPFM_sim3/scripts
ORIG_DIR=`pwd`

worker_name=worker_forces
orig_glob=glob_res/afm_copy.db
result_file=glob_res/afm_forces.db
glob_file=$orig_glob

# copy only if the file does not exists: 
if test -f "$result_file"; then
    echo "$result_file exists."
else
    cp $orig_glob $result_file
fi
if [ ! -d $worker_name ]; then
    mkdir $worker_name
    cd $worker_name
    ln -s ../cp2k_common_params.py ./cp2k_common_params.py
    ln -s ../cp2k_system_params.py ./cp2k_system_params.py
else
    cd $worker_name
fi
#

project_path=./
slurm_id=$SLURM_JOB_ID

# ******************** safety check *************************
echo " *** safety check ***"
echo "PYTHONPATH:" $PYTHONPATH
echo "KPFM_GLOBAL_SCRIPTS:" $KPFM_GLOBAL_SCRIPTS
echo "ORIG_DIR:" $ORIG_DIR
echo "worker_name:" $worker_name
echo "glob_file" $glob_file
echo "slurm_id (2x)" $slurm_id $SLURM_JOB_ID
echo "result_file" $result_file
echo " *** safety check ended ***"
# ***************** safety check ended **********************

# Usage: python calc_atomic_forces.py -p <project_path> -r <result_db_file> -g <global_res_db_file> (-k/--kpts if k-pouints in cp2k calc) (-l/--location <x> <y> <V>);39M;39m
    
python3 $KPFM_GLOBAL_SCRIPTS/calc_atomic_forces.py -p ../ -r $result_file -g $glob_file -l 12.3707 -3.1602 0.0

echo "Done, done"

