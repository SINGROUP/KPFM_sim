#!/bin/bash
#SBATCH -p short             # debug - for testing ; short - up to 4 hours ; batch - up to 5 days #
#SBATCH --time=04:00:00      # dd-hh:mm:ss
#SBATCH --mem-per-cpu=1000   # 1000 - 3500 nromally, but can be even higher
#SBATCH -J k-yCl             # name 
#SBATCH -o slurm-%j.out      # where the outputs & errors are written
#SBATCH --constraint=csl     # require Haswell CPUs with 24 cores per node (ivb, hsw, vdw, skl }= avx;  hsw, vdw, skl }= avx2 )  
#SBATCH -N 1
#SBATCH -n 16                # n processes to run (N x 24 = n); max 192 ; max debug 48
    
module purge
module load cp2k/6.1-openmpi-scalapack-popt anaconda/2020-04-tf2
export OMP_NUM_THREADS=1
ulimit -s unlimited

eval "$(conda shell.bash hook)"
conda activate KPFM1
progdir=/home/krejcio1/wrkdir/Programs/ana-20.04tf2_py3.7.7

export PYTHONPATH=$PYTHONPATH:$progdir/KPFM_sim3:$progdir/CP2k_mtools3:$progdir/DFT_gridIO3
KPFM_GLOBAL_SCRIPTS=/scratch/work/krejcio1/Programs/ana-20.04tf2_py3.7.7/KPFM_sim3/scripts
ORIG_DIR=`pwd`

worker_name=worker_y-Cl
glob_db=glob_res/afm.db
result_file=afm_res_y-Cl.db
task_db_file=afm_workers.db

mkdir $worker_name
cp $glob_db ./$result_file
cd $worker_name
cp ../$task_db_file ./$task_db_file
ln -s ../cp2k_common_params.py ./cp2k_common_params.py
ln -s ../cp2k_system_params.py ./cp2k_system_params.py

# ############### #
# for run_task.py #
#if len(sys.argv) == 4:
#    task_db_file = sys.argv[1]
#    project_path = sys.argv[2]
#    slurm_id = int(sys.argv[3])
#    task_type_constraint = None
#    task_state_constraint = default_state_constraint
# ############### #

project_path=./
slurm_id=$SLURM_JOB_ID

# ******************** safety check *************************
echo " *** safety check ***"
echo "PYTHONPATH:" $PYTHONPATH
echo "KPFM_GLOBAL_SCRIPTS:" $KPFM_GLOBAL_SCRIPTS
echo "ORIG_DIR:" $ORIG_DIR
echo "worker_name:" $worker_name
echo "glob_db" $glob_db
echo "task_db_file" $task_db_file
echo "slurm_id (2x)" $slurm_id $SLURM_JOB_ID
echo "result_file" $result_file
echo " *** safety check ended ***"
# ***************** safety check ended **********************


#python3 afm_opt_y.py

#./cp_opt.sh

python3 ../plan_descend_tip_task_y-Cl.py --kpts ./$task_db_file ../$result_file ../$glob_db

trap "python3 $KPFM_GLOBAL_SCRIPTS/call_error_handler.py $SLURM_JOB_ID $ORIG_DIR $TASK_DB_FILE; exit" ERR TERM

python3 ../run_task.py --kpts ./$task_db_file $project_path $slurm_id

echo "Done, done"

#echo "arguments to srun" $cp2k_bin $output_file $input_file
#echo "Name of the partition in which the job is running: $SLURM_JOB_PARTITION"
#echo "Name of the node running the job script: $SLURMD_NODENAME"
#echo "     The ID of the job allocation: $SLURM_JOB_ID"
#echo "     List of nodes allocated to the job: $SLURM_JOB_NODELIST"
