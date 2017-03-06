#!/bin/bash
#SBATCH -p parallel
#SBATCH -n 144
#SBATCH --exclusive
#SBATCH --time=2-00:00:00 --mem-per-cpu=4000
#SBATCH --constraint=hsw

module load cp2k/2.7-dev
export PYTHONPATH=$PYTHONPATH:/homeappl/home/juritala/.local/lib/python2.7/site-packages/pycp2k-0.1-py2.7.egg/pycp2k
export REAL_ELPA_KERNEL=REAL_ELPA_KERNEL_AVX_BLOCK2

USAGE="Usage: $0 <task_db_file> <task_type> <task_constraint>"

if [ "$#" != "3" ]; then
    echo "$USAGE"
    exit 1
fi

TASK_DB_FILE=$1
TASK_TYPE=$2
TASK_CONSTRAINT=$3

ORIG_DIR=$PWD
cd $TMPDIR/$SLURM_JOB_ID

cp $CP2K_DATA/BASIS_MOLOPT .
cp $CP2K_DATA/GTH_POTENTIALS .
cp $ORIG_DIR/cp2k_system_params.py .

trap "python $KPFM_GLOBAL_SCRIPTS/call_error_handler.py $SLURM_JOB_ID $ORIG_DIR $TASK_DB_FILE; exit" ERR TERM

python $KPFM_GLOBAL_SCRIPTS/run_task.py $TASK_DB_FILE $ORIG_DIR $SLURM_JOB_ID $TASK_TYPE $TASK_CONSTRAINT
