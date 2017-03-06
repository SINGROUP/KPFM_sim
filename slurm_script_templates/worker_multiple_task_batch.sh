#!/bin/bash
#SBATCH -p parallel
#SBATCH -n 144
#SBATCH --exclusive
#SBATCH --time=2-00:00:00 --mem-per-cpu=7000
#SBATCH --constraint=hsw

module load cp2k/2.7-dev
export PYTHONPATH=$PYTHONPATH:/homeappl/home/juritala/.local/lib/python2.7/site-packages/pycp2k-0.1-py2.7.egg/pycp2k
export REAL_ELPA_KERNEL=REAL_ELPA_KERNEL_AVX_BLOCK2

USAGE="Usage: $0 <task_db_file> <task_type> <task_constraint>"

if [ "$#" != "3" ]; then
    echo "$USAGE"
    exit 1
fi

MINIMUM_TASK_TIME=6

TASK_DB_FILE=$1
TASK_TYPE=$2
TASK_CONSTRAINT=$3

ORIG_DIR=$PWD
cd $TMPDIR/$SLURM_JOB_ID

HOURS_LEFT=100
trap "python $KPFM_GLOBAL_SCRIPTS/call_error_handler.py $SLURM_JOB_ID $ORIG_DIR $TASK_DB_FILE; exit" ERR TERM

while (( HOURS_LEFT >= MINIMUM_TASK_TIME )); do
    cp $CP2K_DATA/BASIS_MOLOPT .
    cp $CP2K_DATA/GTH_POTENTIALS .
    cp $ORIG_DIR/cp2k_system_params.py .

    python $KPFM_GLOBAL_SCRIPTS/run_task.py $TASK_DB_FILE $ORIG_DIR $SLURM_JOB_ID $TASK_TYPE $TASK_CONSTRAINT

    TIME_LEFT_STR=$(squeue -h -j $SLURM_JOB_ID -o %L)
    TIME_LEFT=($(echo $TIME_LEFT_STR | tr '-' ' ' | tr ':' ' '))
    echo "Time left: $TIME_LEFT_STR"

    if (( ${#TIME_LEFT[@]} < 3)); then
        HOURS_LEFT=0
    elif (( ${#TIME_LEFT[@]} == 3 )); then
        HOURS_LEFT=10#${TIME_LEFT[0]}
    else
        HOURS_LEFT=$((${TIME_LEFT[0]}*24+10#${TIME_LEFT[1]}))
    fi
    echo "Hours left: $HOURS_LEFT"
    
    rm *
done
