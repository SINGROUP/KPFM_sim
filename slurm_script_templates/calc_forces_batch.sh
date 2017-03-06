#!/bin/bash
#SBATCH -p parallel
#SBATCH -n 144
#SBATCH --exclusive
#SBATCH --time=18:00:00 --mem-per-cpu=4000
#SBATCH --constraint=hsw

module load cp2k/2.7-dev
export PYTHONPATH=$PYTHONPATH:/homeappl/home/juritala/.local/lib/python2.7/site-packages/pycp2k-0.1-py2.7.egg/pycp2k
export REAL_ELPA_KERNEL=REAL_ELPA_KERNEL_AVX_BLOCK2

USAGE="Usage: $0 <result_db_file> <global_res_db_file>"

if [ "$#" != "2" ]; then
    echo "$USAGE"
    exit 1
fi

RESULT_DB_FILE=$1
GLOBAL_RES_DB_FILE=$2

ORIG_DIR=$PWD
cd $TMPDIR/$SLURM_JOB_ID

cp $ORIG_DIR/cp2k_system_params.py .
cp $CP2K_DATA/BASIS_MOLOPT .
cp $CP2K_DATA/GTH_POTENTIALS .

trap "mkdir $ORIG_DIR/$SLURM_JOB_ID; mv * $ORIG_DIR/$SLURM_JOB_ID/; exit" ERR TERM

python $KPFM_GLOBAL_SCRIPTS/calc_atomic_forces.py $ORIG_DIR $RESULT_DB_FILE $GLOBAL_RES_DB_FILE
