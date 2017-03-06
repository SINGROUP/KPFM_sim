#!/bin/bash
#SBATCH -p parallel
#SBATCH -n 120
#SBATCH --exclusive
#SBATCH --time=1:00:00 --mem-per-cpu=4000
#SBATCH --constraint=hsw

module load cp2k/2.7-dev
export PYTHONPATH=$PYTHONPATH:/homeappl/home/juritala/.local/lib/python2.7/site-packages/pycp2k-0.1-py2.7.egg/pycp2k
export REAL_ELPA_KERNEL=REAL_ELPA_KERNEL_AVX_BLOCK2

USAGE="Usage: $0 <x_tip> <y_tip> <s_tip> <V> <result_db_file> <global_res_db_file>"

if [ "$#" != "6" ]; then
    echo "$USAGE"
    exit 1
fi

X=$1
Y=$2
S=$3
V=$4
RESULT_DB_FILE=$5
GLOBAL_RES_DB_FILE=$6

ORIG_DIR=$PWD
cd $TMPDIR/$SLURM_JOB_ID

cp $CP2K_DATA/BASIS_MOLOPT .
cp $CP2K_DATA/GTH_POTENTIALS .
cp $ORIG_DIR/cp2k_system_params.py .

mkdir $ORIG_DIR/$SLURM_JOB_ID
trap "mv * $ORIG_DIR/$SLURM_JOB_ID/; exit" ERR TERM

python $KPFM_GLOBAL_SCRIPTS/calc_potential.py $X $Y $S $V $ORIG_DIR $RESULT_DB_FILE $GLOBAL_RES_DB_FILE
mv * $ORIG_DIR/$SLURM_JOB_ID/
