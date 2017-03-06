#!/bin/bash
#SBATCH -J kpfm_init
#SBATCH -p parallel
#SBATCH -n 144
#SBATCH --exclusive
#SBATCH --time=12:00:00 --mem-per-cpu=7000
#SBATCH --constraint=hsw

module load cp2k/2.7-dev
export PYTHONPATH=$PYTHONPATH:/homeappl/home/juritala/.local/lib/python2.7/site-packages/pycp2k-0.1-py2.7.egg/pycp2k
export REAL_ELPA_KERNEL=REAL_ELPA_KERNEL_AVX_BLOCK2

PROJECT_NAME=nacl_slab_cu_tip

USAGE="Usage: $0 <restart_slurm_id>"

if [ "$#" != "1" ]; then
    echo "$USAGE"
    exit 1
fi

RESTART_SLURM_ID=$1

ORIG_DIR=$PWD
cd $TMPDIR/$SLURM_JOB_ID

cp $ORIG_DIR/$RESTART_SLURM_ID/* .

mkdir $ORIG_DIR/$SLURM_JOB_ID
trap "mv * $ORIG_DIR/$SLURM_JOB_ID/; exit" TERM EXIT

srun cp2k.popt -i $PROJECT_NAME-1.restart -o $PROJECT_NAME.out
