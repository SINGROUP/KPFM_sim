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

PROJECT_NAME=kpfm_init_cu_tip_on_nacl

ORIG_DIR=$PWD
cd $TMPDIR/$SLURM_JOB_ID

cp $ORIG_DIR/$PROJECT_NAME.py .
cp $CP2K_DATA/BASIS_MOLOPT .
cp $CP2K_DATA/GTH_POTENTIALS .
cp $ORIG_DIR/*.xyz .
cp $ORIG_DIR/*.txt .

mkdir $ORIG_DIR/$SLURM_JOB_ID
trap "mv * $ORIG_DIR/$SLURM_JOB_ID/; exit" TERM EXIT

python $PROJECT_NAME.py
